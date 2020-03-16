#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <exception>

#include "colourizer.h"

using namespace std;

double calcPSNR(unsigned char* img1,unsigned char* img2,int length)
{
  float PSNR;
  float MSE=0;

  for(int i=0;i<length;i++)
    {
      MSE+=pow(float(img1[i]-img2[i]),float(2.0))/length;
    }
  PSNR=10*log10(255*255/MSE);
  return PSNR;
}

Colourizer::Colourizer(map<string, string> configMap)
{
  _files = FileManager::getManager();
  // find sequence type in config map
  auto it = configMap.find("SequenceType");
  if (it != configMap.end() && !strcmp(it->second.c_str(), "CIF")) {
   _width = 352;
   _height = 288;
  } else if (it != configMap.end() && !strcmp(it->second.c_str(), "QCIF")) {
   _width = 176;
   _height = 144;
  } else {
    cerr << "Invalid sequence type!" << endl;
  }

  _gopLevel = atoi(configMap["GopLevel"].c_str());
  _gop = 1 << _gopLevel;
  _nframes = atoi(configMap["FramesToBeEncoded"].c_str());

  _files->addFile("original", configMap["InputFile"])->openFile("rb");
  if (!_files->getFile("original")->getFileHandle()) {
    cerr << "No such file: " << configMap["InputFile"] << endl;
    throw invalid_argument("Invalid input video file");
  }

  _files->addFile("key", configMap["KeyFile"])->openFile("rb");
  if (!_files->getFile("key")->getFileHandle()) {
    cerr << "No such file: " << configMap["KeyFile"] << endl;
    throw invalid_argument("Invalid Keyframe file");
  }

  _files->addFile("out", "recoloured.yuv")->openFile("wb");
  _grayFrameSize = _width * _height;
  _yuvFrameSize = 3*(_grayFrameSize)>>1;
}


void Colourizer::colourize()
{
  int frameStep, idx, frameNo, prevIdx;
  double dPSNRY = 0;
  double dPSNRU = 0;
  double dPSNRV = 0;
  double dPSNRAvg = 0;
  FILE* fWritePtr   = _files->getFile("out")->getFileHandle();
  FILE* fKeyReadPtr = _files->getFile("key")->getFileHandle();
  FILE* fOrigReadPtr = _files->getFile("original")->getFileHandle();
  imgpel* tmp_ptr, *prevFrame, *currFrame;
  imgpel* prevKeyFrame = new imgpel[_yuvFrameSize];
  imgpel* origFrame = new imgpel[_yuvFrameSize];
  imgpel** recFrames = new imgpel*[_gop];

  for (int i = 0; i < _gop-1; i++)
    recFrames[i] = new imgpel[_yuvFrameSize];

  for (int keyFrameNo = 0; keyFrameNo < (_nframes-1)/_gop; keyFrameNo++) {
    // skip grayscale version of Key-frame, then write out each of the
    // grayscale frames using the chroma channels of the previous key-frame
    fseek(fKeyReadPtr, keyFrameNo * _gop * _yuvFrameSize, SEEK_SET);
		fread(origFrame, _yuvFrameSize, 1, fOrigReadPtr);
    fread(prevKeyFrame, _yuvFrameSize, 1, fKeyReadPtr);
    dPSNRY = calcPSNR(origFrame, prevKeyFrame, _grayFrameSize);
    dPSNRU = calcPSNR(origFrame+_grayFrameSize, prevKeyFrame+_grayFrameSize, _grayFrameSize>>2);
    dPSNRV = calcPSNR(origFrame+5*(_grayFrameSize>>2),
                      prevKeyFrame+5*(_grayFrameSize>>2),
                      _grayFrameSize>>2);
    fwrite(prevKeyFrame, _yuvFrameSize, 1, fWritePtr);

    for (int idx = 0; idx < _gop-1; idx++) {
      frameNo = keyFrameNo * _gop + idx + 1;
			cout << "Decoding frame " << frameNo << " (Wyner-Ziv frame)" << endl;

      // Setup frame pointers within the GOP
      prevIdx = idx - 1;
      currFrame = recFrames[idx];
      prevFrame = (prevIdx < 0)    ? prevKeyFrame :
                                     recFrames[prevIdx];
      prevFrame = prevKeyFrame;                         
	  	fseek(fKeyReadPtr, frameNo*_yuvFrameSize, SEEK_SET);
			fread(currFrame, _grayFrameSize, 1, fKeyReadPtr);
			fread(origFrame, _yuvFrameSize, 1, fOrigReadPtr);
			addColour(prevFrame, currFrame);
      dPSNRY = calcPSNR(origFrame, currFrame, _grayFrameSize);
      dPSNRU = calcPSNR(origFrame+_grayFrameSize, currFrame+_grayFrameSize, _grayFrameSize>>2);
      dPSNRV = calcPSNR(origFrame+5*(_grayFrameSize>>2),
                        currFrame+5*(_grayFrameSize>>2),
                        _grayFrameSize>>2);
      dPSNRAvg = (6*dPSNRY + dPSNRU + dPSNRV) / 8;
      cout << "PSNR Recoloured Chroma (U): " << dPSNRU << endl;
      cout << "PSNR Recoloured Chroma (V): " << dPSNRV << endl;
      cout << "PSNR Luma: " << dPSNRY << endl;
      cout << "PSNR Frame Average: " << dPSNRAvg << endl;
      fwrite(currFrame, _yuvFrameSize, 1, fWritePtr);
    }
  }
  delete [] prevKeyFrame;
  for (int i = 0; i < _gop-1; i++)
    delete [] recFrames[i];
  delete [] recFrames;
}

void Colourizer::MC(imgpel* prevKeyFrame, imgpel* currFrame)
{
  int cX, cY, mvX, mvY, voffset;
  imgpel* imgPrevUV = prevKeyFrame + _grayFrameSize;
  imgpel* imgCurrUV = currFrame + _grayFrameSize;
  voffset = _width * _height / 4;
  // copy UV pixels from the previous frame given some motion vectors
  for (int i = 0; i < (_nMV); i++) {
    // get the "start" values: coordinates of the top-left pixel of each MB
    // one for both U and V channels
    cX  = _mvs[i].iCx / 2;
    mvX = cX + _mvs[i].iMvx / 2;
    cY  = _mvs[i].iCy / 2;
    mvY = cY + _mvs[i].iMvy / 2;
    
    for (int j = 0; j < _iRange / 2; j++) {
      // copy each row in the MB
      // U channel
      memcpy(imgCurrUV + cX + (cY + j) * _width / 2,
             imgPrevUV + mvX + (mvY + j) * _width / 2, 
             _iRange / 2);

      // V channel
      memcpy(imgCurrUV + voffset + cX + (cY + j) * _width / 2,
             imgPrevUV + voffset + mvX + (mvY + j) * _width / 2, 
             _iRange / 2);
    }
  }
}

/**
 * Description:
 * Recolour the WZ frames using the exact chroma channels from the Key frames.
 * Not REALLY a viable solution, but useful as a lower bounds for the PSNR of
 * a recoloured WZ video.
 * ----------------------------------------------------------------------------
 * Param: prevKeyFrame
 * Return: currFrame
 */
void Colourizer::addColour(imgpel* prevKeyFrame, imgpel* currFrame)
{
  memcpy(currFrame + _grayFrameSize,
         prevKeyFrame + _grayFrameSize,
         _grayFrameSize>>1);
}

map<string, string>&
readConfig(string filename)
{
  string line;
  ifstream cfile(filename);
  static map<string, string> configMap;
  if (cfile)
  {
    while( getline(cfile, line) )
    {
      string key, value;
      istringstream is_line(line);

      // Check if the first non-whitespace is a #
      if( getline(is_line >> ws, key, '=')  && !key.empty() && key[0] != '#')
        if( getline(is_line >> ws, value) )
          configMap[key] = value;
    }
    cfile.close();
    return configMap;
  }
  else
  {
    cerr << "No such file: " << filename << endl;
    throw invalid_argument("Invalid config file");
  }
}
