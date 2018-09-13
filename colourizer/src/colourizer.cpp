#include <sstream>
#include <fstream>
#include <iostream>
#include <exception>

#include "colourizer.h"

using namespace std;

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
  _files->addFile("wz", configMap["WZFile"])->openFile("rb");
  if (!_files->getFile("wz")->getFileHandle()) {
    cerr << "No such file: " << configMap["WZFile"] << endl;
    throw invalid_argument("Invalid WZ file");
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
  int frameStep, idx, wzFrameNo, prevIdx, nextIdx;
  FILE* fWritePtr   = _files->getFile("out")->getFileHandle();
  FILE* fWZReadPtr  = _files->getFile("wz")->getFileHandle();
  FILE* fKeyReadPtr = _files->getFile("key")->getFileHandle();
  imgpel* tmp_ptr, *prevFrame, *currFrame, *nextFrame;
  imgpel* prevKeyFrame = new imgpel[_yuvFrameSize];
  imgpel* nextKeyFrame = new imgpel[_yuvFrameSize];
  imgpel** recFrames = new imgpel*[_gop];

  for (int i = 0; i < _gop-1; i++)
    recFrames[i] = new imgpel[_yuvFrameSize];

  // Read first key frame
  fread(prevKeyFrame, _yuvFrameSize, 1, fKeyReadPtr);
  for (int keyFrameNo = 0; keyFrameNo < (_nframes-1)/_gop; keyFrameNo++) {

    // read in the next key frame (this also moves the file ptr)
    fread(nextKeyFrame, _yuvFrameSize, 1, fKeyReadPtr);

    // skip grayscale version of Key-frame, then write out each of the
    // grayscale frames using the chroma channels of the previous key-frame
    fseek(fWZReadPtr, _grayFrameSize, SEEK_CUR);


    // start with frameStep = 1/2 gop, then half it each iteration
    for (int il = 0; il < _gopLevel; il++) {
      frameStep = _gop / ((il+1)<<1);
      idx = frameStep;

      // Start decoding the WZ frame
      while (idx < _gop) {
        wzFrameNo = keyFrameNo*_gop + idx;
				cout << "Colouring frame " << wzFrameNo << " (Wyner-Ziv frame)" << endl;

        // Setup frame pointers within the GOP
        prevIdx = idx - frameStep;
        nextIdx = idx + frameStep;

        currFrame = recFrames[idx-1];
        prevFrame = (prevIdx == 0)    ? prevKeyFrame :
                                        recFrames[prevIdx-1];
        nextFrame = (nextIdx == _gop) ? nextKeyFrame :
                                        recFrames[nextIdx-1];
				
				fseek(fWZReadPtr, wzFrameNo*_grayFrameSize, SEEK_SET);
				fread(currFrame, _grayFrameSize, 1, fWZReadPtr);
				addColour(prevFrame, nextFrame, currFrame);

        idx += 2*frameStep;
      }
    }

    // ---------------------------------------------------------------------
    // Output decoded frames of the whole GOP
    // ---------------------------------------------------------------------

    // write out prevFrame to the output file
    fwrite(prevKeyFrame, _yuvFrameSize, 1, fWritePtr);

    for (int i = 0; i < _gop-1; i++)
      fwrite(recFrames[i], _yuvFrameSize, 1, fWritePtr);

    // swap prev and next key frame pointers;
    // nextKeyFrame will be overwritten in the next step
    tmp_ptr = prevKeyFrame;
    prevKeyFrame = nextKeyFrame;
    nextKeyFrame = tmp_ptr;
  }
  // write out the last key frame
  fread(currFrame, _yuvFrameSize, 1, fKeyReadPtr);
  fwrite(currFrame, _yuvFrameSize, 1, fWritePtr);
  delete [] prevKeyFrame;
  delete [] nextKeyFrame;
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
void Colourizer::addColour(imgpel* prevKeyFrame, imgpel* nextKeyFrame, imgpel* currFrame)
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

