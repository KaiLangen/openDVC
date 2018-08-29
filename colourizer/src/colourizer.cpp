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

  _gop = atoi(configMap["Gop"].c_str());
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
  FILE* fWritePtr   = _files->getFile("out")->getFileHandle();
  FILE* fWZReadPtr  = _files->getFile("wz")->getFileHandle();
  FILE* fKeyReadPtr = _files->getFile("key")->getFileHandle();
  imgpel* prevKeyFrame = new imgpel[_yuvFrameSize];
  imgpel* currFrame = new imgpel[_yuvFrameSize];
  // Read first key frame
  for (int keyFrameNo = 0; keyFrameNo < (_nframes-1)/_gop; keyFrameNo++) {
    // read out one key frame (this also moves the file ptr), then
    // write it to the output file
    fread(prevKeyFrame, _yuvFrameSize, 1, fKeyReadPtr);
    fwrite(prevKeyFrame, _yuvFrameSize, 1, fWritePtr);

    // skip grayscale version of Key-frame, then write out each of the
    // grayscale frames using the chroma channels of the previous key-frame
    fseek(fWZReadPtr, _grayFrameSize, SEEK_CUR);
    for (int i = 0; i < _gop-1; i++) {
      fread(currFrame, _grayFrameSize, 1, fWZReadPtr);
      addColour(prevKeyFrame, currFrame);
      fwrite(currFrame, _yuvFrameSize, 1, fWritePtr);
    }
  }
  // write out the last key frame
  fread(currFrame, _yuvFrameSize, 1, fKeyReadPtr);
  fwrite(currFrame, _yuvFrameSize, 1, fWritePtr);
  delete [] prevKeyFrame;
  delete [] currFrame;
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

