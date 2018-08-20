#include <sstream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "colourizer.h"
#include "fileManager.h"

using namespace std;

Colourizer::Colourizer(char **argv)
{
  _files = FileManager::getManager();
  map<string, string> configMap = readConfig(argv[1]);
  initialize(configMap);
}

map<string, string>&
Colourizer::readConfig(string filename)
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
    cerr << "Invalid config file!" << endl;
    return configMap;
  }
}

void
Colourizer::initialize(map<string, string>& configMap)
{
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
  _method = atoi(configMap["Method"].c_str());
  _nframes = atoi(configMap["FramesToBeEncoded"].c_str());
  string wz = configMap["WZFile"];
  _files->addFile("wz", wz)->openFile("rb");
  _files->addFile("key", configMap["KeyFile"])->openFile("rb");
  _files->addFile("out", "recoloured.yuv")->openFile("wb");


  _grayFrameSize = _width * _height;
  _yuvFrameSize = 3*(_grayFrameSize)>>1;
}

void
Colourizer::colourize()
{
  switch(_method)
  {
    case SIMPLE:
      cout << "Simple" << endl;
      simpleRecolour();
      break;
    case MVSEARCH:
      cout << "MVSearch" << endl;
      break;
    case ML:
      cout << "ML" << endl;
      break;
    default:
      cout << "Invalid colourizing method" << endl;
      break;
  }
}

/**
 * Description:
 * Recolour the WZ frames using the exact chroma channels from the Key frames.
 * Not REALLY a viable solution, but useful as a lower bounds for the PSNR of
 * a recoloured WZ video.
 * ----------------------------------------------------------------------------
 * Param: None
 * Return: None
 */
void
Colourizer::simpleRecolour()
{
  FILE* fWritePtr   = _files->getFile("out")->getFileHandle();
  FILE* fWZReadPtr   = _files->getFile("wz")->getFileHandle();
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
    memcpy(currFrame + _grayFrameSize,
           prevKeyFrame + _grayFrameSize,
           _grayFrameSize>>1);
    for (int i = 0; i < _gop-1; i++) {
      fread(currFrame, _grayFrameSize, 1, fWZReadPtr);
      fwrite(currFrame, _yuvFrameSize, 1, fWritePtr);
    }
  }
  // write out the last key frame
  fread(currFrame, _yuvFrameSize, 1, fKeyReadPtr);
  fwrite(currFrame, _yuvFrameSize, 1, fWritePtr);
  delete [] prevKeyFrame;
  delete [] currFrame;
}
