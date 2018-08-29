#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <exception>

#include "ME.h"
#include "colourizer.h"

using namespace std;

MvMultiFrame::MvMultiFrame(map<string, string> configMap)
  : Colourizer(configMap)
{
  _iRange = atoi(configMap["BlockSize"].c_str());
  _nCandidates = atoi(configMap["RefFrames"].c_str());
  if (_nCandidates <= 0 ) {
    throw invalid_argument("Invalid number of reference frames");
  }
  _nMV    = _width * _height / (_iRange * _iRange);
  _mvs    = new mvinfo[_nMV];
  _param  = atoi(configMap["Param"].c_str());
}

void MvMultiFrame::addColour(imgpel** keyFrameList, imgpel* currFrame)
{
  int index, pos, param, n;
  mvinfo tempMV;
  unsigned int cost = UINT_MAX;
  unsigned int minSAD;
  // take the nearest power of two to  param to find
  // the step size in pixels
  double L = floor(log2(_param + 1.0));
  double stepMax = pow(2.0, (L-1.0));
  int step = (int)stepMax;
  for (int y = 0; y < _height; y+=_iRange) {
    for (int x = 0; x < _width; x+=_iRange) {
      pos = x + y*_width;
      index = x/_iRange + y*_width/(_iRange * _iRange);
      _mvs[index].iCx = x;
      _mvs[index].iCy = y;
      minSAD = UINT_MAX;

      for (int i = 0; i < _nCandidates; i++) {
        cost = TSS(currFrame, keyFrameList[i], tempMV, 
                   step, pos, _width, _height, _iRange);
        if (cost < minSAD) {
          _mvs[index].iMvx = tempMV.iMvx;
          _mvs[index].iMvy = tempMV.iMvy;
          _mvs[index].frameNo = i;
          minSAD = cost;
        }
      }
    }
  }
  MC(keyFrameList, currFrame);
}

void MvMultiFrame::colourize()
{
  FILE* fWritePtr   = _files->getFile("out")->getFileHandle();
  FILE* fWZReadPtr  = _files->getFile("wz")->getFileHandle();
  FILE* fKeyReadPtr = _files->getFile("key")->getFileHandle();
  int index;

  imgpel* currFrame = new imgpel[_yuvFrameSize];
  // allocate and read in the first N keyframes
  imgpel** candidates = new imgpel*[_nCandidates];
  for (int i = 0; i < _nCandidates; i++) {
    candidates[i] = new imgpel[_yuvFrameSize];
    fread(candidates[i], 1, _yuvFrameSize, fKeyReadPtr);
  }

  for (int keyFrameNo = 0; keyFrameNo < (_nframes-1)/_gop; keyFrameNo++) {
    // skip grayscale version of Key-frame
    fseek(fWZReadPtr, _grayFrameSize, SEEK_CUR);
    // compute the index of the oldest key frame (candidates is a FIFO),
    // ithen write out that key frame to file 
    index = keyFrameNo % _nCandidates; 
    fwrite(candidates[index], _yuvFrameSize, 1, fWritePtr);
    for (int i = 0; i < _gop-1; i++) {
      fread(currFrame, _grayFrameSize, 1, fWZReadPtr);
      addColour(candidates, currFrame);
      //memset(currFrame + _grayFrameSize, 0x80, _grayFrameSize>>1);
      fwrite(currFrame, _yuvFrameSize, 1, fWritePtr);
    }

    // overwrite the oldest keyframe with a new one
    fread(candidates[index], _yuvFrameSize, 1, fKeyReadPtr);
  }
  // write out the last frame to file
  index = ((_nframes-1)/_gop) % _nCandidates; 
  fwrite(candidates[index], _yuvFrameSize, 1, fWritePtr);

  for (int i = 0; i < _nCandidates; i++)
    delete [] candidates[i];
  delete [] candidates;
  delete [] currFrame;
}


void MvMultiFrame::MC(imgpel** keyFrameList, imgpel* currFrame)
{
  int cX, cY, mvX, mvY, voffset;
  imgpel* imgRefUV;
  //= prevKeyFrame + _grayFrameSize;
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
    // get the frame that the motion vector references, then increment it
    // by the Luma channel size
    imgRefUV = keyFrameList[_mvs[i].frameNo] + _grayFrameSize;
    
    for (int j = 0; j < _iRange / 2; j++) {
      // copy each row in the MB
      // U channel
      memcpy(imgCurrUV + cX + (cY + j) * _width / 2,
             imgRefUV + mvX + (mvY + j) * _width / 2, 
             _iRange / 2);

      // V channel
      memcpy(imgCurrUV + voffset + cX + (cY + j) * _width / 2,
             imgRefUV + voffset + mvX + (mvY + j) * _width / 2, 
             _iRange / 2);
    }
  }
}
