#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>

#include "codec.h"
#include "ME.h"
#include "sideInformation.h"

using namespace std;
TwoStage::TwoStage(Codec* codec, FILE* file, int srcHeight, int srcWidth) 
: SideInformation(codec)
{
  // src dimensions
  _srcHeight = srcHeight;
  _srcWidth = srcWidth;
  _srcSize = _srcHeight * _srcWidth;

  // target dimensions
  _trgHeight = codec->getFrameHeight();
  _trgWidth = codec->getFrameWidth();
  _trgSize = _trgHeight * _trgWidth;

  // hard-coded for now
  _param = 7;
  _iRange = 16;

  _file = file;
  if (_trgSize != _srcSize) {
    _prevBuffer = new imgpel[_srcSize];
    _nextBuffer = new imgpel[_srcSize];
    _currBuffer = new imgpel[_srcSize];
  }

  _mvs = new mvinfo[_trgSize / (_iRange * _iRange)];
  _prevKeyFrame = new imgpel[_trgSize];
  _nextKeyFrame = new imgpel[_trgSize];
  _currFrame =  new imgpel[_trgSize];
}

TwoStage::~TwoStage()
{
  fclose(_file);

  if (_trgSize != _srcSize) {
    delete [] _prevBuffer;
    delete [] _currBuffer;
    delete [] _nextBuffer;
  }

  delete [] _mvs;
  delete [] _prevKeyFrame;
  delete [] _nextKeyFrame;
  delete [] _currFrame;
}

/* Needs to be changed from original version*/
void TwoStage::createSideInfo(imgpel* prevTrg, imgpel* nextTrg, imgpel* currTrg,
                              int prevFrameNo, int nextFrameNo, int currFrameNo)
{
  int index, pos, param, n;
//  imgpel* srcBuffer = new imgpel[_srcWidth*_srcHeight];
  /* find the correct frames */
  fseek(_file, _srcSize * prevFrameNo, SEEK_SET);
  fread(_prevBuffer, _srcSize, 1, _file);
  fseek(_file, _srcSize * currFrameNo, SEEK_SET);
  fread(_currBuffer, _srcSize, 1, _file);
  fseek(_file, _srcSize * nextFrameNo, SEEK_SET);
  fread(_nextBuffer, _srcSize, 1, _file);


  if (_srcSize < _trgSize){
    // Src < Trg: up-size src then motion search 
//    lowpassFilter(_currBuffer, srcBuffer, _srcWidth, _srcHeight, 3);
    bilinear(_prevBuffer, _prevKeyFrame, _srcWidth, _srcHeight,
             _srcWidth, _srcHeight, 0, 0);
    bilinear(_currBuffer, _currFrame, _srcWidth, _srcHeight,
             _srcWidth, _srcHeight, 0, 0);
    bilinear(_nextBuffer, _nextKeyFrame, _srcWidth, _srcHeight,
             _srcWidth, _srcHeight, 0, 0);
    ME(_prevKeyFrame, _nextKeyFrame, _currFrame, _trgWidth, _trgHeight);
    MC(prevTrg, nextTrg, currTrg, 1);
  } else {
    ME(_prevBuffer, _nextBuffer, _currBuffer, _srcWidth, _srcHeight);
    if (_srcSize > _trgSize)
      MC(prevTrg, nextTrg, currTrg, 2);
    else
      MC(prevTrg, nextTrg, currTrg, 1);
  }
}
/*
    FILE* fout = fopen("output", "wb");
    fwrite(prevTrg, _trgSize, 1, fout);
    for (int i = 0; i < _trgSize>>1; i++) fputc(127, fout);
    fwrite(currTrg, _trgSize, 1, fout);
    for (int i = 0; i < _trgSize>>1; i++) fputc(127, fout);
    fwrite(nextTrg, _trgSize, 1, fout);
    for (int i = 0; i < _trgSize>>1; i++) fputc(127, fout);
    fclose(fout);
*/

void TwoStage::ME(imgpel* prevFrame, imgpel* nextFrame, imgpel* currFrame,
                  int width, int height)
{
  mvinfo mv1, mv2;
  unsigned int cost1, cost2, index, pos;
/*  double L = floor(log2(_param + 1.0));
  double stepMax = pow(2.0, (L-1.0));
  int step = (int)stepMax;*/
  for (int y = 0; y < height; y+=_iRange) {
    for (int x = 0; x < width; x+=_iRange) {
      pos = x + y*width;
      index = x/_iRange + y*width/(_iRange * _iRange);
      _mvs[index].iCx = x;
      _mvs[index].iCy = y;

      cost1 = ES(currFrame, prevFrame, mv1, 
                 _param, pos, width, height, _iRange);
      cost2 = ES(currFrame, nextFrame, mv2, 
                 _param, pos, width, height, _iRange);
      if (cost1 < cost2) {
        _mvs[index].iMvx = mv1.iMvx;
        _mvs[index].iMvy = mv1.iMvy;
        _mvs[index].frameNo = 0;
      } else {
        _mvs[index].iMvx = mv2.iMvx;
        _mvs[index].iMvy = mv2.iMvy;
        _mvs[index].frameNo = 1;
      }
    }
  }
}

void TwoStage::MC(imgpel* prev, imgpel* next, imgpel* curr, int factor)
{
  imgpel* ref;
  int cX, cY, mvX, mvY, size, width, numMV;
  if (factor == 1) {
    size = _srcSize;
    width = _srcWidth;
  } else {
    size = _srcSize;
    width = _trgWidth;
  }
  numMV = size / (_iRange * _iRange);
  for (int i = 0; i < numMV; i++) {
    // get the "start" values: coordinates of the top-left pixel of each MB
    // get the frame that the motion vector references, then increment it
    ref = (_mvs[i].frameNo == 1) ? next : prev;
    cX  = _mvs[i].iCx / factor;
    mvX = cX + _mvs[i].iMvx / factor;
    cY  = _mvs[i].iCy / factor;
    mvY = cY + _mvs[i].iMvy / factor;
    
    for (int j = 0; j < _iRange / factor; j++) {
      // copy each row in the MB
      memcpy(curr + cX + (cY + j) * width,
             ref + mvX + (mvY + j) * width, 
             _iRange / factor);
    }
  }
}
