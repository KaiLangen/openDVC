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
  _iRange = 8;

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


/*  if (_srcSize < _trgSize){
    // Src < Trg: up-size src then motion search 
    lowpassFilter(_currBuffer, srcBuffer, _srcWidth, _srcHeight, 3);
    bilinear(_prevBuffer, _prevKeyFrame, _srcWidth, _srcHeight,
             _srcWidth,_srcHeight, 0, 0);
    bilinear(srcBuffer, _currFrame, _srcWidth, _srcHeight,
             _srcWidth, _srcHeight, 0, 0);
    bilinear(_nextBuffer, _nextKeyFrame, _srcWidth, _srcHeight,
             _srcWidth,_srcHeight, 0, 0);
    ME(_prevKeyFrame, _nextKeyFrame, _currFrame);
    MC(prevTrg, nextTrg, currTrg);
  } else {
  */
    ME(_prevBuffer, _nextBuffer, _currBuffer);
    MC(prevTrg, nextTrg, currTrg);
 /*   if (_srcSize > _trgSize) {
      bilinear(prevTrg, _prevBuffer, _trgWidth, _trgHeight,
               _trgWidth,_trgHeight, 0, 0);
      bilinear(currTrg, _currBuffer, _trgWidth, _trgHeight,
               _trgWidth, _trgHeight, 0, 0);
      bilinear(nextTrg, _nextBuffer, _trgWidth, _trgHeight,
               _trgWidth, _trgHeight, 0, 0);
      MC(_prevBuffer, _nextBuffer, _currBuffer);
      bilinear(_currBuffer, currTrg, _trgWidth/2, _trgHeight/2,
               _srcWidth, _srcHeight, 0, 0);
    } else {
      MC(prevTrg, nextTrg, currTrg);
    }
  }
  memcpy(currTrg, prevTrg, _trgSize);
  */

/*  lowpassFilter(trgLowPass, currTrg, _trgWidth, _trgHeight, 5);
  FILE* fout = fopen("output", "wb");
  fwrite(prevTrg, _trgSize, 1, fout);
  for (int i = 0; i < _trgSize>>1; i++) fputc(127, fout);
  fwrite(currTrg, _trgSize, 1, fout);
  for (int i = 0; i < _trgSize>>1; i++) fputc(127, fout);
  fwrite(nextTrg, _trgSize, 1, fout);
  for (int i = 0; i < _trgSize>>1; i++) fputc(127, fout);
  fclose(fout);*/
  //delete [] srcBuffer;
}

void TwoStage::ME(imgpel* prevFrame, imgpel* nextFrame, imgpel* currFrame)
{
  mvinfo mv1, mv2;
  unsigned int cost1, cost2, index, pos;
  double L = floor(log2(_param + 1.0));
  double stepMax = pow(2.0, (L-1.0));
  int step = (int)stepMax;
  for (int y = 0; y < _srcHeight; y+=_iRange) {
    for (int x = 0; x < _srcWidth; x+=_iRange) {
      pos = x + y*_srcWidth;
      index = x/_iRange + y*_srcWidth/(_iRange * _iRange);
      _mvs[index].iCx = x;
      _mvs[index].iCy = y;

      cost1 = TSS(currFrame, prevFrame, mv1, 
                 step, pos, _srcWidth, _srcHeight, _iRange);
      cost2 = TSS(currFrame, nextFrame, mv2, 
                 step, pos, _srcWidth, _srcHeight, _iRange);
      //cout << cost1 << "," << cost2 << endl;
      if (cost1 < cost2) {
        _mvs[index].iMvx = mv1.iMvx;
        _mvs[index].iMvy = mv1.iMvy;
        _mvs[index].frameNo = 0;
      } else {
        _mvs[index].iMvx = mv2.iMvx;
        _mvs[index].iMvy = mv2.iMvy;
        _mvs[index].frameNo = 1;
      }
/*      cout << _mvs[index].iCx << " ," << 
              _mvs[index].iCy << " ," <<
              _mvs[index].iMvx << " ," <<
              _mvs[index].iMvy <<  " ," <<
              _mvs[index].frameNo << endl;
*/
    }
  }
}

void TwoStage::MC(imgpel* prev, imgpel* next, imgpel* curr)
{
  imgpel* ref;
  int cX, cY, mvX, mvY;
  int numMV = _srcSize / (_iRange * _iRange);
  for (int i = 0; i < numMV; i++) {
    // get the "start" values: coordinates of the top-left pixel of each MB
    // get the frame that the motion vector references, then increment it
    ref = (_mvs[i].frameNo == 1) ? next : prev;
    cX  = _mvs[i].iCx / 2;
    mvX = cX + _mvs[i].iMvx / 2;
    cY  = _mvs[i].iCy / 2;
    mvY = cY + _mvs[i].iMvy / 2;
    
    for (int j = 0; j < _iRange/2; j++) {
      // copy each row in the MB
      memcpy(curr + cX + (cY + j) * _srcWidth / 2,
             ref + mvX + (mvY + j) * _srcWidth / 2, 
             _iRange / 2);
    }
  }
}
