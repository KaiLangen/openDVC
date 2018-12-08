
#ifndef ENCODER_INC_FRAMEBUFFER_H
#define ENCODER_INC_FRAMEBUFFER_H

#include "config.h"

class FrameBuffer
{
public:
  FrameBuffer(int gop = 0)
  {
    _prevFrame        = new imgpel[FSIZE];
    _currFrame        = new imgpel[FSIZE];
    _nextFrame        = new imgpel[FSIZE];
    _sideInfoFrame    = new imgpel[FSIZE];

    if (gop != 0) {
      _recFrames      = new imgpel*[gop-1];

      for (int i = 0; i < gop-1; i++)
        _recFrames[i] = new imgpel[frameSize];
    }

    _dctFrame         = new int[FSIZE];
    _quantDctFrame    = new int[FSIZE];
    _decFrame         = new int[FSIZE];
    _invQuantDecFrame = new int[FSIZE];
  };

  imgpel*  getPrevFrame()        { return _prevFrame; };
  imgpel*  getCurrFrame()        { return _currFrame; };
  imgpel*  getNextFrame()        { return _nextFrame; };
  imgpel*  getSideInfoFrame()    { return _sideInfoFrame; };
  imgpel** getRecFrames()        { return _recFrames; };
  int*     getDctFrame()         { return _dctFrame; };
  int*     getQuantDctFrame()    { return _quantDctFrame; };
  int*     getDecFrame()         { return _decFrame; };
  int*     getInvQuantDecFrame() { return _invQuantDecFrame; };

private:
  imgpel*  _prevFrame;
  imgpel*  _currFrame;
  imgpel*  _nextFrame;
  imgpel*  _sideInfoFrame;
  imgpel** _recFrames;
  int*     _dctFrame;
  int*     _quantDctFrame;
  int*     _decFrame;
  int*     _invQuantDecFrame;
};

#endif // ENCODER_INC_FRAMEBUFFER_H

