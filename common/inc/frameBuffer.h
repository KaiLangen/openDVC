
#ifndef ENCODER_INC_FRAMEBUFFER_H
#define ENCODER_INC_FRAMEBUFFER_H

#include "config.h"

class FrameBuffer
{
public:
  FrameBuffer(int width, int height, int gop = 0)
  {
    int frameSize = width * height;
    _prevFrame        = new imgpel[frameSize];
    _currFrame        = new imgpel[frameSize];
    _nextFrame        = new imgpel[frameSize];
    _sideInfoFrame    = new imgpel[frameSize];

    if (gop != 0) {
      _recFrames      = new imgpel*[gop];

      for (int i = 0; i < gop; i++)
        _recFrames[i] = new imgpel[frameSize];
    }

    _dctFrame         = new int[frameSize];
    _quantDctFrame    = new int[frameSize];
    _decFrame         = new int[frameSize];
    _invQuantDecFrame = new int[frameSize];
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

