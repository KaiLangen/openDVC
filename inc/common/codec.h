
#ifndef COMMON_INC_CODEC_H
#define COMMON_INC_CODEC_H

#include "config.h"

class Bitstream;

class Codec
{
public:
  Codec() {};

  int getQp()                              { return _qp; };
  int getKeyQp()                           { return _keyQp; };
  int getNumChnCodeBands(int c)            { return _numChnCodeBands[c]; };


  double* getSigma()                       { return _sigma; };

  int getQuantMatrix(int qp, int x, int y) { return QuantMatrix[qp][y][x]; };

  int getQuantStep(int c, int x, int y)    { return _quantStep[c][y][x]; };
  double* getAlpha(int c)                  { return _alpha[c]; };
  Bitstream* getBitstream(int c)           { return _bs[c]; };
  double* getAverage(int c)                { return _average[c]; };

protected:
  const static int  ResidualBlockSize;
  const static int  SkipBlockSize;

  const static int  QuantMatrix[8][4][4];
  const static int  BitPlaneNum[8];
  const static int  MaxBitPlane[4][4];
# if RESIDUAL_CODING
  const static int  MinQStepSize[8][4][4];
# endif

  const static int  ScanOrder[16][2];
  const static int  HuffmanCodeValue[4][3][16];
  const static int  HuffmanCodeLength[4][3][16];


  int               _numFrames;
  int               _qp;
  int               _keyQp;
  int               _gopLevel;
  int               _gop;

  double*           _sigma;

  int               _numChnCodeBands[NCHANS];
  bool*             _parity[NCHANS];
  double*           _dParity[NCHANS]; // TODO temporary for decoder
  double*           _average[NCHANS];
  double*           _alpha[NCHANS];
  unsigned char*    _crc[NCHANS];
  Bitstream*        _bs[NCHANS];
  int               _quantStep[NCHANS][4][4];
};

#endif // COMMON_INC_CODEC_H

