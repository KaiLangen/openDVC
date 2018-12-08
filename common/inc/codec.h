
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

  double* getAverage()                     { return _average; };
  double* getSigma()                       { return _sigma; };

  int getQuantMatrix(int qp, int x, int y) { return QuantMatrix[qp][y][x]; };
  int getQuantStep(int x, int y)           { return _quantStep[y][x]; };

  double* getAlpha(int i)                  { return _alpha[i]; };
  Bitstream* getBitstream(int i)           { return _bs[i]; };

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

  int               _quantStep[4][4];

  int               _numFrames;
  int               _qp;
  int               _keyQp;
  int               _gopLevel;
  int               _gop;

  double*           _dParity; // TODO temporary for decoder
  double*           _average;
  double*           _sigma;

  bool*             _parity[NCHANS];
  double*           _alpha[NCHANS];
  unsigned char*    _crc[NCHANS];
  Bitstream*        _bs[NCHANS];
};

#endif // COMMON_INC_CODEC_H

