
#ifndef DECODER_INC_CAVLCDEC_H
#define DECODER_INC_CAVLCDEC_H

#include "config.h"
#include "cavlc.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class CavlcDec : public Cavlc
{
public:
  CavlcDec(Codec* codec, int blockSize);

  int decode(int* iDCT, int ix, int iy, int c);

  void clearNnz(int index, int c) { _mbs[c][index].nnz[0][0] = 0; };

private:
  int decodeNumCoeffTrailingOnes(int& numCoeff, int& t1s, int vlc, int c);
  int decodeLevel(int iNumCoef,int iTrailingOnes,int *iLevel,int *iOnes, int c);
  int decodeTotalZero(int &iTotalZeros,int iNumCoef, int c);
  int decodeRun(int &iRun,int iZerosLeft, int c);
};

#endif // DECODER_INC_CAVLCDEC_H

