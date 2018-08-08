
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

  int decode(int* iDCT, int ix, int iy);

  void clearNnz(int index) { _mbs[index].nnz[0][0] = 0; };

private:
  int decodeNumCoeffTrailingOnes(int& numCoeff, int& t1s, int vlc);
  int decodeLevel(int iNumCoef,int iTrailingOnes,int *iLevel,int *iOnes);
  int decodeTotalZero(int &iTotalZeros,int iNumCoef);
  int decodeRun(int &iRun,int iZerosLeft);
};

#endif // DECODER_INC_CAVLCDEC_H

