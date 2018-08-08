
#include <cstdio>

#include "ldpcaEnc.h"
#include "codec.h"

LdpcaEnc::LdpcaEnc(const string& fileName, Codec* codec) : Ldpca(fileName, codec)
{
}

void LdpcaEnc::encode(int* source, bool* accumulatedSyndrome)
{
  int m = (_n / _totalNumInc) * _numInc[64];

  for (int k = 0; k < m; k++)
    accumulatedSyndrome[k] = 0;

  // source * H'
  for (int k = 0; k < _n; k++)
    for (int l = _jc[k]; l < _jc[k+1]; l++)
      accumulatedSyndrome[_ir[64][l]] ^= (source[k] & 0x01);

  // accumulate
  for (int k = 1; k < m; k++)
    accumulatedSyndrome[k] ^= accumulatedSyndrome[k-1];
}

