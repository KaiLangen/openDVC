
#ifndef COMMON_INC_CAVLC_H
#define COMMON_INC_CAVLC_H

#include "config.h"

class Codec;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class Cavlc
{
public:
  Cavlc(Codec* codec, int blockSize);

protected:
  const static int ScanOrder[16][2];
  const static byte NumVlcTableL[3][4][17];
  const static byte NumVlcTableC[3][4][17];
  const static byte TotalZerosTableL[15][16];
  const static byte TotalZerosTableC[15][16];
  const static byte RunTableL[15][16];
  const static byte RunTableC[15][16];

  int getNumNonzero(int x,int y);

  Codec*  _codec;

  int     _blockSize;

  mb*     _mbs;
};

#endif // COMMON_INC_CAVLC_H

