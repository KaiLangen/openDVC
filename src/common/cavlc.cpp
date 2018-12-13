
#include "cavlc.h"
#include "codec.h"
#include "defs.h"

const int Cavlc::ScanOrder[16][2] = {
  {0, 0}, {1, 0}, {0, 1}, {0, 2},
  {1, 1}, {2, 0}, {3, 0}, {2, 1},
  {1, 2}, {0, 3}, {1, 3}, {2, 2},
  {3, 1}, {3, 2}, {2, 3}, {3, 3}
};

const byte Cavlc::NumVlcTableL[3][4][17] =
{
  {
    { 1, 6, 8, 9,10,11,13,13,13,14,14,15,15,16,16,16,16},
    { 0, 2, 6, 8, 9,10,11,13,13,14,14,15,15,15,16,16,16},
    { 0, 0, 3, 7, 8, 9,10,11,13,13,14,14,15,15,16,16,16},
    { 0, 0, 0, 5, 6, 7, 8, 9,10,11,13,14,14,15,15,16,16},
  },
  {
    { 2, 6, 6, 7, 8, 8, 9,11,11,12,12,12,13,13,13,14,14},
    { 0, 2, 5, 6, 6, 7, 8, 9,11,11,12,12,13,13,14,14,14},
    { 0, 0, 3, 6, 6, 7, 8, 9,11,11,12,12,13,13,13,14,14},
    { 0, 0, 0, 4, 4, 5, 6, 6, 7, 9,11,11,12,13,13,13,14},
  },
  {
    { 4, 6, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 9,10,10,10,10},
    { 0, 4, 5, 5, 5, 5, 6, 6, 7, 8, 8, 9, 9, 9,10,10,10},
    { 0, 0, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,10},
    { 0, 0, 0, 4, 4, 4, 4, 4, 5, 6, 7, 8, 8, 9,10,10,10},
  },
};

const byte Cavlc::NumVlcTableC[3][4][17] =
{
  {
    { 1, 5, 7, 7, 7, 7,15,11, 8,15,11,15,11,15,11, 7, 4},
    { 0, 1, 4, 6, 6, 6, 6,14,10,14,10,14,10, 1,14,10, 6},
    { 0, 0, 1, 5, 5, 5, 5, 5,13, 9,13, 9,13, 9,13, 9, 5},
    { 0, 0, 0, 3, 3, 4, 4, 4, 4, 4,12,12, 8,12, 8,12, 8},
  },
  {
    { 3,11, 7, 7, 7, 4, 7,15,11,15,11, 8,15,11, 7, 9, 7},
    { 0, 2, 7,10, 6, 6, 6, 6,14,10,14,10,14,10,11, 8, 6},
    { 0, 0, 3, 9, 5, 5, 5, 5,13, 9,13, 9,13, 9, 6,10, 5},
    { 0, 0, 0, 5, 4, 6, 8, 4, 4, 4,12, 8,12,12, 8, 1, 4},
  },
  {
    {15,15,11, 8,15,11, 9, 8,15,11,15,11, 8,13, 9, 5, 1},
    { 0,14,15,12,10, 8,14,10,14,14,10,14,10, 7,12, 8, 4},
    { 0, 0,13,14,11, 9,13, 9,13,10,13, 9,13, 9,11, 7, 3},
    { 0, 0, 0,12,11,10, 9, 8,13,12,12,12, 8,12,10, 6, 2},
  },
};

const byte Cavlc::TotalZerosTableL[15][16] =
{
  { 1, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9},
  { 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6},
  { 4, 3, 3, 3, 4, 4, 3, 3, 4, 5, 5, 6, 5, 6},
  { 5, 3, 4, 4, 3, 3, 3, 4, 3, 4, 5, 5, 5},
  { 4, 4, 4, 3, 3, 3, 3, 3, 4, 5, 4, 5},
  { 6, 5, 3, 3, 3, 3, 3, 3, 4, 3, 6},
  { 6, 5, 3, 3, 3, 2, 3, 4, 3, 6},
  { 6, 4, 5, 3, 2, 2, 3, 3, 6},
  { 6, 6, 4, 2, 2, 3, 2, 5},
  { 5, 5, 3, 2, 2, 2, 4},
  { 4, 4, 3, 3, 1, 3},
  { 4, 4, 2, 1, 3},
  { 3, 3, 1, 2},
  { 2, 2, 1},
  { 1, 1},
};

const byte Cavlc::TotalZerosTableC[15][16] =
{
  { 1, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 1},
  { 7, 6, 5, 4, 3, 5, 4, 3, 2, 3, 2, 3, 2, 1, 0},
  { 5, 7, 6, 5, 4, 3, 4, 3, 2, 3, 2, 1, 1, 0},
  { 3, 7, 5, 4, 6, 5, 4, 3, 3, 2, 2, 1, 0},
  { 5, 4, 3, 7, 6, 5, 4, 3, 2, 1, 1, 0},
  { 1, 1, 7, 6, 5, 4, 3, 2, 1, 1, 0},
  { 1, 1, 5, 4, 3, 3, 2, 1, 1, 0},
  { 1, 1, 1, 3, 3, 2, 2, 1, 0},
  { 1, 0, 1, 3, 2, 1, 1, 1, },
  { 1, 0, 1, 3, 2, 1, 1, },
  { 0, 1, 1, 2, 1, 3},
  { 0, 1, 1, 1, 1},
  { 0, 1, 1, 1},
  { 0, 1, 1},
  { 0, 1},
};

const byte Cavlc::RunTableL[15][16] =
{
  { 1, 1},
  { 1, 2, 2},
  { 2, 2, 2, 2},
  { 2, 2, 2, 3, 3},
  { 2, 2, 3, 3, 3, 3},
  { 2, 3, 3, 3, 3, 3, 3},
  { 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11},
  { 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11},
  { 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11},
  { 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11},
  { 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11},
  { 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11},
  { 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11},
  { 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11},
  { 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11}
};

const byte Cavlc::RunTableC[15][16] =
{
  { 1, 0},
  { 1, 1, 0},
  { 3, 2, 1, 0},
  { 3, 2, 1, 1, 0},
  { 3, 2, 3, 2, 1, 0},
  { 3, 0, 1, 3, 2, 5, 4},
  { 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  { 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  { 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  { 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  { 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  { 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  { 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  { 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  { 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1}
};

Cavlc::Cavlc(Codec* codec, int blockSize)
{
  _codec = codec;

  _blockSize = blockSize;

  int bplen;
  for (int c = 0; c < NCHANS; c++) {
    bplen = (c == 0) ? Y_BPLEN : UV_BPLEN;
    _mbs[c] = new mb[bplen];
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int Cavlc::getNumNonzero(int x, int y, int c)
{
  int frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;

  int mbX = x / _blockSize;
  int mbY = y / _blockSize;

  return _mbs[c][mbX + mbY*(frameWidth/_blockSize)].nnz[0][0];
}

