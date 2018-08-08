
#include <iostream>
#include <sstream>

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "cavlcEnc.h"
#include "codec.h"
#include "encoder.h"
#include "bitstream.h"
#include "fileManager.h"

using namespace std;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
CavlcEnc::CavlcEnc(Codec* codec, int blockSize) : Cavlc(codec, blockSize)
{
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcEnc::encode(int* frame, int* skipMask)
{
  int width    = _codec->getFrameWidth();
  int height   = _codec->getFrameHeight();
  int bitCount = 0;

  _patternFile = FileManager::getManager()->addFile("pattern_cavlc", "pattern_cavlc.dat");
  _patternFile->openFile("w");
  _patternFh = _patternFile->getFileHandle();

  for (int j = 0; j < height/_blockSize; j++)
    for (int i = 0; i < width/_blockSize; i++) {
# if SKIP_MODE
      if (skipMask[i+j*(width/_blockSize)] == 0) {
# endif // SKIP_MODE
        setupMacroBlock(frame, i, j);

        bitCount += encodeMacroBlock(i, j);
# if SKIP_MODE
      }
      else
        _mbs[i+j*(width/_blockSize)].nnz[0][0] = 0;
# endif // SKIP_MODE
    }

  _patternFile->closeFile();

  return bitCount;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void CavlcEnc::setupMacroBlock(int* frame, int mbX, int mbY)
{
  int   width  = _codec->getFrameWidth();
  int** buffer = AllocArray2D<int>(4, 4);

  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 4; i++) {
      int index = (i+mbX*_blockSize) + (j+mbY*_blockSize)*width;
      int qp = _codec->getQp();
      int mask;
      int sign;
      int value;

      mask  = (0x1<<(_codec->getQuantMatrix(qp, i, j)-1))-1;
      sign  = (frame[index]>>(_codec->getQuantMatrix(qp, i, j)-1)) & 0x1;
      value = frame[index] & mask;

      buffer[j][i] = (sign == 1) ? -value : value;

# if !RESIDUAL_CODING
      if (i == 0 && j == 0)
        buffer[j][i] = frame[index];
# endif // !RESIDUAL_CODING
    }
  }

  int nnz = 0;

# if MODE_DECISION
  for (int i = 0; i < (16-_codec->getNumChnCodeBands()); i++) {
    int x = ScanOrder[i+_codec->getNumChnCodeBands()][0];
    int y = ScanOrder[i+_codec->getNumChnCodeBands()][1];
# else // if !MODE_DECISION
  for (int i = 0; i < 16; i++) {
    int x = ScanOrder[i][0];
    int y = ScanOrder[i][1];
# endif // MODE_DECISION

    _mbs[mbX + mbY*(width/_blockSize)].coef_lac[0][0][i] = buffer[y][x];

    if (buffer[y][x] != 0)
      nnz++;
  }

  _mbs[mbX + mbY*(width/_blockSize)].nnz[0][0] = nnz;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcEnc::encodeMacroBlock(int mbX, int mbY)
{
  int         iWidth = _codec->getFrameWidth();
  int*        pLevel = 0;
  int         max_coeff_num;    // Number of bands coded by CAVLC
  int         x, y;
  int         Coef = 0;         // Total Coefficient
  int         trailone = 0;     // Trailing Ones num
  int         totalzeros = 0;   // total zeros
  int         totalbits = 0;    // bits count
  int         trailflag = 0;
  int         zeroflag = 0;

  vector<int> level;            // level
  vector<int> sign;             // sign bit for trailing ones
  vector<int> runs;             // runs
  vector<int> zeroleft;         // zeroleft

# if MODE_DECISION
  max_coeff_num = 16 - _codec->getNumChnCodeBands();
# else // if !MODE_DECISION
  max_coeff_num = 16;
# endif // MODE_DECISION

  pLevel = _mbs[mbX + mbY*(iWidth/_blockSize)].coef_lac[0][0];

  for (int i = max_coeff_num-1; i >= 0; i--) {
    if (pLevel[i] != 0) {
      Coef++;
      zeroflag = 1;

      if (abs(pLevel[i]) == 1 && trailone < 3 && trailflag == 0) {
        trailone++;

        if (pLevel[i] == 1)
          sign.push_back(0);
        else
          sign.push_back(1);
      }
      else {
        trailflag = 1;
        level.push_back(pLevel[i]);
      }

      runs.push_back(0);
    }
    else {
      if (zeroflag == 1) {
        totalzeros++;
        runs.back()++;
      }
    }
  }

  int totzero = totalzeros;

  for (unsigned i = 0; i < runs.size(); i++) {
    zeroleft.push_back(totzero);
    totzero = totzero-runs[i];
  }

  int nc = 0;
  int vlc = 0;

  x = mbX*4;
  y = mbY*4;

  if (x == 0 && y == 0)
    nc = 0;
  else if (x == 0 && y != 0)
    nc = getNumNonzero(x, y-4);
  else if (x != 0 && y == 0)
    nc = getNumNonzero(x-4, y);
  else
    nc = (getNumNonzero(x, y-4) + getNumNonzero(x-4, y) + 1) >> 1;

  if (nc < 2)
    vlc = 0;
  else if (nc < 4)
    vlc = 1;
  else if (nc < 8)
    vlc = 2;
  else
    vlc = 3;

  // encode totalCoef and trailingOnes
  SyntaxElement* tot = new SyntaxElement(Coef, trailone, vlc, 0);

  totalbits += encodeNumTrail(tot);

  delete tot;

  // trailingOnes sign bit
  if (trailone != 0)
    totalbits += encodeSignTrail(sign);

  // encode level
  // here assume using VLC0 initially
  SyntaxElement* lev = new SyntaxElement[Coef-trailone];

  int level_two_or_higher = (Coef > 3 && trailone == 3) ? 0 : 1;

  for (int index = 0; index < (Coef-trailone); index++) {
    lev[index].setValue1(level[index]);

    if (level_two_or_higher) {
      level_two_or_higher = 0;

      if (level[index] < 0)
        lev[index].setValue1(level[index]+1);
      else
        lev[index].setValue1(level[index]-1);
    }
  }

  int vlcnum = 0;

  if (Coef > 10 && trailone < 3)
    vlcnum = 1;

  for (int index = 0; index < (Coef-trailone); index++) {
    if (vlcnum == 0)
      totalbits += encodeLevelsVlc0((lev+index));
    else
      totalbits += encodeLevelsVlcN((lev+index),vlcnum);

    if (vlcnum == 0)
      vlcnum = 1;

    if (abs(lev[index].getValue1()) > int(3*pow(2.0, double(vlcnum-1))))
      ++vlcnum;

    if ((index == 0) && abs(level[index]) > 3)
      vlcnum = 2;
  }

  delete [] lev;

  //encode totalZeros
  if (Coef != 0 && Coef != max_coeff_num) {
    SyntaxElement* totzeros = new SyntaxElement(totalzeros, Coef, 0, 0);

    totalbits += encodeTotalZeros(totzeros);

    delete totzeros;
  }

  //encode runbefore
  SyntaxElement* run = new SyntaxElement[Coef];

  for (int index = 0; index < Coef-1; index++) {
    if (zeroleft[index] == 0)
      break;

    run[index].set(runs[index], Min(zeroleft[index]-1, 6), 0, 0);

    totalbits += encodeRuns(run+index);
  }

  delete [] run;

  return totalbits;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcEnc::symbol2vlc(SyntaxElement* sym)
{
  int info_len = sym->getLength();

  // vlc coding
  stringstream str;

  while (--info_len >= 0)
    str << (char)((0x01 & (sym->getInfo() >> info_len))+'0');

  _codec->getBitstream()->write(sym->getInfo(), sym->getLength());

# if TESTPATTERN
  for (int idx = 1; idx <= 8; idx++) {
    int value = ((sym->getInfo()) >> (32-idx*4)) & 0xf;
    fprintf(_patternFh, "%x", value);
  }
  fprintf(_patternFh, "\n");
# endif

  return info_len;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcEnc::encodeNumTrail(SyntaxElement* se)
{
  int numCoeff = se->getValue1();
  int t1s      = se->getValue2();
  int tableSel = se->getLength();

  if (tableSel != 3) {
    // Select one of the three Num-VLC tables
    se->setLength(NumVlcTableL[tableSel][t1s][numCoeff]);
    se->setInfo  (NumVlcTableC[tableSel][t1s][numCoeff]);
  }
  else {
    // Use 6-bit fixed length coding (FLC) xxxxyy
    // xxxx for NumCoeff-1 and yy for NumT1s
    // Codeword 000011 is used when NumCoeff = 0
    se->setLength(6);

    if (numCoeff > 0)
      se->setInfo(((numCoeff-1) << 2) | t1s);
    else
      se->setInfo(3);
  }

  symbol2vlc(se);

  return se->getLength();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcEnc::encodeSignTrail(vector<int>& sign)
{
  for (unsigned i = 0; i < sign.size(); i++) {
    SyntaxElement* se = new SyntaxElement(0, 0, 1, sign[i]);

    symbol2vlc(se);

    delete se;
  }

  return sign.size();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcEnc::encodeLevelsVlc0(SyntaxElement* se)
{
  int level  = se->getValue1();
  int sign   = (level < 0 ? 1 : 0);
  int levabs = abs(level);

  if (levabs < 8) {
    se->setLength(2*levabs - 1 + sign);
    se->setInfo(1);
  }
  else if (levabs < 16) {
    // Escape code 1: 000000000000001xxxx
    se->setLength(19);
    se->setInfo(16 | ((levabs << 1) - 16) | sign);
  }
  else {
    // Escape code 2: 0000000000000001xxxxxxxxxxxx
    int escapeBase   = 4096;
    int escapeOffset = levabs + 2048 - 16;
    int numPrefix    = 0;

    if (escapeOffset >= 4096) {
      numPrefix++;

      while (escapeOffset >= (4096 << numPrefix))
        numPrefix++;
    }

    escapeBase <<= numPrefix;

    se->setLength(28 + (numPrefix << 1));
    se->setInfo(escapeBase | ((escapeOffset << 1) - escapeBase) | sign);
  }

  symbol2vlc(se);

  return se->getLength();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcEnc::encodeLevelsVlcN(SyntaxElement* se, int vlc)
{
  int level  = se->getValue1();
  int sign   = (level < 0 ? 1 : 0);
  int levabs = abs(level) - 1;

  int shift  = vlc - 1;
  int escape = (15 << shift);

  if (levabs < escape) {
    int sufmask = ~((0xffffffff) << shift);
    int suffix  = (levabs) & sufmask;

    se->setLength(((levabs) >> shift) + 1 + vlc);
    se->setInfo((2 << shift) | (suffix << 1) | sign);
  }
  else {
    // Escape code: 0000000000000001xxxxxxxxxxxx
    int escapeBase = 4096;
    int escapeOffset = levabs + 2048 - escape;
    int numPrefix = 0;

    if ((escapeOffset) >= 4096) {
      numPrefix++;

      while ((escapeOffset) >= (4096 << numPrefix))
        numPrefix++;
    }

    escapeBase <<= numPrefix;

    se->setLength(28 + (numPrefix << 1));
    se->setInfo(escapeBase | ((escapeOffset << 1) - escapeBase) | sign);
  }

  symbol2vlc(se);

  return se->getLength();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcEnc::encodeTotalZeros(SyntaxElement* se)
{
  int totZeros = se->getValue1();
  int numCoeff = se->getValue2();

  se->setLength(TotalZerosTableL[numCoeff-1][totZeros]);
  se->setInfo  (TotalZerosTableC[numCoeff-1][totZeros]);

  symbol2vlc(se);

  return se->getLength();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcEnc::encodeRuns(SyntaxElement *se)
{
  int runLength = se->getValue1();
  int zerosLeft = se->getValue2();

  se->setLength(RunTableL[zerosLeft][runLength]);
  se->setInfo  (RunTableC[zerosLeft][runLength]);

  symbol2vlc(se);

  return se->getLength();
}

