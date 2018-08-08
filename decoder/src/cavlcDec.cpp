
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "cavlcDec.h"
#include "fileManager.h"
#include "codec.h"
#include "bitstream.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
CavlcDec::CavlcDec(Codec* codec, int blockSize) : Cavlc(codec, blockSize)
{
  for (int i = 0; i < _codec->getBitPlaneLength(); i++)
    _mbs[i].nnz[0][0] = 0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcDec::decode(int* iDCT, int ix, int iy)
{
  int width = _codec->getFrameWidth();
  int qp = _codec->getQp();
  int iBitCount=0;
  int numCoeff,t1s,totalZeros,zerosLeft;
  int nc,vlc;

  int index = (ix/4)+(iy/4)*(width/4);

  numCoeff = t1s = totalZeros = zerosLeft = 0;

  if (ix == 0 && iy == 0)
    nc = 0;
  else if (ix == 0 && iy != 0)
    nc = getNumNonzero(ix, iy-4);
  else if (ix != 0 && iy == 0)
    nc = getNumNonzero(ix-4, iy);
  else
    nc = (getNumNonzero(ix, iy-4) + getNumNonzero(ix-4, iy) + 1) >> 1;

  if (nc < 2)
    vlc = 0;
  else if (nc < 4)
    vlc = 1;
  else if (nc < 8)
    vlc = 2;
  else
    vlc = 3;

  iBitCount += decodeNumCoeffTrailingOnes(numCoeff, t1s, vlc);

  _mbs[index].nnz[0][0] = numCoeff;

  if (numCoeff == 0)
    return iBitCount;

  iBitCount += decodeLevel(numCoeff, t1s, _mbs[index].level, _mbs[index].Ones);

  _mbs[index].nnz[0][0] = numCoeff;

  if (numCoeff == 16) {
    totalZeros = 0;
    zerosLeft  = 0;
  }
  else {
    iBitCount += decodeTotalZero(totalZeros, numCoeff);
    zerosLeft  = totalZeros;
  }

  for (int i = 0; i < (numCoeff-1); i++) {
    if (zerosLeft > 0)
      iBitCount += decodeRun(_mbs[index].iRun[i], zerosLeft);
    else
      _mbs[index].iRun[i] = 0;

    zerosLeft -= _mbs[index].iRun[i];
  }

  _mbs[index].iRun[numCoeff-1] = zerosLeft;

  int iSign;
  int iCoeffNum = -1;

#if MODE_DECISION
  iCoeffNum += _codec->getNumChnCodeBands();
#endif

  for (int i = numCoeff-1; i >= 0; i--) {
    iCoeffNum += _mbs[index].iRun[i] + 1;

    int x = ScanOrder[iCoeffNum][0];
    int y = ScanOrder[iCoeffNum][1];

  //  iDCT[ix+x + (iy+y)*width] = _mbs[index].level[i];
    if (x == 0 && y == 0) {
#if RESIDUAL_CODING
      iSign = (_mbs[index].level[i] >= 0) ? 0 : 1;
      iDCT[ix+x + (iy+y)*width] = abs(_mbs[index].level[i]);

      if (iSign == 1)
        iDCT[ix+x + (iy+y)*width] |= (0x1 << (_codec->getQuantMatrix(qp, x, y)-1));
#else
      iDCT[ix+x + (iy+y)*width] = _mbs[index].level[i];
#endif
    }
    else {
      iSign = (_mbs[index].level[i] >= 0) ? 0 : 1 ;
      iDCT[ix+x + (iy+y)*width] = abs(_mbs[index].level[i]);

      if (iSign == 1)
        iDCT[ix+x + (iy+y)*width] |= (0x1 << (_codec->getQuantMatrix(qp, x, y)-1));
    }
  }

  return iBitCount;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcDec::decodeNumCoeffTrailingOnes(int& numCoeff, int& t1s, int vlc)
{
  int value;

  int length = 0;

  if (vlc < 3) {
    value = 0;

    for (length = 1; length < 17; length++) {
      value <<= 1;
      value |= _codec->getBitstream()->read(1);

      for (int j = 0; j < 4; j++)
        for (int i = 0; i < 17; i++) {
          if (NumVlcTableC[vlc][j][i] == value &&
              NumVlcTableL[vlc][j][i] == length) {
            numCoeff = i;
            t1s = j;
            goto DecodeNumTrailDone;
          }
        }
    }

    DecodeNumTrailDone:
      ;
  }
  else {
    length = 6;
    value = _codec->getBitstream()->read(6);

    if (value == 3)
      numCoeff = t1s = 0;
    else {
      numCoeff = (value>>2) + 1;
      t1s = value & 0x3;
    }
  }

  return length;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcDec::decodeLevel(int iNumCoef, int iTrailingOnes, int* iLevel, int* iOnes)
{
  int iValue = 0;
  int iLevelPrefix,iSuffixLength,iLevelSuffixSize;
  unsigned int iLevelSuffix;
  int iLevelCode;
  int iIndex = 0;

  int length = 0;

  //decode trailingOnes
  for (int i = 0; i < iTrailingOnes; i++) {
    int b = _codec->getBitstream()->read(1);
    length++;

    if (b == 0)
      iLevel[iIndex] = 1;
    else
      iLevel[iIndex] = -1;

    iIndex++;
  }

  //initialize SuffixLength
  if (iNumCoef > 10 && iTrailingOnes < 3)
    iSuffixLength = 1;
  else
    iSuffixLength = 0;

  for(int i = 0; i < iNumCoef-iTrailingOnes; i++) {
    iLevel[iIndex] = 0;

    //get level prefix
    iLevelPrefix = -1;

    for (int b = 0; !b; iLevelPrefix++) {
      b = _codec->getBitstream()->read(1);
      length++;
    }

    if (iLevelPrefix == 14 && iSuffixLength == 0)
      iLevelSuffixSize = 4;
    else if (iLevelPrefix >= 15)
      iLevelSuffixSize = iLevelPrefix - 3;
    else
      iLevelSuffixSize = iSuffixLength;

    if (iLevelSuffixSize > 0) {
      iLevelSuffix = (unsigned int)_codec->getBitstream()->read(iLevelSuffixSize);
      length += iLevelSuffixSize;
    }
    else
      iLevelSuffix = 0;

    iLevelCode = (Min(15,iLevelPrefix)<<iSuffixLength) + iLevelSuffix;

    if (iLevelPrefix >= 15 && iSuffixLength == 0)
      iLevelCode += 15;
    if (iLevelPrefix >= 16)
      iLevelCode += ((1<<( iLevelPrefix-3 ))-4096);
    if (iIndex == iTrailingOnes && iTrailingOnes < 3)
      iLevelCode += 2;

    if (iLevelCode % 2 == 0) //even number
      iLevel[iIndex] = (iLevelCode + 2 ) >> 1;
    else
      iLevel[iIndex] = (-iLevelCode-1) >> 1;

    if (iSuffixLength == 0)
      iSuffixLength = 1;
    if((abs(iLevel[iIndex]) > (3<<(iSuffixLength-1))) && (iSuffixLength < 6))
      iSuffixLength++;

    iIndex++;
  }

  return length;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcDec::decodeTotalZero(int& iTotalZeros, int iNumCoef)
{
  int iValue = 0;

  int length = 0;

  for (length = 1; length < 10; length++) {
    int success = 0;
    iValue <<= 1;
    iValue |= _codec->getBitstream()->read(1);

    for (int i = 0; i < 16; i++) {
      if (TotalZerosTableC[iNumCoef-1][i] == iValue && TotalZerosTableL[iNumCoef-1][i] == length) {
        success = 1;
        iTotalZeros = i;
        break;
      }
    }

    if (success == 1)
      break;
  }

  return length;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int CavlcDec::decodeRun(int& iRun, int iZerosLeft)
{
  int iValue = 0;

  int length = 0;

  for (length = 1; length < 10; length++) {
    int success = 0;
    iValue <<= 1;
    iValue |= _codec->getBitstream()->read(1);

    for (int i = 0; i < 15; i++) {
      if (RunTableC[iZerosLeft-1][i] == iValue && RunTableL[iZerosLeft-1][i] == length) {
        success = 1;
        iRun = i;
        break;
      }
    }

    if (success == 1)
      break;
  }

  return length;
}

