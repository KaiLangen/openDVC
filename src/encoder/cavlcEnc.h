
#ifndef ENCODER_INC_CAVLCENC_H
#define ENCODER_INC_CAVLCENC_H

#include <vector>

#include <cstdio>

#include "defs.h"
#include "cavlc.h"

using std::vector;

class File;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class SyntaxElement
{
public:
  SyntaxElement(int value1 = 0, int value2 = 0, int length = 0, int info = 0) {
    _value1 = value1;
    _value2 = value2;
    _length = length;
    _info   = info;
  };

  void set(int value1, int value2, int length, int info) {
    _value1 = value1;
    _value2 = value2;
    _length = length;
    _info   = info;
  };

  int getValue1() { return _value1; };
  int getValue2() { return _value2; };
  int getLength() { return _length; };
  int getInfo()   { return _info; };

  void setValue1(int value1) { _value1 = value1; };
  void setValue2(int value2) { _value2 = value2; };
  void setLength(int length) { _length = length; };
  void setInfo  (int info)   { _info   = info; };

private:
  int _value1;  //!< numerical value of syntax element
  int _value2;  //!< for blocked symbols, e.g. run/level
  int _length;  //!< length of code
  int _info;    //!< info part of UVLC code
};

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class CavlcEnc : public Cavlc
{
public:
  CavlcEnc(Codec* codec, int blockSize);

  int encode(int* frame, int* skipMask);

private:
  void setupMacroBlock(int* frame, int mbX, int mbY, int c);
  int encodeMacroBlock(int mbX, int mbY, int c);

  int symbol2vlc(SyntaxElement* sym, int c);

  int encodeNumTrail(SyntaxElement* se, int c);
  int encodeSignTrail(vector<int>& sign, int c);
  int encodeLevelsVlc0(SyntaxElement* se, int c);
  int encodeLevelsVlcN(SyntaxElement* se, int vlc, int c);
  int encodeTotalZeros(SyntaxElement* se, int c);
  int encodeRuns(SyntaxElement* se, int c);

  File*     _patternFile;
  FILE*     _patternFh;
};

#endif // ENCODER_INC_CAVLCENC_H

