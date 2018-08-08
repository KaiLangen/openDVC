
#ifndef COMMON_INC_LDPCA_H
#define COMMON_INC_LDPCA_H

#include <string>

#include "config.h"

using std::string;

class Codec;

class Ldpca
{
public:
  Ldpca(const string& fileName, Codec* codec);

protected:
  Codec*  _codec;

  int     _numCodes;
  int     _n;
  int     _nzmax;
  int     _totalNumInc;
  int**   _ir;
  int*    _jc;
  int**   _txSeq;
  int*    _numInc;
};

#endif // COMMON_INC_LDPCA_H

