
#include <cstdio>

#include "ldpca.h"
#include "codec.h"

Ldpca::Ldpca(const string& fileName, Codec* codec)
{
  _codec = codec;

  FILE* fh;

  fh = fopen(fileName.c_str(), "r");

  fscanf(fh, "%d", &_numCodes);
  fscanf(fh, "%d", &_n);
  fscanf(fh, "%d", &_nzmax);
  fscanf(fh, "%d", &_totalNumInc);

  _ir     = new int*[_numCodes];
  _jc     = new int [_n+1];
  _numInc = new int [_numCodes];
  _txSeq  = new int*[_numCodes];

  for (int i = 0; i < _numCodes; i++) {
    _ir[i]    = new int[_nzmax];
    _txSeq[i] = new int[_totalNumInc];
  }

  for (int k = 0; k < (_n+1); k++)
    fscanf(fh, "%d", _jc+k);

  for (int k = 0; k < _numCodes; k++) {
    fscanf(fh, "%d", &_numInc[k]);

    for (int l = 0; l < _numInc[k]; l++)
      fscanf(fh, "%d", _txSeq[k]+l);

    for (int l = 0; l < _nzmax; l++)
      fscanf(fh, "%d", _ir[k]+l);
  }

  fclose(fh);
}

Ldpca::~Ldpca()
{
delete [] _jc;
delete [] _numInc;

for (int i = 0; i < _numCodes; i++) {
  delete [] _ir[i];
  delete [] _txSeq[i];
}
delete [] _ir;
delete [] _txSeq;
}
