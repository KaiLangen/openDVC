
#ifndef DECODER_INC_LDPCADEC_H
#define DECODER_INC_LDPCADEC_H

#include <string>

#include "ldpca.h"

using std::string;

class LdpcaDec : public Ldpca
{
public:
  LdpcaDec(const string& fileName, Codec* codec);

  void decode(double *LLR_intrinsic, double *accumulatedSyndrome, double *source,
                double *decoded, double *rate, double *numErrors,unsigned char crccode,int numcode);

private:
  int beliefPropagation(int *ir, int *jc, int m, int n, int nzmax,
                         double *LLR_intrinsic, double *syndrome,
                         double *decoded);
  bool checkCRC(double * source,const int length,unsigned char crc);

  int* _invMatrix;
};

#endif // DECODER_INC_LDPCADEC_H

