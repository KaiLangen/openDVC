
#ifndef ENCODER_INC_LDPCAENC_H
#define ENCODER_INC_LDPCAENC_H

#include <string>

#include "ldpca.h"

using std::string;

class LdpcaEnc : public Ldpca
{
public:
  LdpcaEnc(const string& fileName, Codec* codec);

  void encode(int* source, bool* accumulatedSyndrome);

private:
};

#endif // ENCODER_INC_LDPCAENC_H

