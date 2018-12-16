
#ifndef DECODER_INC_DECODER_H
#define DECODER_INC_DECODER_H

#include <iostream>

#include "config.h"
#include "codec.h"

using namespace std;

class FileManager;
class Transform;
class CorrModel;
class SideInformation;
class CavlcDec;
class FrameBuffer;
class LdpcaDec;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class Decoder : public Codec
{
public:
  Decoder(char** argv);
  ~Decoder() { /* TODO Remember to free memory space */ };

  void decodeWZframe();

  int*  getSpiralSearchX()     { return _spiralSearchX; };
  int*  getSpiralSearchY()     { return _spiralSearchY; };
  int*  getSpiralHpelSearchX() { return _spiralHpelSearchX; };
  int*  getSpiralHpelSearchY() { return _spiralHpelSearchY; };

private:
  void initialize();

  void decodeWzHeader();

  void parseKeyStat(const char* filename, double& rate, double& psnr, int& QP);

  int getSyndromeData();
  int decodeSkipMask();

  void getSourceBit(int* dct_q, double* source, int q_i, int q_j, int curr_pos);
  double decodeLDPC(int* iQuantDCT, int* iDCT, int* iDecoded, int x, int y, int iOffset);

  void motionSearchInit(int maxsearch_range);

private:
  FileManager*      _files;

  FrameBuffer*      _fb;

  char              _channel;

  Transform*        _trans;

  CorrModel*        _model;
  SideInformation*  _si;

  CavlcDec*         _cavlc;
  LdpcaDec*         _ldpca;

  int               _maxValue[4][4];
  int*              _skipMask;

# if RESIDUAL_CODING | MODE_DECISION
  int               _rcBitPlaneNum;
  int*              _rcList;
  int               _rcQuantMatrix[4][4];
# endif

  int*              _spiralSearchX;
  int*              _spiralSearchY;
  int*              _spiralHpelSearchX;
  int*              _spiralHpelSearchY;
};

void decodeBits(double *LLR_intrinsic, double *accumulatedSyndrome, double *source,
                double *decoded, double *rate, double *numErrors,unsigned char crccode,int numcode);
int  beliefPropagation(int *ir, int *jc, int m, int n, int nzmax,
                       double *LLR_intrinsic, double *syndrome,
                       double *decoded);

bool checkCRC(double * source,const int length,unsigned char crc);

double calcPSNR(unsigned char* img1,unsigned char* img2,int length);

int getSymbol(int len,int &curr_pos,char *buffer);

#endif // DECODER_INC_DECODER_H

