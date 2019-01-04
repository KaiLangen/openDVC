#include <cstdio>

#ifndef DECODER_INC_SIDEINFORMATION_H
#define DECODER_INC_SIDEINFORMATION_H

#include "config.h"

class CorrModel;
class Codec;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class SideInformation
{
public:
  SideInformation(Codec* codec);
  virtual ~SideInformation() {};
  virtual void createSideInfo(imgpel*, imgpel*, imgpel*,
                              int=0, int=0, int=0) = 0;

# if RESIDUAL_CODING
  void getResidualFrame(imgpel* bRefFrame, imgpel* fRefFrame,
                        imgpel* currFrame, int* residue, int* dirList);
  void getRecFrame(imgpel *imgFReference, imgpel *imgBReference, int *iResidue,
                   imgpel *imgRec,int *iList);
# endif

protected:
  void lowpassFilter(imgpel* src, imgpel* dst, const int boxSize);
  void lowpassFilter(imgpel* src, imgpel* dst, int width, int height, const int boxSize);

  void spatialSmooth(imgpel* imgPrev, imgpel* imgNext, mvinfo* varCandidate,
                     const int iBlockSize, const int iPadSize);

  void pad(imgpel *src,imgpel *dst, const int iPadSize);


protected:
  Codec*      _codec;
};

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class SI_MC : public SideInformation
{
public:
  SI_MC(Codec* codec, FILE* file, int srcHeight, int srcWidth);
  ~SI_MC();

  void createSideInfo(imgpel* prevTrg, imgpel* nextTrg, imgpel* currTrg,
                      int prevFrameno, int nextFrameNo, int currFrameNo);

protected:
  void ME(imgpel* prevFrame, imgpel* currFrame, imgpel* nextFrame,
          int width, int height);

  void MC(imgpel* prevTrg, imgpel* currTrg, imgpel* nextTrg, int factor);

protected:
  int _srcHeight;
  int _srcWidth;
  int _srcSize;

  int _trgHeight;
  int _trgWidth;
  int _trgSize;

  FILE* _file;
  mvinfo* _mvs;
  int _param;
  int _iRange;
  int _nMV;

  imgpel* _prevBuffer;
  imgpel* _currBuffer;
  imgpel* _nextBuffer;
  imgpel* _prevKeyFrame;
  imgpel* _nextKeyFrame;
  imgpel* _currFrame;
};

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class SI_MCI : public SideInformation
{
public:
  SI_MCI(Codec* codec, CorrModel* model, FILE* mvs = 0);

  ~SI_MCI();

# if SI_REFINEMENT
  void getRefinedSideInfo(imgpel *imgPrevKey,imgpel *imgNextKey,
                          imgpel *imgCurrFrame, imgpel *imgTmpRec,
                          imgpel *imgRefined, int iMode);
# endif

  void createSideInfo(imgpel* imgPreKey,imgpel* imgNextKey,
                      imgpel* imgCurrFrame, int a = 0 ,int b = 0, int c = 0);

protected:
  void readMVFromFile(mvinfo* varCand);
  void forwardME(imgpel* prev, imgpel* curr, mvinfo* candidate,const int range);
  void bidirectME(imgpel* prev, imgpel* next, mvinfo* candidate,
                  const int iPadSize, const int range);
  void MC(imgpel* imgPrev, imgpel* imgNext, imgpel* imgDst ,
          imgpel* mc1,imgpel* mc2, mvinfo* candidate, mvinfo* candidate2,
          const int iPadSize, const int range, const int mode);

  void getSkippedRecFrame(imgpel* imgPrevKey,imgpel * imgWZFrame, int* skipMask);

# if SI_REFINEMENT
  void createSideInfoProcess(imgpel* imgPrevKey,imgpel* imgNextKey,
                             imgpel* imgMCForward, imgpel* imgMCBackward,
                             int iMode);
  void getRefinedSideInfoProcess(imgpel* imgPrevBuffer,imgpel* imgTmpRec,
                                 imgpel* imgSI, imgpel* imgRefined,
                                 mvinfo* varList,int iMode);
# else
  void createSideInfoProcess(imgpel* imgPrevKey,imgpel* imgNextKey,
                             imgpel* imgMCForward, imgpel* imgMCBackward);
# endif

protected:
  CorrModel*  _model;
  FILE*       _mvFile;
  int         _nMV;
  mvinfo*     _varList0;
  mvinfo*     _varList1;

# if SI_REFINEMENT
  int*        _refinedMask;
# endif
};

void bilinear(imgpel *source, imgpel *buffer, int buffer_w, int buffer_h,
              int picwidth, int picheight, int px, int py);

#if SI_REFINEMENT
float getDCValue(imgpel* img,int iStep,int iStripe,int iBlock);
int calcDist(imgpel* blk1, imgpel* blk2, int width1, int width2,
             int s1, int s2, int blocksize);
#endif

#endif // DECODER_INC_SIDEINFORMATION_H

