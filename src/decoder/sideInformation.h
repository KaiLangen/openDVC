
#ifndef DECODER_INC_SIDEINFORMATION_H
#define DECODER_INC_SIDEINFORMATION_H

#include "defs.h"

class CorrModel;
class Codec;

class SideInformation
{
public:
  SideInformation(Codec* codec, CorrModel* model);

# if SI_REFINEMENT
  void getRefinedSideInfo(imgpel *imgPrevKey, imgpel *imgNextKey,
                          imgpel *imgCurrFrame, imgpel *imgTmpRec,
                          imgpel *imgRefined,int iMode, int c);
# endif

  void createSideInfo(imgpel* imgPreKey,imgpel* imgNextKey,
                      imgpel* imgCurrFrame, std::FILE* mvFilePtr, int c);

# if RESIDUAL_CODING
  void getResidualFrame(imgpel* bRefFrame, imgpel* fRefFrame,
                        imgpel* currFrame, int* residue, int* dirList, int c);
  void getRecFrame(imgpel *imgFReference, imgpel *imgBReference,
                   int *iResidue,imgpel *imgRec, int *iList, int c);
# endif

private:
  void lowpassFilter(imgpel* src, imgpel* dst, const int boxSize, int c);
  void forwardME(imgpel* prev, imgpel* curr, mvinfo* candidate,
                 const int range, int c);
  void bidirectME(imgpel* prev, imgpel* next, mvinfo* candidate,
                  const int iPadSize, const int range, int c);
  void MC(imgpel* imgPrev, imgpel* imgNext, imgpel* imgDst ,imgpel* mc1,
          imgpel* mc2, mvinfo* candidate, mvinfo* candidate2, const int iPadSize,
          const int range, const int mode, int c);

  void spatialSmooth(imgpel* imgPrev, imgpel* imgNext, mvinfo* varCandidate,
                     const int iBlockSize, const int iPadSize, int c);

  void pad(imgpel *src,imgpel *dst, const int iPadSize, int c);

  void getSkippedRecFrame(imgpel* imgPrevKey, imgpel* imgWZFrame,
                          int* skipMask, int c);

# if SI_REFINEMENT
  void createSideInfoProcess(imgpel* imgPrevKey, imgpel* imgNextKey,
                             imgpel* imgMCForward, imgpel* imgMCBackward,
                             int iMode, std::FILE* mvFilePtr, int c);
  void getRefinedSideInfoProcess(imgpel* imgPrevBuffer,imgpel* imgTmpRec,
                                 imgpel* imgSI,imgpel* imgRefined,
                                 mvinfo* varList,int iMode, int c);
# else
  void createSideInfoProcess(imgpel* imgPrevKey, imgpel* imgNextKey, 
                             imgpel* imgMCForward, imgpel* imgMCBackward,
                             std::FILE* mvFilePtr, int c);
# endif

  Codec*      _codec;

  CorrModel*  _model;

  mvinfo*     _varList0[NCHANS];
  mvinfo*     _varList1[NCHANS];

# if SI_REFINEMENT
  int*        _refinedMask[NCHANS];
# endif
};

void bilinear(imgpel *source, imgpel *buffer, int buffer_w, int buffer_h,
              int picwidth, int picheight, int px, int py, int c);

#if SI_REFINEMENT
float getDCValue(imgpel* img, int iStep, int iStripe, int iBlock);
int calcDist(imgpel* blk1, imgpel* blk2, int width1, int width2,
             int s1,int s2,int blocksize);
#endif

#endif // DECODER_INC_SIDEINFORMATION_H

