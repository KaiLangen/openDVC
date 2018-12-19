#include <vector>
#include <cstdlib>

#include "codec.h"
#include "ME.h"
#include "sideInformation.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
SideInformation::SideInformation(Codec* codec)
{
  _codec = codec;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::lowpassFilter(imgpel* src, imgpel* dst, const int boxSize)
{
  int width  = _codec->getFrameWidth();
  int height = _codec->getFrameHeight();

  const int left  = boxSize/2;
  const int right = boxSize - (left + 1);

  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++) {
      int sum   = 0;
      int total = 0;

      // Sum up pixels within the box
      for (int j = y-left; j <= y+right; j++)
        for (int i = x-left; i <= x+right; i++)
          if (i >= 0 && i < width &&
              j >= 0 && j < height) {
            sum += src[i+j*width];
            total++;
          }

      dst[x+y*width] = (imgpel)(sum/total);
    }
}

void SideInformation::lowpassFilter(imgpel* src, imgpel* dst,
                                    int width, int height, const int boxSize)
{
  const int left  = boxSize/2;
  const int right = boxSize - (left + 1);

  for (int y = 0; y < height; y++)
    for (int x = 0; x < width; x++) {
      int sum   = 0;
      int total = 0;

      // Sum up pixels within the box
      for (int j = y-left; j <= y+right; j++)
        for (int i = x-left; i <= x+right; i++)
          if (i >= 0 && i < width &&
              j >= 0 && j < height) {
            sum += src[i+j*width];
            total++;
          }

      dst[x+y*width] = (imgpel)(sum/total);
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*
* Spatial Smoothing
*/
void SideInformation::spatialSmooth(imgpel* imgPrev, imgpel* imgNext, mvinfo* varCandidate, const int iBlockSize, const int iPadSize)
{
  int iWidth,iHeight;
  iWidth  = _codec->getFrameWidth();
  iHeight = _codec->getFrameHeight();
  int iIndex[9];
  double dWeight[9];
  int iSAD[9];
  double dMinWeight;
  int iBestMVx=0,iBestMVy=0;
  int iPx[2],iPy[2];
  int iBestIdx=0;

  imgpel *imgPrevBuffer,*imgBufferNext;
  mvinfo *varRefine = new mvinfo[iWidth*iHeight/(iBlockSize*iBlockSize)];
  imgPrevBuffer     = new imgpel[(iWidth+2*iPadSize)*(iHeight+2*iPadSize)*4];
  imgBufferNext     = new imgpel[(iWidth+2*iPadSize)*(iHeight+2*iPadSize)*4];

  bilinear(imgPrev,imgPrevBuffer,iWidth+2*iPadSize,iHeight+2*iPadSize,iWidth+2*iPadSize,iHeight+2*iPadSize,0,0);
  bilinear(imgNext,imgBufferNext,iWidth+2*iPadSize,iHeight+2*iPadSize,iWidth+2*iPadSize,iHeight+2*iPadSize,0,0);

  for (int j = 0; j < iHeight/iBlockSize; j++)
    for (int i = 0; i < iWidth/iBlockSize; i++) {
      for (int k = 0; k < 9; k++)
        iIndex[k] = -1;

      if(j>0)                           iIndex[0]=i+(j-1)*(iWidth/iBlockSize);
      if(i>0)                           iIndex[1]=i-1+j*(iWidth/iBlockSize);
      if(i<iWidth/iBlockSize-1)         iIndex[2]=i+1+j*(iWidth/iBlockSize);
      if(j<iHeight/iBlockSize-1)        iIndex[3]=i+(j+1)*(iWidth/iBlockSize);
      if(i>0 && j>0)                    iIndex[5]=(i-1)+(j-1)*(iWidth/iBlockSize);
      if(i>0 && j<iHeight/iBlockSize-1) iIndex[6]=(i-1)+(j+1)*(iWidth/iBlockSize);
      if(i<iWidth/iBlockSize-1 && j>0)  iIndex[7]=(i+1)+(j-1)*(iWidth/iBlockSize);
      if(i<iWidth/iBlockSize-1 && j<iHeight/iBlockSize-1) iIndex[8]=(i+1)+(j+1)*(iWidth/iBlockSize);

      iIndex[4]=i+j*(iWidth/iBlockSize);

      varRefine[iIndex[4]].iMvx = varCandidate[iIndex[4]].iMvx;
      varRefine[iIndex[4]].iMvy = varCandidate[iIndex[4]].iMvy;
      varRefine[iIndex[4]].iCx  = varCandidate[iIndex[4]].iCx;
      varRefine[iIndex[4]].iCy  = varCandidate[iIndex[4]].iCy;

      for (int k = 0; k < 9; k++) {
        if (iIndex[k] != -1) {
          iPx[0]=2*(varCandidate[iIndex[4]].iCx+iPadSize)+varCandidate[iIndex[k]].iMvx;
          iPy[0]=2*(varCandidate[iIndex[4]].iCy+iPadSize)+varCandidate[iIndex[k]].iMvy;
          iPx[1]=2*(varCandidate[iIndex[4]].iCx+iPadSize)-varCandidate[iIndex[k]].iMvx;
          iPy[1]=2*(varCandidate[iIndex[4]].iCy+iPadSize)-varCandidate[iIndex[k]].iMvy;

          iSAD[k] = 0;
          iSAD[k] = calcSAD(imgPrevBuffer+iPx[0]+iPy[0]*2*(2*iPadSize+iWidth),
                            imgBufferNext+iPx[1]+iPy[1]*2*(2*iPadSize+iWidth),
                            2*(iWidth+2*iPadSize), 2*(iWidth+2*iPadSize), 2, 2, iBlockSize);
        }
      }

      dMinWeight = -1;

      for (int il = 0; il < 9; il++) {
        dWeight[il] = 0;

        if (iIndex[il] != -1) {
          for (int im = 0; im < 9; im++) {
            if (il != im && iIndex[im] != -1)
              dWeight[il] += (double)(iSAD[il]/std::max<double>(iSAD[im],0.0001))
                               *( abs(varCandidate[iIndex[il]].iMvx-varCandidate[iIndex[im]].iMvx)
                                 +abs(varCandidate[iIndex[il]].iMvy-varCandidate[iIndex[im]].iMvy));
          }
        }
        else
          dWeight[il]=-1;

        if ((dMinWeight<0 || dWeight[il]<=dMinWeight) && dWeight[il]>=0) {
          iBestMVx = varCandidate[iIndex[il]].iMvx;
          iBestMVy = varCandidate[iIndex[il]].iMvy;
          dMinWeight = dWeight[il];
          iBestIdx = il;
        }
      }

      varRefine[iIndex[4]].iMvx = iBestMVx;
      varRefine[iIndex[4]].iMvy = iBestMVy;
      varRefine[iIndex[4]].fDist = (float)iSAD[iBestIdx];
    }

  for (int iIndex = 0; iIndex < iWidth*iHeight/(iBlockSize*iBlockSize); iIndex++) {
    varCandidate[iIndex].iMvx = varRefine[iIndex].iMvx;
    varCandidate[iIndex].iMvy = varRefine[iIndex].iMvy;
    varCandidate[iIndex].iCx  = varRefine[iIndex].iCx;
    varCandidate[iIndex].iCy  = varRefine[iIndex].iCy;
    varCandidate[iIndex].fDist= varRefine[iIndex].fDist;
  }

  delete [] imgPrevBuffer;
  delete [] imgBufferNext;
  delete [] varRefine;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::pad(imgpel* src, imgpel* dst, const int iPadSize)
{
  int width  = _codec->getFrameWidth();
  int height = _codec->getFrameHeight();

  int paddedWidth  = width  + 2*iPadSize;
  int paddedHeight = height + 2*iPadSize;

  // Upper left
  for (int y = 0; y < iPadSize; y++)
    for (int x = 0; x < iPadSize; x++)
      dst[x+y*paddedWidth] = src[0];

  // Upper
  for (int y = 0; y < iPadSize; y++)
    for (int x = iPadSize; x < iPadSize+width; x++)
      dst[x+y*paddedWidth] = src[x-iPadSize];

  // Upper right
  for (int y = 0; y < iPadSize; y++)
    for (int x = iPadSize+width; x < paddedWidth; x++)
      dst[x+y*paddedWidth] = src[width-1];

  // Left
  for (int y = iPadSize; y < iPadSize+height; y++)
    for (int x = 0; x < iPadSize; x++)
      dst[x+y*paddedWidth] = src[(y-iPadSize)*width];

  // Middle
  for (int y = iPadSize; y < iPadSize+height; y++)
    for (int x = iPadSize; x < iPadSize+width; x++)
      dst[x+y*paddedWidth] = src[x-iPadSize+(y-iPadSize)*width];

  // Right
  for (int y = iPadSize; y < iPadSize+height; y++)
    for (int x = iPadSize+width; x < paddedWidth; x++)
      dst[x+y*paddedWidth] = src[(y-iPadSize+1)*width-1];

  // Bottom left
  for (int y = iPadSize+height; y < paddedHeight; y++)
    for (int x = 0; x < iPadSize; x++)
      dst[x+y*paddedWidth] = src[(height-1)*width];

  // Bottom
  for (int y = iPadSize+height; y < paddedHeight; y++)
    for (int x = iPadSize; x < iPadSize+width; x++)
      dst[x+y*paddedWidth] = src[(x-iPadSize)+(height-1)*width];

  // Bottom right
  for (int y = iPadSize+height; y < paddedHeight; y++)
    for (int x = iPadSize+width; x < paddedWidth; x++)
      dst[x+y*paddedWidth] = src[height*width-1];

}

#if RESIDUAL_CODING
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::getResidualFrame(imgpel* bRefFrame, imgpel* fRefFrame, imgpel* currFrame, int* residue, int* dirList)
{
  int width  = _codec->getFrameWidth();
  int height = _codec->getFrameHeight();
  int iBlock = 8;
  int blockCount = 0;

  for (int j = 0; j < height; j += iBlock)
    for (int i = 0; i < width; i += iBlock) {
      imgpel* refFrame = (dirList[blockCount] == 0) ? bRefFrame : fRefFrame;

      for (int y = 0; y < iBlock; y++)
        for (int x = 0; x < iBlock; x++) {
          int idx = (i+x) + (j+y)*width;

          residue[idx] = currFrame[idx] - refFrame[idx];
        }

      blockCount++;
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SideInformation::getRecFrame(imgpel *imgFReference,imgpel *imgBReference,int *iResidue,imgpel *imgRec,int *iList)
{
  int iWidth,iHeight;
  iWidth  = _codec->getFrameWidth();
  iHeight = _codec->getFrameHeight();
  int iBlock = 8;
  int iIndex=0;
  int iPos;
  imgpel* imgRef;
  for(int j=0;j<iHeight;j+=iBlock)
    for(int i=0;i<iWidth;i+=iBlock)
    {
      imgRef = (iList[iIndex]==0)? imgFReference:imgBReference;

      for(int y=0;y<iBlock;y++)
        for(int x=0;x<iBlock;x++)
        {
          iPos = (i+x) + (j+y)*iWidth;
          imgRec[iPos]=Clip(0,255,iResidue[iPos]+imgRef[iPos]);
        }
      iIndex++;
    }

}

#endif



// FRIEND FUNCTIONS
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void bilinear(imgpel *source,imgpel *buffer,int buffer_w,int buffer_h,int picwidth,int picheight,int px,int py){
  for(int j=0;j<buffer_h;j++)
    for(int i=0;i<buffer_w;i++)
    {
      int buffer_r=2*buffer_w;
      int a,b,c,d;

      int x=px+i;
        int y=py+j;
      if(x>picwidth-1)x=picwidth-1;
      if(x<0)x=0;
      if(y>picheight-1)y=picheight-1;
      if(y<0)y=0;
      a=source[(x)+(y)*picwidth];

      if((x+1)<picwidth) b=source[(x+1)+(y)*picwidth];
      else b=a;

      if((y+1)<picheight) c=source[(x)+(y+1)*picwidth];
      else c=a;

      if((x+1)<picwidth && (y+1)<picheight) d=source[(x+1)+(y+1)*picwidth];
      else d=a;

      buffer[2*i+(2*j)*buffer_r]=a;
      buffer[(2*i+1)+(2*j)*buffer_r]=(a+b)/2;
      buffer[2*i+(2*j+1)*buffer_r]=(a+c)/2;
      buffer[(2*i+1)+(2*j+1)*buffer_r]=(a+b+c+d)/4;
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
float getDCValue(imgpel* img,int iStep,int iStripe,int iBlock)
{
  float fDCValue;
  int sum = 0;
  for(int j=0;j<iBlock;j++)
    for(int i=0;i<iBlock;i++)
    {
      sum += *(img + i*iStep + j*iStep*iStripe);
    }
  fDCValue = (float)sum/float(iBlock*iBlock);
  return fDCValue;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int calcDist(imgpel* blk1, imgpel* blk2, int width1, int width2, int s1,int s2,int blocksize){
  int iDist=0;
  for(int y=0;y<blocksize;y++)
    for(int x=0;x<blocksize;x++)
    {
      imgpel pel1=*(blk1+s1*x+s1*y*width1);
      imgpel pel2=*(blk2+s2*x+s2*y*width2);
      iDist+=(pel1-pel2)*(pel1-pel2);
    }
  return iDist;
}
