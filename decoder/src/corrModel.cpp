
#include <cmath>

#include "corrModel.h"
#include "transform.h"
#include "codec.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
double getLaplacian(int iLowerBound,int iUpperBound,double dAlpha,int iSideInfo,int iQuantSize,int iMode){
  double dProb = 0.0;
  double dx = iLowerBound;
  if(iMode==1)
  {
    do
    {
      dProb+=exp(-1*dAlpha*iQuantSize*fabs(dx-iSideInfo));
      dx+=1.0;
    }while(dx<iUpperBound);
  }
  else if(iMode==2)
  {
    if(iUpperBound!=iLowerBound)
    {
      if(iUpperBound>= iSideInfo && iSideInfo>= iLowerBound)
      {
        dProb=1.0-(exp(-1.0*dAlpha*(double(iSideInfo-iLowerBound)))+exp(-1.0*dAlpha*(double(iUpperBound-iSideInfo))))/2.0;
      }
      else if(iLowerBound>=iSideInfo)
      {
        dProb=(exp(-1.0*dAlpha*(double(iLowerBound-iSideInfo)))-exp(-1.0*dAlpha*(double(iUpperBound-iSideInfo))))/2.0;
      }
      else if(iSideInfo>=iUpperBound)
      {
        dProb=(exp(dAlpha*(double(iUpperBound-iSideInfo)))-exp(dAlpha*(double(iLowerBound-iSideInfo))))/2.0;
      }
    }
    else
    {
      dProb=0.5*dAlpha*exp(-1*dAlpha*fabs(double(iUpperBound-iSideInfo)));
    }
  }
  return dProb;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
double getCondProb(int iBit, int iBand, int iBitLength, int iCurrPos,
                   int iDecodedBit, double dAlpha, int iQuantSize,
                   int iMode){

  int iSign;
  int c;
  int lower;
  int upper;
  unsigned int value=0;

  if(iCurrPos==iBitLength-1)//msb
  {
    iSign=(iBit==1)?-1:1;
    c=(iSign==1)?0:-1*int(pow(2.0,double(iBitLength-1)));
    if(iSign==1)
    {
      upper=int(pow(2.0,double(iBitLength-1)))-1;
      lower=0;
    }
    else
    {
      upper=0;
      lower=-1*int(pow(2.0,double(iBitLength-1)))+1;
    }
    return getLaplacian(lower,upper,dAlpha,iBand,iQuantSize,iMode);
  }
  else
  {
    for(int i=iBitLength-2;i>iCurrPos;i--)
    {
      int b=(iDecodedBit>>i)&0x01;
      value|=(b*(0x01<<i));
    }
    value|=(iBit<<iCurrPos);
    iSign=((iDecodedBit>>(iBitLength-1))&0x01)?-1:1;
    c=(iSign==1)?0:-1*int(pow(2.0,double(iBitLength-1)));

    if(iSign==1)//positive
    {
      upper=(value+int(pow(2.0,double(iCurrPos)))-1);
      lower=(value);
    }
    else
    {
      lower=-1*(value+int(pow(2.0,double(iCurrPos)))-1);
      upper=-1*(value);
    }

    if(iMode==2)
    {
      if(iSign == 1)
      {
        upper=(value+int(pow(2.0,double(iCurrPos))))*iQuantSize-1;
        lower=(value)*iQuantSize;
      }
      else
      {
        lower=-1*((value+int(pow(2.0,double(iCurrPos))))*iQuantSize-1);
        upper=-1*(value)*iQuantSize;
      }
    }
    return getLaplacian(lower,upper,dAlpha,iBand,iQuantSize,iMode);
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
double CorrModel::getSoftInput(int* si, int* skipMask, int iCurrPos,
                               int *iDecoded, double* dLLR, int x, int y,
                               int iMode, int c){
  double p0,p1;
  double p=0.0;
  double tp=0.0;
  double entropy=0.0;
  int iIndex,iQuantSize,iBitLength;
  int iSideInfo,iQP;
  double* dAlpha;

  dAlpha = _codec->getAlpha(c);
  iQP = _codec->getQp();
  int mask;
  int   frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
  int  frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;
  int        bplen = (c == 0) ? Y_BPLEN : UV_BPLEN;

  for(int j=0;j<frameHeight/4;j++)
    for(int i=0;i<frameWidth/4;i++)
    {
      iIndex=i+j*(frameWidth/4);
      iQuantSize=_codec->getQuantStep(x, y, c);
#if RESIDUAL_CODING
      iBitLength= _codec->getQuantMatrix(iQP, x, y);

      if(iMode==2)
        iSideInfo = si[(i*4+x)+(j*4+y)*frameWidth];
      else
      {
        mask = (0x1<<(_codec->getQuantMatrix(iQP, x, y)-1))-1;
        int sign = (si[(i*4+x)+(j*4+y)*frameWidth]>>(_codec->getQuantMatrix(iQP, x, y)-1))&0x1;
        int value = si[(i*4+x)+(j*4+y)*frameWidth] & mask;

        iSideInfo = (sign==0)?value:-1*value;
      }
#else
      iBitLength=(x==0 && y==0) ? _codec->getQuantMatrix(iQP, x, y)+1 : _codec->getQuantMatrix(iQP, x, y);

      if(iMode==2 || (x==0 && y==0))
        iSideInfo = si[(i*4+x)+(j*4+y)*frameWidth];
      else
      {
        mask = (0x1<<(_codec->getQuantMatrix(iQP, x, y)-1))-1;
        int sign = (si[(i*4+x)+(j*4+y)*frameWidth]>>(_codec->getQuantMatrix(iQP, x, y)-1))&0x1;
        int value = si[(i*4+x)+(j*4+y)*frameWidth] & mask;

        iSideInfo = (sign==0)?value:-1*value;
      }
#endif

      p0=getCondProb(0,iSideInfo,iBitLength,iCurrPos,iDecoded[(i)+(j)*frameWidth/4],dAlpha[(i*4+x)+(j*4+y)*frameWidth],iQuantSize,iMode);
      p1=getCondProb(1,iSideInfo,iBitLength,iCurrPos,iDecoded[(i)+(j)*frameWidth/4],dAlpha[(i*4+x)+(j*4+y)*frameWidth],iQuantSize,iMode);

#if SKIP_MODE
      if(skipMask[iIndex]==1)//skip
        dLLR[iIndex]=500;
      else
#endif
      dLLR[iIndex]=log(p0/p1);

      p=p1/(p1+p0);
      entropy+= p*(log10(1.0/p)/log10(2.0))+(1.0-p)*(log10(1.0/(1-p))/log10(2.0));
      tp+=(p>0.5)?(1.0-p):(p);
    }

  entropy/=(bplen);
  tp/=(bplen);
  return 0.5*(entropy)*exp(entropy)+sqrt(tp*(1-tp));
}

/*
*update the sigma values for all dct bands
*Param
*/
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void CorrModel::correlationNoiseModeling(imgpel *imgMCFoward,
                                         imgpel* imgMCBackward, int c){
  double e,e2,r;
  int iIndex;
  double* dAlpha = _codec->getAlpha(c);
  double* dSigma = _codec->getSigma();
  double* dAverage = _codec->getAverage(c);
  int   frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
  int  frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;
  int        bplen = (c == 0) ? Y_BPLEN : UV_BPLEN;
  float *residue  = new float[frameWidth*frameHeight];
  float *residue_b  = new float[frameWidth*frameHeight];

  for(int j=0;j<frameHeight;j++)
    for(int i=0;i<frameWidth;i++)
    {
      residue[i+j*frameWidth] =
        float(imgMCFoward[i+j*frameWidth]-imgMCBackward[i+j*frameWidth])/2;
    }
  _trans->dctTransform(residue,residue_b, c);

  for(int j=0;j<4;j++)
    for(int i=0;i<4;i++)
    {
      e=e2=0.0;
      for(int y=0;y<frameHeight;y=y+4)
        for(int x=0;x<frameWidth;x=x+4)
        {
          iIndex=(x+i)+(y+j)*frameWidth;
          r=fabs((double)residue_b[iIndex]);
          e+=r;
          e2+=r*r;
        }

      e  /=bplen;
      e2 /=bplen;
      dSigma[i+j*4]=(e2-e*e);
      dAverage[i+j*4]=e;
    }
  for(int j=0;j<frameHeight;j=j+4)
    for(int i=0;i<frameWidth;i=i+4)
    {
      for(int y=0;y<4;y++)
        for(int x=0;x<4;x++)
        {
          iIndex=(x+i)+(y+j)*frameWidth;
          double d=fabs((double)residue_b[iIndex])-dAverage[x+y*4];
          int iQuantStep=_codec->getQuantStep(x, y, c);
          if(d*d<=dSigma[x+y*4])
          {
            dAlpha[iIndex]=sqrt(2/(dSigma[x+y*4]+0.1*iQuantStep*iQuantStep));
            //dAlpha[iIndex]=sqrt(2/(dSigma[x+y*4]));
          }
          else
          {
            dAlpha[iIndex]=sqrt(2/(d*d+0.1*iQuantStep*iQuantStep));
            //dAlpha[iIndex]=sqrt(2/(d*d));
          }
        }
    }
  delete [] residue;
  delete [] residue_b;

}


#if SI_REFINEMENT
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void CorrModel::updateCNM(imgpel* imgForward, imgpel* imgBackward,
                          int *refinedMask, int c)
{
  int iIndex;
  double* dAlpha = _codec->getAlpha(c);
  double* dSigma = _codec->getSigma();
  double* dAverage = _codec->getAverage(c);
  int   frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
  int  frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;
  int        bplen = (c == 0) ? Y_BPLEN : UV_BPLEN;

  float *residue  = new float[frameWidth*frameHeight];
  float *residue_b  = new float[frameWidth*frameHeight];

  for(int j=0;j<frameHeight;j=j+4)
    for(int i=0;i<frameWidth;i=i+4)
    {
      if(refinedMask[i/4+(j/4)*frameWidth/4])
      {
        for(int y=0;y<4;y++)
          for(int x=0;x<4;x++)
          {
            if(refinedMask[i/4+(j/4)*frameWidth/4]==2)
              residue[i+x+(j+y)*frameWidth] =
                float(imgForward[i+x+(j+y)*frameWidth]
                     -imgBackward[i+x+(j+y)*frameWidth])/2;
            else
              residue[i+x+(j+y)*frameWidth] =
                float(imgForward[i+x+(j+y)*frameWidth]
                    -imgBackward[i+x+(j+y)*frameWidth]);
          }

        _trans->dct4x4(residue,residue_b,i,j,c);

        for(int y=0;y<4;y++)
          for(int x=0;x<4;x++)
          {
            iIndex=(x+i)+(y+j)*frameWidth;
            double d=fabs((double)residue_b[iIndex])-dAverage[x+y*4];
            int iQuantStep=_codec->getQuantStep(x, y,c);
            if(d*d<=dSigma[x+y*4])
            {
              dAlpha[iIndex]=sqrt(2/(dSigma[x+y*4]+0.1*iQuantStep*iQuantStep));
            }
            else
            {
              dAlpha[iIndex]=sqrt(2/(d*d+0.1*iQuantStep*iQuantStep));
            }
          }
      }
    }
  delete [] residue;
  delete [] residue_b;
}
#endif

