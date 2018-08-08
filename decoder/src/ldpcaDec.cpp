
#include <cstdio>
#include <cmath>
#include <cstring>

#include "ldpcaDec.h"
#include "codec.h"

LdpcaDec::LdpcaDec(const string& fileName, Codec* codec) : Ldpca(fileName, codec)
{
# if INVERSE_MATRIX
  FILE* fh;

  int width          = _codec->getFrameWidth();
  int height         = _codec->getFrameHeight();
  int bitPlaneLength = _codec->getBitPlaneLength();

  if (width == 352 && height == 288)
#   if HARDWARE_LDPC
    fh = fopen("ldpca/Inverse_Matrix_H_Reg1584.dat", "r");
#   else // if !HARDWARE_LDPC
    fh = fopen("ldpca/Inverse_Matrix_H_Reg6336.dat", "r");
#   endif // HARDWARE_LDPC
  else
    fh = fopen("ldpca/Inverse_Matrix_H_Reg1584.dat", "r");

  _invMatrix = new int[1584*1584];

  for (int i = 0; i < 1584*1584; i++)
    fscanf(fh, "%d", _invMatrix+i);

  fclose(fh);
# endif // INVERSE_MATRIX
}

//decode() finds the minimum rate for which the decoded varBitstream matches
//the transmitted portion of the accumulated syndrome.
//The number of residual bit errors is also calculated.
void LdpcaDec::decode(double *LLR_intrinsic, double *accumulatedSyndrome, double *source,
                      double *decoded, double *rate, double *numErrors,unsigned char crccode,int numcode)
{
   // FILE *fp;
    int n, m, nzmax, *ir, *jc;
    int numCodes, totalNumInc, numInc, *txSeq;
    int code, k, currIndex, prevIndex;
    double *syndrome;

    syndrome = new double[_n]; //actual length: m

  numCodes    = _numCodes;
  n           = _n;
  nzmax       = _nzmax;
  totalNumInc = _totalNumInc;
  jc          = _jc;
    //iterate through codes of increasing rate
    for(code=numcode-2; code<numCodes; code++)
    {
    numInc= _numInc[code];
    txSeq = _txSeq[code];
    ir    = _ir[code];
        m = (n/totalNumInc)*numInc;

        rate[0] = ((double) m)/((double) n);

        currIndex = txSeq[0];
        syndrome[0] = accumulatedSyndrome[currIndex];
        for(k=1; k<m; k++)
        {
            prevIndex = currIndex;
            currIndex = txSeq[k%numInc] + (k/numInc)*totalNumInc;
            syndrome[k] = (double) (((int) (accumulatedSyndrome[currIndex] + accumulatedSyndrome[prevIndex])) % 2);
        }
#if INVERSE_MATRIX
    //INVERSE_MATRIX
    if(rate[0]==1)
        {
      //translate accumulatedSyndrome into ori syndrome bits
      for(int j=n-1;j>0;j--)
      {
        accumulatedSyndrome[j]=accumulatedSyndrome[j]-accumulatedSyndrome[j-1];
        if(accumulatedSyndrome[j]!=0)accumulatedSyndrome[j]=1;
      }
      //use inverse matrix to reconstruct current g_iBitplaneNum
            for(int j=0;j<n;j++)
      {
        int tmp=0;
        for(int i=0;i<n;i++)
          tmp+=(_invMatrix[i*n+j]*(int)accumulatedSyndrome[i]);
        decoded[j]=tmp%2;
      }
      numErrors[0] = 0;
            return;
        }
#endif

        if((beliefPropagation(ir, jc, m, n, nzmax, LLR_intrinsic, syndrome, decoded) ))
    {
#if !HARDWARE_FLOW
      if(checkCRC(decoded,n,crccode))
#endif
      {
        numErrors[0] = 0;
        for(k=0; k<n; k++)
          numErrors[0] += (double) (decoded[k]!=source[k]);
        delete [] syndrome;
        return ;
      }
        }
    }
    delete [] syndrome;
}

//For implementation outline of beliefPropagation(), refer to
//W. E. Ryan, "An Introduction to LDPC Codes," in CRC Handbook for Coding
//and Signal Processing for Recording Systems (B. Vasic, ed.) CRC Press, 2004.
//available online (as of May 8, 2006) at
//http://www.ece.arizona.edu/~ryan/New%20Folder/ryan-crc-ldpc-chap.pdf

//beliefPropagation() runs several iterations belief propagation until
//either the decoded varBitstream agrees with the transmitted portion of
//accumulated syndrome or convergence or the max number of iterations.
//Returns 1 if decoded varBitstream agrees with
//transmitted portion of accumulated syndrome.
int LdpcaDec::beliefPropagation(int *ir, int *jc, int m, int n, int nzmax,
                       double *LLR_intrinsic, double *syndrome,
                       double *decoded)
{
    int iteration, k, l, sameCount;
    double *LLR_extrinsic, *check_LLR, *check_LLR_mag, *rowTotal, *LLR_overall;

    LLR_extrinsic = new double[nzmax];
    check_LLR = new double[nzmax];
    check_LLR_mag = new double[nzmax];
    rowTotal = new double[m];
    LLR_overall = new double[n];

    sameCount = 0;
    for(k=0; k<n; k++)
        decoded[k] = 0;

    //initialize variable-to-check messages
    for(k=0; k<n; k++)
        for(l=jc[k]; l<jc[k+1]; l++)
            LLR_extrinsic[l] = LLR_intrinsic[k];

    for(iteration=0; iteration<100; iteration++)
    {
        //Step 1: compute check-to-variable messages

        for(k=0; k<nzmax; k++)
        {
            check_LLR[k] = (double) ((LLR_extrinsic[k]<0) ? -1 : 1);
            check_LLR_mag[k] = ((LLR_extrinsic[k]<0) ? -LLR_extrinsic[k] : LLR_extrinsic[k]);
        }

        for(k=0; k<m; k++)
            rowTotal[k] = (double) ((syndrome[k]==1) ? -1 : 1);
        for(k=0; k<nzmax; k++)
            rowTotal[ir[k]] *= check_LLR[k];
        for(k=0; k<nzmax; k++)
            check_LLR[k] = check_LLR[k] * rowTotal[ir[k]];
            //sign of check-to-variable messages

        for(k=0; k<nzmax; k++)
            check_LLR_mag[k] = -log( tanh( Max(check_LLR_mag[k], 0.000000001)/2 ) );
        for(k=0; k<m; k++)
            rowTotal[k] = (double) 0;
        for(k=0; k<nzmax; k++)
            rowTotal[ir[k]] += check_LLR_mag[k];
        for(k=0; k<nzmax; k++)
            check_LLR_mag[k] = -log( tanh( Max(rowTotal[ir[k]] - check_LLR_mag[k], 0.000000001)/2 ) );
            //magnitude of check-to-variable messages

        for(k=0; k<nzmax; k++)
            check_LLR[k] = check_LLR[k] * check_LLR_mag[k];
            //check-to-variable messages

        //Step 2: compute variable-to-check messages

        for(k=0; k<n; k++)
        {
            LLR_overall[k] = LLR_intrinsic[k];
            for(l=jc[k]; l<jc[k+1]; l++)
                LLR_overall[k] += check_LLR[l];
        }

        for(k=0; k<n; k++)
            for(l=jc[k]; l<jc[k+1]; l++)
                LLR_extrinsic[l] = LLR_overall[k] - check_LLR[l];
                //variable-to-check messages

        //Step 3: test convergence and syndrome condition

        l = 0;
        for(k=0; k<n; k++)
            if(decoded[k] == ((LLR_overall[k]<0) ? 1 : 0))
                l++;
            else
                decoded[k] = ((LLR_overall[k]<0) ? 1 : 0);

        sameCount = ((l==n) ? sameCount+1 : 0);

        if(sameCount==5)
    {
      delete [] LLR_extrinsic;
      delete [] check_LLR;
      delete [] check_LLR_mag;
      delete [] rowTotal;
      delete [] LLR_overall;
            return 0; //convergence (to wrong answer)
    }
        for(k=0; k<m; k++)
            rowTotal[k] = syndrome[k];
        for(k=0; k<n; k++)
            for(l=jc[k]; l<jc[k+1]; l++)
                rowTotal[ir[l]] += decoded[k];

        for(k=0; k<m; k++)
            if(((int) rowTotal[k] % 2) != 0)
                break;
            else if(k==m-1)
      {
        delete [] LLR_extrinsic;
        delete [] check_LLR;
        delete [] check_LLR_mag;
        delete [] rowTotal;
        delete [] LLR_overall;
                return 1; //all syndrome checks satisfied
      }
    }

  delete [] LLR_extrinsic;
    delete [] check_LLR;
    delete [] check_LLR_mag;
    delete [] rowTotal;
    delete [] LLR_overall;

    return 0;
}
/*
*Do CRC checking
*Param
*/
bool LdpcaDec::checkCRC(double * source,const int length,unsigned char crc){


  //crc8 110011011
  const int code[9]={1,1,0,0,1,1,0,1,1};
  int *remainder = new int[length+8];
  unsigned char tmp=0;

  memset(remainder,0x0,4*(length+8));
  for(int i=0;i<length;i++) remainder[i]=(int)source[i];
  for(int i=0;i<8;i++)
  {
    remainder[i+length]=int(((crc)>>(7-i))&0x1);
  }
  for(int i=0;i<length;i++)
  {
    if(remainder[i]==1)
    {
      for(int j=0;j<9;j++)
      {
        remainder[i+j]=code[j]^remainder[i+j];
      }
    }
  }

  for(int i=0;i<length+8;i++)
  {
    if(remainder[i]!=0)
    {
      delete [] remainder;
      return false;
    }
  }

  delete [] remainder;
  return true;
}

