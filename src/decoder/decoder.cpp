
#include <sstream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// common headers
#include "bitstream.h"
#include "defs.h"
#include "fileManager.h"
#include "frameBuffer.h"
#include "regExp.h"
#include "time.h"
#include "transform.h"

// decoder headers
#include "cavlcDec.h"
#include "corrModel.h"
#include "decoder.h"
#include "ldpcaDec.h"
#include "sideInformation.h"

// cmake config header
#include "config.h"

using namespace std;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
Decoder::Decoder(char **argv)
{
  _files = FileManager::getManager();

  string wzFileName = argv[1];
  string recFileName = wzFileName.substr(0, wzFileName.find(".bin")) + ".yuv";

  _files->addFile("key",    argv[2])->openFile("rb");
  _files->addFile("origin", argv[3])->openFile("rb");
  _files->addFile("rec",    recFileName.c_str())->openFile("wb");
  _files->addFile("mv",     "mv.csv")->openFile("w");

  char chan[] = "yuv";
  for(int c = 0; c < NCHANS; c++) {
    string wzChanFile = wzFileName.substr(0, wzFileName.find(".bin"))
                         + "_" + chan[c] + ".bin";
    stringstream wzKey; 
    wzKey << "wz_" << chan[c];
    _files->addFile(wzKey.str(), wzChanFile)->openFile("rb");
    _bs[c] = new Bitstream(1024, _files->getFile(wzKey.str())->getFileHandle());
  }

  decodeWzHeader();

  initialize();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::initialize()
{
  _gop                   = 1 << _gopLevel;
  _fb                    = new FrameBuffer(_gop);
  _trans                 = new Transform(this);
  _model                 = new CorrModel(this, _trans);
  _si                    = new SideInformation(this, _model);
  _cavlc                 = new CavlcDec(this, 4);
  _sigma                 = new double[NBANDS];


  _dParity[0]            = new double[Y_BPLEN * BitPlaneNum[_qp]];
  _dParity[1]            = new double[UV_BPLEN * BitPlaneNum[_qp]];
  _dParity[2]            = new double[UV_BPLEN * BitPlaneNum[_qp]];
  _alpha[0]              = new double[Y_FSIZE];
  _alpha[1]              = new double[UV_FSIZE];
  _alpha[2]              = new double[UV_FSIZE];
  _skipMask[0]           = new int[Y_BPLEN];
  _skipMask[1]           = new int[UV_BPLEN];
  _skipMask[2]           = new int[UV_BPLEN];

  stringstream ladstream;
  stringstream datstream;

# if HARDWARE_LDPC
  _crc[0]                = new unsigned char[BitPlaneNum[_qp] * 4];
# else
  _crc[0]                = new unsigned char[BitPlaneNum[_qp]];
  ladstream << TOP_LVL_DIR << "/data/ldpca/6336_regDeg3.lad";
  datstream << TOP_LVL_DIR << "/data/ldpca/Inverse_Matrix_H_Reg6336.dat";
  _ldpca_cif             = new LdpcaDec(ladstream.str(), datstream.str(), this);
  ladstream.str("");
  datstream.str("");
# endif
  _crc[1]                = new unsigned char[BitPlaneNum[_qp]];
  _crc[2]                = new unsigned char[BitPlaneNum[_qp]];
  ladstream << TOP_LVL_DIR << "/data/ldpca/1584_regDeg3.lad";
  datstream << TOP_LVL_DIR << "/data/ldpca/Inverse_Matrix_H_Reg1584.dat";
  _ldpca                 = new LdpcaDec(ladstream.str(), datstream.str(), this);

  for (int c = 0; c < NCHANS; c++) {
    _numChnCodeBands[c]  = NBANDS;
    _average[c]          = new double[NBANDS];

# if RESIDUAL_CODING
  _rcList[c]             = new int[FSIZE/64];
  memset(_rcList[c], 0, FSIZE/64);
# endif
  }

  motionSearchInit(64);

}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
Decoder::~Decoder()
{
  delete [] _sigma;
  delete _fb;
  delete _trans;
  delete _cavlc;
  delete _ldpca;
# if !HARDWARE_LDPC
  delete _ldpca_cif;
#endif

  for (int c = 0; c < NCHANS; c++) {
    delete [] _parity[c];
    delete [] _alpha[c];
// TODO: fix this! _crc pointers are incremented, so unable to free. 
//                 Make an extra set of pointers, or use indices instead.
//                 delete _skipMask gets "double free or corruption" error.
//                 Investigate.
//    delete [] _crc[c];
//    delete [] _skipMask[c];
    delete [] _average[c];
    delete _bs[c];
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::decodeWzHeader()
{
  _qp           = _bs[0]->read(8);
  _numFrames    = _bs[0]->read(16);
  _gopLevel     = _bs[0]->read(2);

  cout << "--------------------------------------------------" << endl;
  cout << "WZ frame parameters" << endl;
  cout << "--------------------------------------------------" << endl;
  cout << "QP:     " << _qp << endl;
  cout << "Frames: " << _numFrames << endl;
  cout << "GOP:    " << (1<<_gopLevel) << endl;
  cout << "--------------------------------------------------" << endl << endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::decodeWZframe()
{
  double dPSNRAvg=0;
  double dPSNRSIAvg=0;

  clock_t timeStart, timeEnd;
  double cpuTime;
  int frameSize, bplen;

#if RESIDUAL_CODING
  int* iDCTBuffer      = new int [FSIZE];
  int* iDCTResidual    = new int [FSIZE];
#endif
  imgpel* oriCurrFrame = _fb->getCurrFrame();
  imgpel* imgSI        = _fb->getSideInfoFrame();

  int* iDCT            = _fb->getDctFrame();
  int* iDCTQ           = _fb->getQuantDctFrame();
  int* iDecoded        = _fb->getDecFrame();
  int* iDecodedInvQ    = _fb->getInvQuantDecFrame();

#if SI_REFINEMENT
  imgpel* imgRefinedSI[NCHANS];
  for (int c = 0; c < NCHANS; c++) {
    frameSize = (c == 0) ? Y_FSIZE : UV_FSIZE;
    imgRefinedSI[c]    = new imgpel[frameSize];
  }
#endif // SI_REFINEMENT

  int x,y;
  double totalrate[NCHANS] = {0., 0., 0.};
  double dMVrate=0;

  double dKeyCodingRate=0;
  double dKeyPSNR=0;

  FILE* fReadPtr       = _files->getFile("origin")->getFileHandle();
  FILE* fWritePtr      = _files->getFile("rec")->getFileHandle();
  FILE* fKeyReadPtr    = _files->getFile("key")->getFileHandle();
  FILE* mvFilePtr      = _files->getFile("mv")->getFileHandle();

  parseKeyStat("stats.dat", dKeyCodingRate, dKeyPSNR, _keyQp);

  timeStart = clock();

  // Main loop
  // ---------------------------------------------------------------------------
  for (int keyFrameNo = 0; keyFrameNo < (_numFrames-1)/_gop; keyFrameNo++) {
    // Read previous key frame
    fseek(fKeyReadPtr, (keyFrameNo)*FSIZE, SEEK_SET);
    fread(_fb->getPrevFrame(), FSIZE, 1, fKeyReadPtr);

    // Read next key frame
    fseek(fKeyReadPtr, (keyFrameNo+1)*FSIZE, SEEK_SET);
    fread(_fb->getNextFrame(), FSIZE, 1, fKeyReadPtr);

    for (int il = 0; il < _gopLevel; il++) {
      int frameStep = _gop / ((il+1)<<1);
      int idx = frameStep;

      // Start decoding the WZ frame
      while (idx < _gop) {
        int wzFrameNo = keyFrameNo*_gop + idx;

        cout << "Decoding frame " << wzFrameNo << " (Wyner-Ziv frame)" << endl;

        // Read current frame from the original file
        fseek(fReadPtr, wzFrameNo * FSIZE, SEEK_SET);
        fread(oriCurrFrame, FSIZE, 1, fReadPtr);

        // Setup frame pointers within the GOP
        int prevIdx = idx - frameStep;
        int nextIdx = idx + frameStep;

        imgpel* currFrame = _fb->getRecFrames()[idx-1];
        imgpel* prevFrame = (prevIdx == 0)    ? _fb->getPrevFrame() :
                                                _fb->getRecFrames()[prevIdx-1];
        imgpel* nextFrame = (nextIdx == _gop) ? _fb->getNextFrame() :
                                                _fb->getRecFrames()[nextIdx-1];
        imgpel* prevKeyFrame  = _fb->getPrevFrame();
        imgpel* nextKeyFrame  = _fb->getNextFrame();

        // ---------------------------------------------------------------------
        // Perform the following stages for each channel
        // ---------------------------------------------------------------------
        for (int c = 0; c < NCHANS; c++) {
          int  frameSize = (c == 0) ? Y_FSIZE : UV_FSIZE;
          // ---------------------------------------------------------------------
          // STAGE 1 - Create side information
          // ---------------------------------------------------------------------
          _si->createSideInfo(CHOFFSET(prevFrame,c), CHOFFSET(nextFrame,c),
                             CHOFFSET(imgSI,c), mvFilePtr, c);

          // ---------------------------------------------------------------------
          // STAGE 2 -
          // ---------------------------------------------------------------------
          int tmp = getSyndromeData(c);

          double dTotalRate = (double)tmp/1024/8;

          _trans->dctTransform(CHOFFSET(imgSI,c), CHOFFSET(iDCT,c), c);

          memset(CHOFFSET(iDecoded,c), 0, frameSize*sizeof(int));
          memset(CHOFFSET(iDecodedInvQ,c), 0, frameSize*sizeof(int));

# if RESIDUAL_CODING
          _si->getResidualFrame(CHOFFSET(prevKeyFrame,c),
                                CHOFFSET(nextKeyFrame,c),
                                CHOFFSET(imgSI,c),
                                CHOFFSET(iDCTBuffer,c), _rcList[c], c);

          _trans->dctTransform(CHOFFSET(iDCTBuffer,c),
                               CHOFFSET(iDCTResidual,c), c);
          _trans->quantization(CHOFFSET(iDCTResidual,c),
                               CHOFFSET(iDCTQ,c), c);

          int iOffset = 0;
          int iDC;

#   if SI_REFINEMENT
          memcpy(CHOFFSET(iDecodedInvQ,c),
                 CHOFFSET(iDCTResidual,c), frameSize*sizeof(int));
#   endif

          for (int i = 0; i < 16; i++) {
            x = ScanOrder[i][0];
            y = ScanOrder[i][1];

#   if MODE_DECISION
            if (i < _numChnCodeBands[c])
              dTotalRate += decodeLDPC(CHOFFSET(iDCTQ,c),
                                       CHOFFSET(iDCTResidual,c),
                                       CHOFFSET(iDecoded,c),
                                       x, y, iOffset, c);
#   else
            dTotalRate += decodeLDPC(CHOFFSET(iDCTQ,c),
                                     CHOFFSET(iDCTResidual,c),
                                     CHOFFSET(iDecoded,c),
                                     x, y, iOffset, c);
#   endif

#   if SI_REFINEMENT
            //temporal reconstruction
            _trans->invQuantization(CHOFFSET(iDecoded,c),
                                    CHOFFSET(iDecodedInvQ,c),
                                    CHOFFSET(iDCTResidual,c), x, y, c);
            _trans->invDctTransform(CHOFFSET(iDecodedInvQ,c),
                                    CHOFFSET(iDCTBuffer,c), c);

            _si->getRecFrame(CHOFFSET(prevKeyFrame,c),
                             CHOFFSET(nextKeyFrame,c),
                             CHOFFSET(iDCTBuffer,c),
                             CHOFFSET(currFrame,c), _rcList[c], c);

            iDC = (x == 0 && y == 0) ? 0 : 1;

            _si->getRefinedSideInfo(CHOFFSET(prevFrame,c),
                                    CHOFFSET(nextFrame,c),
                                    CHOFFSET(imgSI,c),
                                    CHOFFSET(currFrame,c),
                                    imgRefinedSI[c], iDC, c);

            memcpy(CHOFFSET(imgSI,c),
                   CHOFFSET(imgRefinedSI,c),
                   frameSize*sizof(int));

            _si->getResidualFrame(CHOFFSET(prevKeyFrame,c),
                                  CHOFFSET(nextFrame,c),
                                  CHOFFSET(imgSI,c),
                                  CHOFFSET(iDCTBuffer,c), _rcList[c], c);

            _trans->dctTransform(CHOFFSET(iDCTBuffer,c),
                                 CHOFFSET(iDCTResidual,c), c);
            _trans->quantization(CHOFFSET(iDCTResidual,c),
                                 CHOFFSET(iDCTQ,c), c);
#   endif

            iOffset += QuantMatrix[_qp][y][x];
          }

#   if !SI_REFINEMENT
          _trans->invQuantization(CHOFFSET(iDecoded,c),
                                  CHOFFSET(iDecodedInvQ,c),
                                  CHOFFSET(iDCTResidual,c), c);
          _trans->invDctTransform(CHOFFSET(iDecodedInvQ,c),
                                  CHOFFSET(iDCTBuffer,c), c);

          _si->getRecFrame(CHOFFSET(prevFrame,c),
                           CHOFFSET(nextFrame,c),
                           CHOFFSET(iDCTBuffer,c),
                           CHOFFSET(currFrame,c), _rcList[c], c);
#   endif

# else // if !RESIDUAL_CODING

          _trans->quantization(CHOFFSET(iDCT,c), CHOFFSET(iDCTQ,c), c);

          int iOffset = 0;
          int iDC;

#   if SI_REFINEMENT
          memcpy(CHOFFSET(iDecodedInvQ,c),
                 CHOFFSET(iDCT,c),
                 frameSize*sizeof(int));
#   endif

          for (int i = 0; i < 16; i++) {
            x = ScanOrder[i][0];
            y = ScanOrder[i][1];

#   if MODE_DECISION
            if (i < _numChnCodeBands)
              dTotalRate += decodeLDPC(CHOFFSET(iDCTQ,c), CHOFFSET(iDCT,c),
                                       CHOFFSET(iDecoded,c), x, y, iOffset, c);
#   else
            dTotalRate += decodeLDPC(CHOFFSET(iDCTQ,c), CHOFFSET(iDCT,c),
                                     CHOFFSET(iDecoded,c), x, y, iOffset, c);
#   endif

#   if SI_REFINEMENT
            _trans->invQuantization(CHOFFSET(iDecoded,c),
                                    CHOFFSET(iDecodedInvQ,c),
                                    CHOFFSET(iDCT,c), x, y, c);
            _trans->invDctTransform(CHOFFSET(iDecodedInvQ,c),
                                    CHOFFSET(currFrame,c), c);

#     if SKIP_MODE
            //reconstruct skipped part of wyner-ziv frame
            getSkippedRecFrame(CHOFFSET(prevKeyFrame,c),
                               CHOFFSET(currFrame,c), _skipMask[c]);
#     endif
            iDC = (x == 0 && y == 0) ? 0 : 1;

            _si->getRefinedSideInfo(CHOFFSET(prevFrame,c),
                                    CHOFFSET(nextFrame,c),
                                    CHOFFSET(imgSI,c),
                                    CHOFFSET(currFrame,c),
                                    CHOFFSET(imgRefinedSI,c), iDC);

            memcpy(CHOFFSET(imgSI,c),
                   CHOFFSET(imgRefinedSI,c),
                   frameSize*sizeof(int));

            _trans->dctTransform(CHOFFSET(imgSI,c), CHOFFSET(iDCT,c),c);
            _trans->quantization(CHOFFSET(iDCT,c), CHOFFSET(iDCTQ,c),c);
#   endif
            iOffset += QuantMatrix[_qp][y][x];
          }

#   if !SI_REFINEMENT
          _trans->invQuantization(CHOFFSET(iDecoded,c), 
                                  CHOFFSET(iDecodedInvQ,c), 
                                  CHOFFSET(iDCT,c), c);
          _trans->invDctTransform(CHOFFSET(iDecodedInvQ,c),
                                  CHOFFSET(currFrame,c), c);

#     if SKIP_MODE
          getSkippedRecFrame(CHOFFSET(prevKeyFrame,c),
                             CHOFFSET(currFrame,c), _skipMask[c]);
#     endif

#   endif

# endif // RESIDUAL_CODING
        totalrate[c] += dTotalRate;

        }
        cout << endl;
        //cout << "total bits (Y/frame): " << dTotalRate << " Kbytes" << endl;

        //cout << "side information quality" << endl;
        dPSNRSIAvg += calcPSNR(oriCurrFrame, imgSI, FSIZE);

        //cout << "wyner-ziv frame quality" << endl;
        dPSNRAvg += calcPSNR(oriCurrFrame, currFrame, FSIZE);

        if (wzFrameNo == 3) {
          char fname1[] = "sideinfo";  
          char fname2[] = "invquant";  
          FILE* f1, *f2;
          f1 = fopen(fname1, "wb");
          fwrite(imgSI, 1, FSIZE, f1);
          fclose(f1);
          f2 = fopen(fname2, "wb");
          fwrite(currFrame, 1, FSIZE, f2);
          fclose(f2);
        }

        idx += 2*frameStep;
      }
    }


    // ---------------------------------------------------------------------
    // Output decoded frames of the whole GOP
    // ---------------------------------------------------------------------

    // First output the key frame
    fwrite(_fb->getPrevFrame(), FSIZE, 1, fWritePtr);

    // Then output the rest WZ frames
    for (int i = 0; i < _gop-1; i++) {
      fwrite(_fb->getRecFrames()[i], FSIZE, 1, fWritePtr);
    }
  }

  timeEnd = clock();
  cpuTime = (timeEnd - timeStart) / CLOCKS_PER_SEC;

  int iDecodeWZFrames = ((_numFrames-1)/_gop)*(_gop-1);
  int iNumGOP = (_numFrames-1)/_gop;
  int iTotalFrames = iDecodeWZFrames + iNumGOP;

  cout<<endl;
  cout<<"--------------------------------------------------"<<endl;
  cout<<"Decode statistics"<<endl;
  cout<<"--------------------------------------------------"<<endl;
  cout<<"Total Frames        :   "<<iTotalFrames<<endl;
  float framerate = 30.0;
  dPSNRAvg   /= iDecodeWZFrames;
  dPSNRSIAvg /= iDecodeWZFrames;
  cout<<"Total Bytes         :   "<<totalrate<<endl;
  // TODO: double check these rate calculations
  int sumrate = totalrate[0] + totalrate[1] + totalrate[2];
  cout<<"WZ Avg Rate  (kbps) :   "<<sumrate/double(iDecodeWZFrames)*framerate*(iDecodeWZFrames)/(double)iTotalFrames*8.0<<endl;
  cout<<"Key Avg Rate (kbps) :   "<<dKeyCodingRate*framerate*(iNumGOP)/(double)iTotalFrames<<endl;
  cout<<"Avg Rate (Key+WZ)   :   "<<sumrate/double(iDecodeWZFrames)*framerate*(iDecodeWZFrames)/(double)iTotalFrames*8.0+dKeyCodingRate*framerate*(iNumGOP)/(double)iTotalFrames<<endl;
  cout<<"Key Frame Quality   :   "<<dKeyPSNR<<endl;
  cout<<"SI Avg PSNR         :   "<<dPSNRSIAvg<<endl;
  cout<<"WZ Avg PSNR         :   "<<dPSNRAvg<<endl;
  cout<<"Avg    PSNR         :   "<<(dPSNRAvg+dKeyPSNR)/2<<endl;
  cout<<"Total Decoding Time :   "<<cpuTime<<"(s)"<<endl;
  cout<<"Avg Decoding Time   :   "<<cpuTime/(iDecodeWZFrames)<<endl;
  cout<<"--------------------------------------------------"<<endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::parseKeyStat(const char* filename, double &rate, double &psnr, int &QP)
{
  ifstream stats(filename, ios::in);

  if (!stats.is_open())
    return;

  char buf[1024];
  double iSlice_chroma = 0.0;
  double iSliceRate = 0.0;

  RegExp* rgx = RegExp::getInst();

  while (stats.getline(buf, 1024)) {
    string result;

    if (rgx->match(buf, "\\s*SNR Y\\(dB\\)[ |]*([0-9\\.]+)", result)) {
      psnr = atof(result.c_str());
      continue;
    }

    if (rgx->match(buf, "\\s*Average quant[ |]*([0-9\\.]+)", result)) {
      QP = atoi(result.c_str());
      continue;
    }

    if (rgx->match(buf, "\\s*QP[ |]*([0-9\\.]+)", result)) {
      QP = atoi(result.c_str());
      continue;
    }

    if (rgx->match(buf, "\\s*Coeffs\\. C[ |]*([0-9\\.]+)", result)) {
      iSlice_chroma = atof(result.c_str());
      continue;
    }

    if (rgx->match(buf, "\\s*average bits/frame[ |]*([0-9\\.]+)", result)) {
      iSliceRate = atof(result.c_str());
      continue;
    }
  }

  int count = 0;
  double totalRate = 0;

  if (iSliceRate != 0) {
    totalRate += iSliceRate - iSlice_chroma;
    count++;
  }

  if (count)
    rate = totalRate/(1024*count);
  else
    rate = 0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int Decoder::getSyndromeData(int c)
{
  int* iDecoded = _fb->getDecFrame();
  int  decodedBits = 0;
  int  frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
  int  frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;
  int  frameSize = (c == 0) ? Y_FSIZE : UV_FSIZE;
  int  bplen = (c == 0) ? Y_BPLEN : UV_BPLEN;

# if RESIDUAL_CODING

#   if !HARDWARE_FLOW
  // Decode motion vector
  for (int i = 0; i < frameSize/64; i++)
    _rcList[c][i] = _bs[c]->read(1);
#   endif

# endif // RESIDUAL_CODING

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 4; i++) {
# if RESIDUAL_CODING

#   if !HARDWARE_FLOW
      _maxValue[c][j][i] = _bs[c]->read(11);
#   endif

      if (QuantMatrix[_qp][j][i] != 0) {
#   if HARDWARE_QUANTIZATION
        _quantStep[c][j][i] = 1 << (MaxBitPlane[j][i]+1-QuantMatrix[_qp][j][i]);
#   else
        int iInterval = 1 << QuantMatrix[_qp][j][i];

        _quantStep[c][j][i] =
          (int)(ceil(double(2*abs(_maxValue[c][j][i]))/double(iInterval-1)));
        _quantStep[c][j][i] = Max(_quantStep[c][j][i], MinQStepSize[_qp][j][i]);
#   endif
      }
      else
        _quantStep[c][j][i] = 1;

# else // if !RESIDUAL_CODING

      if (i != 0 || j != 0) {
        _maxValue[c][j][i] = _bs[c]->read(11);

#   if HARDWARE_QUANTIZATION
        if (QuantMatrix[_qp][j][i] != 0)
          _quantStep[c][j][i] = 1 << (MaxBitPlane[j][i]+1-QuantMatrix[_qp][j][i]);
        else
          _quantStep[c][j][i] = 0;
#   else
        int iInterval = 1 << QuantMatrix[_qp][j][i];

#     if AC_QSTEP
        if (QuantMatrix[_qp][j][i] != 0) {
          _quantStep[c][j][i] = (int)(ceil(double(2*abs(_maxValue[c][j][i]))/double(iInterval-1)));

          if (_quantStep[j][i] < 0)
            _quantStep[c][j][i] = 0;
        }
        else
          _quantStep[c][j][i] = 1;
#     else
        if (QuantMatrix[_qp][j][i] != 0)
          _quantStep[c][j][i] = ceil(double(2*abs(_maxValue[c][j][i]))/double(iInterval));
        else
          _quantStep[c][j][i] = 1;
#     endif

#   endif // HARDWARE_QUANTIZATION
      }
      else {
        _maxValue[c][j][i] = DC_BITDEPTH;

        _quantStep[c][j][i] = 1 << (_maxValue[j][i]-QuantMatrix[_qp][j][i]);
      }
# endif // RESIDUAL_CODING
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
# if MODE_DECISION
  int codingMode = _bs[c]->read(2);

  if (codingMode == 0) _numChnCodeBands[c] = 16; else
  if (codingMode == 1) _numChnCodeBands[c] =  3; else
  if (codingMode == 2) _numChnCodeBands[c] =  6; else
                       _numChnCodeBands[c] =  0;

  int numBands = 0;

  _rcBitPlaneNum = 0;

  for (int bandNo = 0; bandNo < 16; bandNo++) {
    int x = ScanOrder[bandNo][0];
    int y = ScanOrder[bandNo][1];

    if (bandNo < _numChnCodeBands[c]) {
      _rcQuantMatrix[y][x] = QuantMatrix[_qp][y][x];
      _rcBitPlaneNum += _rcQuantMatrix[y][x];
    }

    if (QuantMatrix[_qp][y][x] > 0)
      numBands++;
  }

# else // if !MODE_DECISION

  _rcBitPlaneNum = 0;

  for (int j = 0; j < 4; j++)
    for (int i = 0; i < 4; i++) {
      _rcQuantMatrix[j][i] = QuantMatrix[_qp][j][i];
      _rcBitPlaneNum += _rcQuantMatrix[j][i];
    }
# endif // MODE_DECISION

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
# if HARDWARE_FLOW
  // Discard padded bits
  int dummy = 20;
  _bs[c]->read(dummy);

  int bitCount = _bs[c]->getBitCount();
# endif

# if SKIP_MODE
  decodedBits += decodeSkipMask(c);
# endif

# if HARDWARE_FLOW
  bitCount = _bs[c]->getBitCount() - bitCount;

  if (bitCount%32 != 0) {
    dummy = 32 - (bitCount%32);
    _bs[c]->read(dummy);
  }

  bitCount = _bs[c]->getBitCount();
# endif

# if MODE_DECISION
  if (numBands > _numChnCodeBands[c]) {
    for (int j = 0; j < frameHeight; j += 4)
      for (int i = 0; i < frameWidth; i += 4) {
        if (_skipMask[c][i/4+(j/4)*(frameWidth/4)] == 0) //not skip
          decodedBits += _cavlc->decode(iDecoded, i, j, c);
        else
          _cavlc->clearNnz(i/4+(j/4)*(frameWidth/4), c);
      }
  }
# endif

# if HARDWARE_FLOW
  bitCount = _bs[c]->getBitCount() - bitCount;

  if (bitCount%32 != 0) {
    dummy = 32 - (bitCount%32);
    _bs[c]->read(dummy);
  }
# endif

  //TODO: FIX THIS!!
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Read parity and CRC bits from the bitstream
# if RESIDUAL_CODING
  for (int i = 0; i < _rcBitPlaneNum; i++)
# else
  for (int i = 0; i < BitPlaneNum[_qp]; i++)
# endif
  {
    for (int j = 0; j < bplen; j++)
      _dParity[c][j+i*bplen] = (double)_bs[c]->read(1);

# if !HARDWARE_FLOW
#   if HARDWARE_LDPC
    if (c == 0)
      for (int n = 0; n < 4; n++)
        _crc[i*4+n] = _bs[c]->read(8);
    else
      _crc[c][i] = _bs[c]->read(8);
#   else
    _crc[c][i] = _bs[c]->read(8);
#   endif
# endif
  }

  return decodedBits;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int Decoder::decodeSkipMask(int c)
{
  int type   = 0;
  int sign   = 0;
  int index  = 0;
  int length = 0;
  int  bplen = (c == 0) ? Y_BPLEN : UV_BPLEN;

  type = _bs[c]->read(2);
  sign = _bs[c]->read(1);

  memset(_skipMask[c], 0, bplen*sizeof(int));

  while (index < bplen) {
    int code = 0;
    int run  = 0;

    for (length = 1; length < 15; length++) {
      code <<= 1;
      code |= _bs[c]->read(1);

      for (int i = 0; i < 16; i++) {
        int table = _qp / 2;

        if (HuffmanCodeValue [table][type][i] == code &&
            HuffmanCodeLength[table][type][i] == length) {
          run = i+1;
          goto DecodeSkipMaskHuffmanCodeDone;
        }
      }
    }

    DecodeSkipMaskHuffmanCodeDone:

    // Reconstruct skip mask
    if (run == 16)
      for (int i = 0; i < 15; i++)
        _skipMask[c][index++] = sign;
    else {
      for (int i = 0; i < run; i++)
        _skipMask[c][index++] = sign;

      sign = (sign == 0) ? 1 : 0;
    }
  }

  return length;
}

/*
*Decoding Process of LDPC
*Param
*/
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
double Decoder::decodeLDPC(int* iQuantDCT, int* iDCT, int* iDecoded,
                           int x, int y, int iOffset, int c)
{
  int iCurrPos;
  int* iDecodedTmp;
  double* dLLR;
  double* dAccumulatedSyndrome;
  double* dLDPCDecoded;
  double* dSource;
  double dRate,dTotalRate;
  double dErr;
  double dParityRate;
  unsigned char* ucCRCCode;
  int    iNumCode;
  int  bplen = (c == 0) ? Y_BPLEN : UV_BPLEN;
  int  frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
  int  frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;

# if HARDWARE_LDPC
  ucCRCCode             = _crc[c] + iOffset*4;
# else
  ucCRCCode             = _crc + iOffset;
# endif
  dAccumulatedSyndrome  = _dParity[c] + bplen*iOffset;

  dLLR         = new double[bplen];
  iDecodedTmp  = new int   [bplen];
  dLDPCDecoded = new double[bplen];
  dSource      = new double[bplen];
  dTotalRate   = 0;

  memset(iDecodedTmp, 0, bplen*sizeof(int));

  dParityRate = 0;

# if RESIDUAL_CODING
  for (iCurrPos = _rcQuantMatrix[y][x]-1; iCurrPos >= 0; iCurrPos--)
# else
  for (iCurrPos = QuantMatrix[_qp][y][x]-1; iCurrPos >= 0; iCurrPos--)
# endif
  {
# if RESIDUAL_CODING
    if (iCurrPos == _rcQuantMatrix[y][x]-1)
      dParityRate = _model->getSoftInput(iQuantDCT, _skipMask[c],
                                         iCurrPos, iDecodedTmp,
                                         dLLR, x, y, 1, c);
    else
      dParityRate = _model->getSoftInput(iDCT, _skipMask[c],
                                         iCurrPos, iDecodedTmp,
                                         dLLR, x, y, 2, c);
# else
    if (x == 0 && y == 0)
      dParityRate = _model->getSoftInput(iDCT, _skipMask[c],
                                         iCurrPos, iDecodedTmp,
                                         dLLR, 0, 0, 2);
    else {
      if (iCurrPos == QuantMatrix[_qp][y][x]-1)
        dParityRate = _model->getSoftInput(iQuantDCT, _skipMask[c],
                                           iCurrPos, iDecodedTmp,
                                           dLLR, x, y, 1);
      else
        dParityRate = _model->getSoftInput(iDCT, _skipMask[c],
                                           iCurrPos, iDecodedTmp,
                                           dLLR, x, y, 2);
    }
# endif
    iNumCode = int(dParityRate*66);

    if (iNumCode <= 2)
      iNumCode = 2;
    if (iNumCode >= 66)
      iNumCode = 66;
    iNumCode = 2;

# if HARDWARE_LDPC
    if (bplen == 6336) {
      double dRateTmp = 0;

      for (int n = 0; n < 4; n++) {
        iNumCode = 2;

        _ldpca->decode(dLLR+n*1584, dAccumulatedSyndrome+n*1584,
                       dSource+n*1584, dLDPCDecoded+n*1584,
                       &dRate, &dErr, *(ucCRCCode+n), iNumCode);
        dRateTmp += (dRate/4.0);
        //cout<<dRate<<endl;
        dRate = 0;
      }
      ucCRCCode += 4;
      cout << ".";

      dTotalRate += dRateTmp;
      dRate = 0;
    }
    else {
      _ldpca->decode(dLLR, dAccumulatedSyndrome, dSource, dLDPCDecoded,
                     &dRate, &dErr, *ucCRCCode, iNumCode);
      cout << ".";

      dTotalRate += dRate;
      dRate = 0;
      ucCRCCode++;
    }
# else
    _ldpca->decode(dLLR, dAccumulatedSyndrome,
                   dSource, dLDPCDecoded, &dRate,
                   &dErr, *ucCRCCode, iNumCode);
    cout << ".";

    dTotalRate += dRate;
    dRate = 0;
    ucCRCCode++;
# endif

    for (int iIndex = 0; iIndex < bplen; iIndex++)
      if (dLDPCDecoded[iIndex] == 1)
        iDecodedTmp[iIndex] |= 0x1<<iCurrPos;

    dAccumulatedSyndrome += bplen;

    memset(dLDPCDecoded, 0, bplen*sizeof(double));
  }

  for (int j = 0; j < frameHeight; j = j+4)
    for (int i = 0; i < frameWidth; i = i+4) {
      int tmp = i/4 + j/4*(frameWidth/4);
      iDecoded[(i+x)+(j+y)*frameWidth] = iDecodedTmp[tmp];
    }

  delete [] dLLR;
  delete [] iDecodedTmp;
  delete [] dLDPCDecoded;
  delete [] dSource;

  return (dTotalRate*bplen/8/1024);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::getSourceBit(int *dct_q,double *source,
                           int q_i,int q_j,int curr_pos, int c){
  int  frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
  int  frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;
  for(int y=0;y<frameHeight;y=y+4)
    for(int x=0;x<frameWidth;x=x+4)
    {
      source[(x/4)+(y/4)*(frameWidth/4)] =
        (dct_q[(x+q_i)+(y+q_j)*frameWidth]>>curr_pos)&(0x1);
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::motionSearchInit(int maxsearch_range)
{
  _spiralHpelSearchX = new int[(2*maxsearch_range+1)*(2*maxsearch_range+1)];
  _spiralSearchX     = new int[(2*maxsearch_range+1)*(2*maxsearch_range+1)];
  _spiralHpelSearchY = new int[(2*maxsearch_range+1)*(2*maxsearch_range+1)];
  _spiralSearchY     = new int[(2*maxsearch_range+1)*(2*maxsearch_range+1)];

  int k,i,l;

  _spiralSearchX[0] = _spiralSearchY[0] = 0;
  _spiralHpelSearchX[0] = _spiralHpelSearchY[0] = 0;

  for (k=1, l=1; l <= std::max<int>(1,maxsearch_range); l++) {
    for (i=-l+1; i< l; i++) {
      _spiralSearchX[k] =  i;
      _spiralSearchY[k] = -l;
      _spiralHpelSearchX[k] =  i<<1;
      _spiralHpelSearchY[k++] = -l<<1;
      _spiralSearchX[k] =  i;
      _spiralSearchY[k] =  l;
      _spiralHpelSearchX[k] =  i<<1;
      _spiralHpelSearchY[k++] =  l<<1;
    }
    for (i=-l;   i<=l; i++) {
      _spiralSearchX[k] = -l;
      _spiralSearchY[k] =  i;
      _spiralHpelSearchX[k] = -l<<1;
      _spiralHpelSearchY[k++] = i<<1;
      _spiralSearchX[k] =  l;
      _spiralSearchY[k] =  i;
      _spiralHpelSearchX[k] =  l<<1;
      _spiralHpelSearchY[k++] = i<<1;
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
double calcPSNR(unsigned char* img1,unsigned char* img2,int length)
{
  float PSNR;
  float MSE=0;

  for(int i=0;i<length;i++)
    {
      MSE+=pow(float(img1[i]-img2[i]),float(2.0))/length;
    }
  PSNR=10*log10(255*255/MSE);
  //cout<<"PSNR: "<<PSNR<<" dB"<<endl;
  return PSNR;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int getSymbol(int len, int& curr_pos, char* buffer)
{
  int temp = 0;

  for (int count = 0; count < len; count++) {
    int pos = count + curr_pos;

    temp <<= 1;
    temp |= 0x1 & (buffer[pos/8]>>(7-(pos%8)));
  }

  curr_pos += len;

  return temp;
}

