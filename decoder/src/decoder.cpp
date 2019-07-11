
#include <sstream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <exception>

#include "decoder.h"
#include "fileManager.h"
#include "sideInformation.h"
#include "transform.h"
#include "corrModel.h"
#include "time.h"
#include "cavlcDec.h"
#include "frameBuffer.h"
#include "bitstream.h"
#include "ldpcaDec.h"
#include "regExp.h"

using namespace std;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
Decoder::Decoder(char **argv)
{
  _files = FileManager::getManager();

  string wzFileName = argv[1];
  string recFileName = wzFileName.substr(0, wzFileName.find(".bin"));

  _files->addFile("wz",     argv[1])->openFile("rb");
  _files->addFile("key",    argv[2])->openFile("rb");
  _files->addFile("origin", argv[3])->openFile("rb");
  _channel = argv[4][0];
  _SIMethod = (SIMethod)atoi(argv[5]);

  if (_SIMethod == SAME)
    if (_channel == 'y')
      _files->addFile("helper", argv[6])->openFile("w");
    else
      _files->addFile("helper", argv[6])->openFile("r");
  else if (_SIMethod == MC_RE)
    _files->addFile("helper", argv[6])->openFile("rb");
  

  _bs = new Bitstream(1024, _files->getFile("wz")->getFileHandle());
  _files->addFile("rec", recFileName.c_str())->openFile("wb");

  decodeWzHeader();

  initialize();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::initialize()
{
  _frameSize        = _frameWidth * _frameHeight;
  _bitPlaneLength   = _frameSize / 16;
  _gop              = 1 << _gopLevel;

  _numChnCodeBands  = 16;

  _dParity          = new double[_bitPlaneLength * BitPlaneNum[_qp]];

# if HARDWARE_LDPC
  if (_bitPlaneLength == 6336)
    _crc            = new unsigned char[BitPlaneNum[_qp] * 4];
  else
    _crc            = new unsigned char[BitPlaneNum[_qp]];
# else
  _crc              = new unsigned char[BitPlaneNum[_qp]];
# endif

  _average          = new double[16];
  _alpha            = new double[_frameSize];
  _sigma            = new double[16];

# if RESIDUAL_CODING
  _rcList           = new int[_frameSize/64];

  for (int i = 0; i < _frameSize/64; i++)
    _rcList[i] = 0;
# endif

  _skipMask         = new int[_bitPlaneLength];

  _fb = new FrameBuffer(_frameWidth, _frameHeight, _gop);

  _trans = new Transform(this);

  _model = new CorrModel(this, _trans);
  if (_SIMethod == MCI_RE)
    _si = new SI_MCI(this, _model);
  else if (_SIMethod == SAME)
    _si = new SI_MCI(this,_model, _files->getFile("helper")->getFileHandle());
  else if (_SIMethod == MC_RE)
    _si = new SI_MC(this, _files->getFile("helper")->getFileHandle(),
                    _frameHeight*2, _frameWidth*2);
  else
    throw invalid_argument("Unsupported SI generation method.");


  _cavlc = new CavlcDec(this, 4);

  motionSearchInit(64);

  // Initialize LDPC
  string ladderFile;

# if HARDWARE_LDPC
  ladderFile = "ldpca/1584_regDeg3.lad";
# else
  if (_frameWidth == 352 && _frameHeight == 288)
    ladderFile = "ldpca/6336_regDeg3.lad";
  else
    ladderFile = "ldpca/1584_regDeg3.lad";
# endif

  _ldpca = new LdpcaDec(ladderFile, this);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::decodeWzHeader()
{
  _frameWidth   = _bs->read(8) * 16;
  _frameHeight  = _bs->read(8) * 16;
  _qp           = _bs->read(8);
  _numFrames    = _bs->read(16);
  _gopLevel     = _bs->read(3);

  cout << "--------------------------------------------------" << endl;
  cout << "WZ frame parameters" << endl;
  cout << "--------------------------------------------------" << endl;
  cout << "Width:  " << _frameWidth << endl;
  cout << "Height: " << _frameHeight << endl;
  cout << "Frames: " << _numFrames << endl;
  cout << "QP:     " << _qp << endl;
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

  imgpel* oriCurrFrame = _fb->getCurrFrame();
  imgpel* imgSI      = _fb->getSideInfoFrame();
  imgpel* imgRefinedSI = new imgpel[_frameSize];
  imgpel* keyColour = new imgpel[_frameSize>>1];

  int* iDCT         = _fb->getDctFrame();
  int* iDCTQ        = _fb->getQuantDctFrame();
  int* iDecoded     = _fb->getDecFrame();
  int* iDecodedInvQ = _fb->getInvQuantDecFrame();

#if RESIDUAL_CODING
  int* iDCTBuffer   = new int [_frameSize];
  int* iDCTResidual = new int [_frameSize];
#endif

  int x,y;
  double totalrate=0;
  double dMVrate=0;

  double dKeyCodingRate=0;
  double dKeyPSNR=0;

  FILE* fReadPtr    = _files->getFile("origin")->getFileHandle();
  FILE* fWritePtr   = _files->getFile("rec")->getFileHandle();
  FILE* fKeyReadPtr = _files->getFile("key")->getFileHandle();

  parseKeyStat("stats.dat", dKeyCodingRate, dKeyPSNR, _keyQp);

  timeStart = clock();

  // Main loop
  // ---------------------------------------------------------------------------
  for (int keyFrameNo = 0; keyFrameNo < (_numFrames-1)/_gop; keyFrameNo++) {
    // Read previous key frame
    fseek(fKeyReadPtr, keyFrameNo*FSIZE, SEEK_SET);
    if (_channel == 'y') {
      fread(_fb->getPrevFrame(), _frameSize, 1, fKeyReadPtr);
      fseek(fKeyReadPtr, _frameSize>>1, SEEK_CUR);
    } else if (_channel == 'u') {
      fseek(fKeyReadPtr, U_OFFSET, SEEK_CUR);
      fread(_fb->getPrevFrame(), _frameSize, 1, fKeyReadPtr);
      fseek(fKeyReadPtr, _frameSize*5, SEEK_CUR);
    } else if (_channel == 'v') {
      fseek(fKeyReadPtr, V_OFFSET, SEEK_CUR);
      fread(_fb->getPrevFrame(), _frameSize, 1, fKeyReadPtr);
      fseek(fKeyReadPtr, _frameSize*5, SEEK_CUR);
    }
    // read nextFrame from exactly one frame ahead, regardless of offset
    fread(_fb->getNextFrame(), _frameSize, 1, fKeyReadPtr);

    for (int il = 0; il < _gopLevel; il++) {
      int frameStep = _gop / (1<<(il+1));
      int idx = frameStep;
      int prevNo = keyFrameNo*_gop;

      // Start decoding the WZ frame
      while (idx < _gop) {
        int wzFrameNo = prevNo + idx;

        cout << "Decoding frame " << wzFrameNo << " (Wyner-Ziv frame)" << endl;

        // Read current frame from the original file
        fseek(fReadPtr, wzFrameNo*FSIZE, SEEK_SET);
        if (_channel == 'y') {
          fread(oriCurrFrame, _frameSize, 1, fReadPtr);
        } else if (_channel == 'u') {
          fseek(fReadPtr, U_OFFSET, SEEK_CUR);
          fread(oriCurrFrame, _frameSize, 1, fReadPtr);
        } else if (_channel == 'v') {
          fseek(fReadPtr, V_OFFSET, SEEK_CUR);
          fread(oriCurrFrame, _frameSize, 1, fReadPtr);
        }

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
        // STAGE 1 - Create side information
        // ---------------------------------------------------------------------
        if (_SIMethod == SAME || _SIMethod == MCI_RE) {
          _si->createSideInfo(prevFrame, nextFrame, imgSI);
        } else if (_SIMethod == MC_RE) {
          _si->createSideInfo(prevFrame, nextFrame, imgSI,
                              prevNo + prevIdx,
                              prevNo + nextIdx,
                              prevNo + idx);
        } else if (_SIMethod == MC_NEW) {
          _si->createSideInfo(prevFrame, nextFrame, imgSI,
                              prevNo + prevIdx,
                              prevNo + nextIdx,
                              prevNo + idx);
        }

        if (_channel == 'u' || _channel == 'v') {
          memcpy(currFrame, imgSI, _frameSize);
          idx += 2*frameStep;
          dPSNRSIAvg += calcPSNR(oriCurrFrame, imgSI, _frameSize);
          dPSNRAvg += calcPSNR(oriCurrFrame, currFrame, _frameSize);
          continue;
        }

        // ---------------------------------------------------------------------
        // STAGE 2 -
        // ---------------------------------------------------------------------
        int tmp = getSyndromeData();

        //cout << _numChnCodeBands << endl;

        double dTotalRate = (double)tmp/1024/8;

        _trans->dctTransform(imgSI, iDCT);

        memset(iDecoded, 0, _frameSize*4);
        memset(iDecodedInvQ, 0, _frameSize*4);

# if RESIDUAL_CODING
        _si->getResidualFrame(prevKeyFrame, nextKeyFrame, imgSI, iDCTBuffer, _rcList);

        _trans->dctTransform(iDCTBuffer, iDCTResidual);
        _trans->quantization(iDCTResidual, iDCTQ);

        int iOffset = 0;
        int iDC;

#   if SI_REFINEMENT
        memcpy(iDecodedInvQ, iDCTResidual, 4*_frameSize);
#   endif

        for (int i = 0; i < 16; i++) {
          x = ScanOrder[i][0];
          y = ScanOrder[i][1];

#   if MODE_DECISION
          if (i < _numChnCodeBands)
            dTotalRate += decodeLDPC(iDCTQ, iDCTResidual, iDecoded, x, y, iOffset);
#   else
          dTotalRate += decodeLDPC(iDCTQ, iDCTResidual, iDecoded, x, y, iOffset);
#   endif

#   if SI_REFINEMENT
          //temporal reconstruction
          _trans->invQuantization(iDecoded, iDecodedInvQ, iDCTResidual, x, y);
          _trans->invDctTransform(iDecodedInvQ, iDCTBuffer);

          _si->getRecFrame(prevKeyFrame, nextKeyFrame, iDCTBuffer, currFrame, _rcList);

          iDC = (x == 0 && y == 0) ? 0 : 1;

          _si->getRefinedSideInfo(prevFrame, nextFrame, imgSI, currFrame, imgRefinedSI, iDC);

          memcpy(imgSI, imgRefinedSI, _frameSize);

          _si->getResidualFrame(prevKeyFrame, nextFrame, imgSI, iDCTBuffer, _rcList);

          _trans->dctTransform(iDCTBuffer, iDCTResidual);
          _trans->quantization(iDCTResidual, iDCTQ);
#   endif

          iOffset += QuantMatrix[_qp][y][x];
        }

#   if !SI_REFINEMENT
        _trans->invQuantization(iDecoded, iDecodedInvQ, iDCTResidual);
        _trans->invDctTransform(iDecodedInvQ, iDCTBuffer);

        _si->getRecFrame(prevFrame, nextFrame, iDCTBuffer, currFrame, _rcList);
#   endif

# else // if !RESIDUAL_CODING

        _trans->quantization(iDCT, iDCTQ);

        int iOffset = 0;
        int iDC;

#   if SI_REFINEMENT
        memcpy(iDecodedInvQ, iDCT, 4*_frameSize);
#   endif

        for (int i = 0; i < 16; i++) {
          x = ScanOrder[i][0];
          y = ScanOrder[i][1];

#   if MODE_DECISION
          if (i < _numChnCodeBands)
            dTotalRate += decodeLDPC(iDCTQ, iDCT, iDecoded, x, y, iOffset);
#   else
          dTotalRate += decodeLDPC(iDCTQ, iDCT, iDecoded, x, y, iOffset);
#   endif

#   if SI_REFINEMENT
          _trans->invQuantization(iDecoded, iDecodedInvQ, iDCT, x, y);
          _trans->invDctTransform(iDecodedInvQ, currFrame);

#     if SKIP_MODE
          //reconstruct skipped part of wyner-ziv frame
          getSkippedRecFrame(prevKeyFrame, currFrame, _skipMask);
#     endif
          iDC = (x == 0 && y == 0) ? 0 : 1;

          _si->getRefinedSideInfo(prevFrame, nextFrame, imgSI, currFrame, imgRefinedSI, iDC);

          memcpy(imgSI, imgRefinedSI, _frameSize);

          _trans->dctTransform(imgSI, iDCT);
          _trans->quantization(iDCT, iDCTQ);
#   endif
          iOffset += QuantMatrix[_qp][y][x];
        }

#   if !SI_REFINEMENT
        _trans->invQuantization(iDecoded, iDecodedInvQ, iDCT);
        _trans->invDctTransform(iDecodedInvQ, currFrame);

#     if SKIP_MODE
        getSkippedRecFrame(prevKeyFrame, currFrame, _skipMask);
#     endif

#   endif

# endif // RESIDUAL_CODING

        totalrate += dTotalRate;
        //cout << endl;
        //cout << "total bits (Y/frame): " << dTotalRate << " Kbytes" << endl;

        //cout << "side information quality" << endl;
        dPSNRSIAvg += calcPSNR(oriCurrFrame, imgSI, _frameSize);

        //cout << "wyner-ziv frame quality" << endl;
        dPSNRAvg += calcPSNR(oriCurrFrame, currFrame, _frameSize);

        idx += 2*frameStep;
      }
    }

    // ---------------------------------------------------------------------
    // Output decoded frames of the whole GOP
    // ---------------------------------------------------------------------

    // First output the key frame
    fwrite(_fb->getPrevFrame(), _frameSize, 1, fWritePtr);

    // Then output the rest WZ frames
    for (int i = 0; i < _gop-1; i++)
      fwrite(_fb->getRecFrames()[i], _frameSize, 1, fWritePtr);
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
  cout<<"WZ Avg Rate  (kbps) :   "<<totalrate/double(iDecodeWZFrames)*framerate*(iDecodeWZFrames)/(double)iTotalFrames*8.0<<endl;
  cout<<"Key Avg Rate (kbps) :   "<<dKeyCodingRate*framerate*(iNumGOP)/(double)iTotalFrames<<endl;
  cout<<"Avg Rate (Key+WZ)   :   "<<totalrate/double(iDecodeWZFrames)*framerate*(iDecodeWZFrames)/(double)iTotalFrames*8.0+dKeyCodingRate*framerate*(iNumGOP)/(double)iTotalFrames<<endl;
  cout<<"Key Frame Quality   :   "<<dKeyPSNR<<endl;
  cout<<"SI Avg PSNR         :   "<<dPSNRSIAvg<<endl;
  cout<<"WZ Avg PSNR         :   "<<dPSNRAvg<<endl;
  cout<<"Avg    PSNR         :   "<<((_gop-1)*dPSNRAvg+dKeyPSNR)/_gop<<endl;
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

    if ((_channel == 'y') &&
        (rgx->match(buf, "\\s*SNR Y\\(dB\\)[ |]*([0-9\\.]+)", result))) {
      psnr = atof(result.c_str());
      continue;
    } else if ((_channel == 'u') &&
        (rgx->match(buf, "\\s*SNR U\\(dB\\)[ |]*([0-9\\.]+)", result))) {
      psnr = atof(result.c_str());
      continue;
    } else if ((_channel == 'v') &&
        (rgx->match(buf, "\\s*SNR V\\(dB\\)[ |]*([0-9\\.]+)", result))) {
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

  if ((_channel == 'y') && (iSliceRate != 0)) {
    totalRate += iSliceRate - iSlice_chroma;
    count++;
  } else if ((_channel == 'u' || _channel == 'v') && (iSlice_chroma != 0)) {
    totalRate += iSlice_chroma / 2.0;
    count++;
  }

  if (count)
    rate = totalRate/(1024*count);
  else
    rate = 0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int Decoder::getSyndromeData()
{
  int* iDecoded = _fb->getDecFrame();
  int  decodedBits = 0;

# if RESIDUAL_CODING

#   if !HARDWARE_FLOW
  // Decode motion vector
  for (int i = 0; i < _frameSize/64; i++)
    _rcList[i] = _bs->read(1);
#   endif

# endif // RESIDUAL_CODING

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 4; i++) {
# if RESIDUAL_CODING

#   if !HARDWARE_FLOW
      _maxValue[j][i] = _bs->read(11);
#   endif

      if (QuantMatrix[_qp][j][i] != 0) {
#   if HARDWARE_QUANTIZATION
        _quantStep[j][i] = 1 << (MaxBitPlane[j][i]+1-QuantMatrix[_qp][j][i]);
#   else
        int iInterval = 1 << QuantMatrix[_qp][j][i];

        _quantStep[j][i] = (int)(ceil(double(2*abs(_maxValue[j][i]))/double(iInterval-1)));
        _quantStep[j][i] = Max(_quantStep[j][i], MinQStepSize[_qp][j][i]);
#   endif
      }
      else
        _quantStep[j][i] = 1;

# else // if !RESIDUAL_CODING

      if (i != 0 || j != 0) {
        _maxValue[j][i] = _bs->read(11);

#   if HARDWARE_QUANTIZATION
        if (QuantMatrix[_qp][j][i] != 0)
          _quantStep[j][i] = 1 << (MaxBitPlane[j][i]+1-QuantMatrix[_qp][j][i]);
        else
          _quantStep[j][i] = 0;
#   else
        int iInterval = 1 << QuantMatrix[_qp][j][i];

#     if AC_QSTEP
        if (QuantMatrix[_qp][j][i] != 0) {
          _quantStep[j][i] = (int)(ceil(double(2*abs(_maxValue[j][i]))/double(iInterval-1)));

          if (_quantStep[j][i] < 0)
            _quantStep[j][i] = 0;
        }
        else
          _quantStep[j][i] = 1;
#     else
        if (QuantMatrix[_qp][j][i] != 0)
          _quantStep[j][i] = ceil(double(2*abs(_maxValue[j][i]))/double(iInterval));
        else
          _quantStep[j][i] = 1;
#     endif

#   endif // HARDWARE_QUANTIZATION
      }
      else {
        _maxValue[j][i] = DC_BITDEPTH;

        _quantStep[j][i] = 1 << (_maxValue[j][i]-QuantMatrix[_qp][j][i]);
      }
# endif // RESIDUAL_CODING
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
# if MODE_DECISION
  int codingMode = _bs->read(2);

  if (codingMode == 0) _numChnCodeBands = 16; else
  if (codingMode == 1) _numChnCodeBands =  3; else
  if (codingMode == 2) _numChnCodeBands =  6; else
                       _numChnCodeBands =  0;

  int numBands = 0;

  _rcBitPlaneNum = 0;

  for (int bandNo = 0; bandNo < 16; bandNo++) {
    int x = ScanOrder[bandNo][0];
    int y = ScanOrder[bandNo][1];

    if (bandNo < _numChnCodeBands) {
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
  _bs->read(dummy);

  int bitCount = _bs->getBitCount();
# endif

# if SKIP_MODE
  decodedBits += decodeSkipMask();
# endif

# if HARDWARE_FLOW
  bitCount = _bs->getBitCount() - bitCount;

  if (bitCount%32 != 0) {
    dummy = 32 - (bitCount%32);
    _bs->read(dummy);
  }

  bitCount = _bs->getBitCount();
# endif

# if MODE_DECISION
  if (numBands > _numChnCodeBands) {
    for (int j = 0; j < _frameHeight; j += 4)
      for (int i = 0; i < _frameWidth; i += 4) {
        if (_skipMask[i/4+(j/4)*(_frameWidth/4)] == 0) //not skip
          decodedBits += _cavlc->decode(iDecoded, i, j);
        else
          _cavlc->clearNnz(i/4+(j/4)*(_frameWidth/4));
      }
  }
# endif

# if HARDWARE_FLOW
  bitCount = _bs->getBitCount() - bitCount;

  if (bitCount%32 != 0) {
    dummy = 32 - (bitCount%32);
    _bs->read(dummy);
  }
# endif

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Read parity and CRC bits from the bitstream
# if RESIDUAL_CODING
  for (int i = 0; i < _rcBitPlaneNum; i++)
# else
  for (int i = 0; i < BitPlaneNum[_qp]; i++)
# endif
  {
    for (int j = 0; j < _bitPlaneLength; j++)
      _dParity[j+i*_bitPlaneLength] = (double)_bs->read(1);

# if !HARDWARE_FLOW
#   if HARDWARE_LDPC
    if (_bitPlaneLength == 6336)
      for (int n = 0; n < 4; n++)
        _crc[i*4+n] = _bs->read(8);
    else
      _crc[i] = _bs->read(8);
#   else
    _crc[i] = _bs->read(8);
#   endif
# endif
  }

  return decodedBits;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int Decoder::decodeSkipMask()
{
  int type   = 0;
  int sign   = 0;
  int index  = 0;
  int length = 0;

  type = _bs->read(2);
  sign = _bs->read(1);

  memset(_skipMask, 0, _bitPlaneLength*sizeof(int));

  while (index < _bitPlaneLength) {
    int code = 0;
    int run  = 0;

    for (length = 1; length < 15; length++) {
      code <<= 1;
      code |= _bs->read(1);

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
        _skipMask[index++] = sign;
    else {
      for (int i = 0; i < run; i++)
        _skipMask[index++] = sign;

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
double Decoder::decodeLDPC(int* iQuantDCT, int* iDCT, int* iDecoded, int x, int y, int iOffset)
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

# if HARDWARE_LDPC
  ucCRCCode             = _crc + iOffset*4;
# else
  ucCRCCode             = _crc + iOffset;
# endif
  dAccumulatedSyndrome  = _dParity + _bitPlaneLength*iOffset;

  dLLR         = new double[_bitPlaneLength];
  iDecodedTmp  = new int   [_bitPlaneLength];
  dLDPCDecoded = new double[_bitPlaneLength];
  dSource      = new double[_bitPlaneLength];
  dTotalRate   = 0;

  memset(iDecodedTmp, 0, _bitPlaneLength*4);

  dParityRate = 0;

# if RESIDUAL_CODING
  for (iCurrPos = _rcQuantMatrix[y][x]-1; iCurrPos >= 0; iCurrPos--)
# else
  for (iCurrPos = QuantMatrix[_qp][y][x]-1; iCurrPos >= 0; iCurrPos--)
# endif
  {
# if RESIDUAL_CODING
    if (iCurrPos == _rcQuantMatrix[y][x]-1)
      dParityRate = _model->getSoftInput(iQuantDCT, _skipMask, iCurrPos, iDecodedTmp, dLLR, x, y, 1);
    else
      dParityRate = _model->getSoftInput(iDCT, _skipMask, iCurrPos, iDecodedTmp, dLLR, x, y, 2);
# else
    if (x == 0 && y == 0)
      dParityRate = _model->getSoftInput(iDCT, _skipMask, iCurrPos, iDecodedTmp, dLLR, 0, 0, 2);
    else {
      if (iCurrPos == QuantMatrix[_qp][y][x]-1)
        dParityRate = _model->getSoftInput(iQuantDCT, _skipMask, iCurrPos, iDecodedTmp, dLLR, x, y, 1);
      else
        dParityRate = _model->getSoftInput(iDCT, _skipMask, iCurrPos, iDecodedTmp, dLLR, x, y, 2);
    }
# endif
    iNumCode = int(dParityRate*66);

    if (iNumCode <= 2)
      iNumCode = 2;
    if (iNumCode >= 66)
      iNumCode = 66;
    iNumCode = 2;

# if HARDWARE_LDPC
    if (_bitPlaneLength == 6336) {
      double dRateTmp = 0;

      for (int n = 0; n < 4; n++) {
        iNumCode = 2;

        _ldpca->decode(dLLR+n*1584, dAccumulatedSyndrome+n*1584, dSource+n*1584, dLDPCDecoded+n*1584, &dRate, &dErr, *(ucCRCCode+n), iNumCode);
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
      _ldpca->decode(dLLR, dAccumulatedSyndrome, dSource, dLDPCDecoded, &dRate, &dErr, *ucCRCCode, iNumCode);
      cout << ".";

      dTotalRate += dRate;
      dRate = 0;
      ucCRCCode++;
    }
# else
    _ldpca->decode(dLLR, dAccumulatedSyndrome, dSource, dLDPCDecoded, &dRate, &dErr, *ucCRCCode, iNumCode);
    cout << ".";

    dTotalRate += dRate;
    dRate = 0;
    ucCRCCode++;
# endif

    for (int iIndex = 0; iIndex < _bitPlaneLength; iIndex++)
      if (dLDPCDecoded[iIndex] == 1)
        iDecodedTmp[iIndex] |= 0x1<<iCurrPos;

    dAccumulatedSyndrome += _bitPlaneLength;

    memset(dLDPCDecoded, 0, _bitPlaneLength*sizeof(double));
  }

  for (int j = 0; j < _frameHeight; j = j+4)
    for (int i = 0; i < _frameWidth; i = i+4) {
      int tmp = i/4 + j/4*(_frameWidth/4);
      iDecoded[(i+x)+(j+y)*_frameWidth] = iDecodedTmp[tmp];
    }

  delete [] dLLR;
  delete [] iDecodedTmp;
  delete [] dLDPCDecoded;
  delete [] dSource;

  return (dTotalRate*_bitPlaneLength/8/1024);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Decoder::getSourceBit(int *dct_q,double *source,int q_i,int q_j,int curr_pos){
  int iWidth,iHeight;
  iWidth  = _frameWidth;
  iHeight = _frameHeight;
  for(int y=0;y<iHeight;y=y+4)
    for(int x=0;x<iWidth;x=x+4)
    {
      source[(x/4)+(y/4)*(iWidth/4)]=(dct_q[(x+q_i)+(y+q_j)*iWidth]>>curr_pos)&(0x1);
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

