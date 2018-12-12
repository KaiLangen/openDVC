#include <iostream>
#include <sstream>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "encoder.h"
#include "transform.h"
#include "fileManager.h"
#include "frameBuffer.h"
#include "cavlcEnc.h"
#include "ldpcaEnc.h"
#include "bitstream.h"

using namespace std;

# if RESIDUAL_CODING
const int Encoder::Scale[3][8] = {
  {8, 6, 6, 4, 4, 3, 2, 1},
  {8, 8, 8, 4, 4, 4, 2, 1},
  {4, 4, 4, 4, 3, 2, 2, 1}
};
# endif // RESIDUAL_CODING

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
Encoder::Encoder(char** argv)
{
  _files = FileManager::getManager();

  _qp          = atoi(argv[1]);
  _keyQp       = atoi(argv[2]);
  _numFrames   = atoi(argv[3]);
  _gopLevel    = atoi(argv[4]);
  _files->addFile("src", argv[5])->openFile("rb");
  string wzFileName = argv[6];
  _files->addFile("key", argv[7]);


  char chan[NCHANS] = {'y','u','v'};
  for(int c = 0; c < NCHANS; c++) {
    string wzChanFile = wzFileName.substr(0, wzFileName.find(".bin"))
                         + "_" + chan[c] + ".bin";
    stringstream wzKey; 
    wzKey << "wz_" << chan[c];
    _files->addFile(wzKey.str(), wzChanFile)->openFile("wb");
    _bs[c] = new Bitstream(1024, _files->getFile(wzKey.str())->getFileHandle());
  }

  initialize();
}

Encoder::~Encoder()
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
    delete [] _alpha[c];
// TODO: fix this! _crc pointers are incremented, so unable to free. 
//                 Make an extra set of pointers, or use indices instead.
//                 delete _skipMask gets "double free or corruption" error.
//                 Investigate.
//    delete [] _crc[c];
//    delete [] _skipMask[c];
//    delete [] _parity[c];
    delete [] _average[c];
    delete _bs[c];
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::initialize()
{

  _gop                   = 1 << _gopLevel;

  _sigma                 = new double[NBANDS];
  _fb                    = new FrameBuffer();
  _trans                 = new Transform(this);
  _cavlc                 = new CavlcEnc(this, 4);

  _parity[0]             = new bool[Y_BPLEN * BitPlaneNum[_qp]];
  _parity[1]             = new bool[UV_BPLEN * BitPlaneNum[_qp]];
  _parity[2]             = new bool[UV_BPLEN * BitPlaneNum[_qp]];
  _alpha[0]              = new double[Y_FSIZE];
  _alpha[1]              = new double[UV_FSIZE];
  _alpha[2]              = new double[UV_FSIZE];

# if HARDWARE_LDPC
  _crc[0]                = new unsigned char[BitPlaneNum[_qp] * 4];
# else
  _crc[0]                = new unsigned char[BitPlaneNum[_qp]];
#endif
  _crc[1]                = new unsigned char[BitPlaneNum[_qp]];
  _crc[2]                = new unsigned char[BitPlaneNum[_qp]];

  _skipMask[0]           = new int[Y_BPLEN];
  _skipMask[1]           = new int[UV_BPLEN];
  _skipMask[2]           = new int[UV_BPLEN];

  for (int c = 0; c < NCHANS; c++) {
    _average[c]          = new double[NBANDS];
    _prevType[c]         = 0;
    _numChnCodeBands[c]  = 16;
    for (int i = 0; i < 4; i++)
      _modeCounter[c][i] = 0;
  }

  // Initialize LDPC
# if !HARDWARE_LDPC
  _ldpca_cif             = new LdpcaEnc("ldpca/6336_regDeg3.lad", this);
#endif
  _ldpca                 = new LdpcaEnc("ldpca/1584_regDeg3.lad", this);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::encodeKeyFrame()
{
  string srcFileName = _files->getFile("src")->getFileName();
  string keyFileName = _files->getFile("key")->getFileName();

  cout << "Running JM to encode key frames" << endl;

  stringstream cmd(stringstream::in | stringstream::out);

  cmd << "cd jm; ";
  cmd << "./lencod.exe -d encoder_intra_main.cfg ";
  cmd << "-p InputFile= \"../../bin/" << srcFileName << "\" ";
  cmd << "-p ReconFile=" << "\"../../bin/" << keyFileName << "\" ";
  cmd << "-p FramesToBeEncoded=" << ((_numFrames + _gop/2)/_gop) << " ";
  cmd << "-p QPISlice=" << _keyQp << " ";
  cmd << "-p FrameSkip=" << _gop-1 << " ";
  cmd << "-p SourceWidth=" << Y_WIDTH << " ";
  cmd << "-p SourceHeight=" << Y_HEIGHT << " ";
  cmd << " > jm.log";

  system(cmd.str().c_str());

  cout << "Done encoding key frames" << endl << endl;

  _files->getFile("key")->openFile("rb");
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::encodeWzHeader()
{
  // Encode width and height in terms of macroblock
  _bs[0]->write(_qp, 8);
  _bs[0]->write(_numFrames, 16);
  _bs[0]->write(_gopLevel, 2);

  // Read from file and verify that it is correct
//  _bs[0]->restart();

  cout << "--------------------------------------------------" << endl;
  cout << "WZ frame parameters" << endl;
  cout << "--------------------------------------------------" << endl;
//  cout << "QP:     " << _bs[0]->read(8) << endl;
//  cout << "Frames: " << _bs[0]->read(16) << endl;
//  cout << "GOP:    " << (1<<_bs[0]->read(2)) << endl;
  cout << "QP:     " << _qp << endl;
  cout << "Frames: " << _numFrames << endl;
  cout << "GOP:    " << (1<<_gopLevel) << endl;
  cout << "--------------------------------------------------" << endl << endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::encodeWzFrame()
{
  FILE* fReadPtr    = _files->getFile("src")->getFileHandle();
  FILE* fKeyReadPtr = _files->getFile("key")->getFileHandle();

  File* patternFile;
  FILE* patternFh;

  clock_t timeStart, timeEnd;
  double cpuTime;
  int bplen;

  imgpel* currFrame     = _fb->getCurrFrame();
  int*    dctFrame      = _fb->getDctFrame();
  int*    quantDctFrame = _fb->getQuantDctFrame();

  int*    residue       = new int[FSIZE];

  timeStart = clock();

  encodeWzHeader();

  // Main loop
  // ---------------------------------------------------------------------------
  for (int keyFrameNo = 0; keyFrameNo < _numFrames/_gop; keyFrameNo++) {
    // Read previous key frame from the reconstructed key frame file
    fseek(fKeyReadPtr, keyFrameNo * FSIZE, SEEK_SET);
    fread(_fb->getPrevFrame(), FSIZE, 1, fKeyReadPtr);

    // Read next key frame from the reconstructed key frame file
    fseek(fKeyReadPtr, (keyFrameNo+1) * FSIZE, SEEK_SET);
    fread(_fb->getNextFrame(), FSIZE, 1, fKeyReadPtr);

    for (int il = 0; il < _gopLevel; il++) {
      int frameStep = _gop / ((il+1)<<1);
      int idx = frameStep;

      // Start encoding the WZ frame
      while (idx < _gop) {
        int wzFrameNo = keyFrameNo*_gop + idx;

        cout << "Encoding frame " << wzFrameNo << " (Wyner-Ziv frame)" << endl;

        // Read current frame from the source file
        fseek(fReadPtr, wzFrameNo * FSIZE, SEEK_SET);
        fread(currFrame, FSIZE, 1, fReadPtr);

        // ---------------------------------------------------------------------
        // STAGE 1 - Residual coding & DCT
        // ---------------------------------------------------------------------
# if RESIDUAL_CODING
        computeResidue(residue);

        for (int c = 0; c < NCHANS; c++) {
          _trans->dctTransform(CHOFFSET(residue,c), CHOFFSET(dctFrame,c), c);
# else // if !RESIDUAL_CODING
        for (int c = 0; c < NCHANS; c++) {
          _trans->dctTransform(CHOFFSET(currFrame,c), CHOFFSET(dctFrame,c), c);
# endif // RESIDUAL_CODING

          /* TODO: check what this does*/
          updateMaxValue(CHOFFSET(dctFrame,c), c);

          // ---------------------------------------------------------------------
          // STAGE 2 - Calculate quantization step size
          // ---------------------------------------------------------------------
          /* TODO: check what this does*/
          computeQuantStep(c);

          // ---------------------------------------------------------------------
          // STAGE 3 - Quantization
          // ---------------------------------------------------------------------
          _trans->quantization(CHOFFSET(dctFrame,c),
                               CHOFFSET(quantDctFrame,c), c);
          // ---------------------------------------------------------------------
          // STAGE 4 - Mode decision
          // ---------------------------------------------------------------------
# if MODE_DECISION
          selectCodingMode(CHOFFSET(quantDctFrame,c), c);
# endif // MODE_DECISION

          // ---------------------------------------------------------------------
          // STAGE 5 - Skip mode
          // ---------------------------------------------------------------------
# if SKIP_MODE
          generateSkipMask(c);

          encodeSkipMask(c);
# endif // SKIP_MODE

          // ---------------------------------------------------------------------
          // STAGE 6 - Encode (Channel/Entropy)
          // ---------------------------------------------------------------------
          bplen = (c == 0) ? Y_BPLEN : UV_BPLEN;
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

          int bits = 0;

          // Entropy encode
          if (numBands > _numChnCodeBands[c])
            bits = _cavlc->encode(quantDctFrame, _skipMask[c]);

# if HARDWARE_FLOW
          if (bits%32 != 0) {
            int dummy = 32 - (bits%32);
            _bs[c]->write(0, dummy);
          }
# endif // HARDWARE_FLOW

          // Channel encode
          encodeFrameLdpca(quantDctFrame, c);

          // ---------------------------------------------------------------------
          // STAGE 7 - Write parity and CRC bits to the bitstream
          // ---------------------------------------------------------------------
# if RESIDUAL_CODING
          for (int i = 0; i < _rcBitPlaneNum; i++)
# else // if !RESIDUAL_CODING
          for (int i = 0; i < BitPlaneNum[_qp]; i++)
# endif // RESIDUAL_CODING
          {
            for (int j = 0; j < bplen; j++)
              _bs[c]->write(int(_parity[c][j+i*bplen]), 1);

# if !HARDWARE_FLOW
#   if HARDWARE_LDPC
            if (c == 0)
              for (int n = 0; n < 4; n++)
                _bs[c]->write(_crc[c][i*4+n], 8);
            else
              _bs[c]->write(_crc[c][i], 8);
#   else // if !HARDWARE_LDPC
              _bs[c]->write(_crc[c][i], 8);
#   endif // HARDWARE_LDPC
# endif // !HARDWARE_FLOW
          }
        }  // Finish encoding each channel separately
        idx += 2*frameStep;
      } // Finish encoding the WZ frame
    }
  }

  // Flush all bitstreams
  for (int c = 0; c < NCHANS; c++)
    _bs[c]->flush();

  timeEnd = clock();
  cpuTime = (timeEnd - timeStart) / CLOCKS_PER_SEC;

  cout << endl;
  cout << "--------------------------------------------------" << endl;
  cout << "Encode statistics" << endl;
  cout << "--------------------------------------------------" << endl;

  report();

  cout << "Total   encoding time: " << cpuTime << "(s)" << endl;
  cout << "Average encoding time: " << cpuTime/_numFrames << "(s)" << endl;
  cout << "--------------------------------------------------" << endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::computeResidue(int* residue)
{
  int     blockCount;
  int     frameHeight, frameWidth, frameSize;
  imgpel* refFrame;
  imgpel* backwardRefFrame = _fb->getPrevFrame();
  imgpel* forwardRefFrame  = _fb->getNextFrame();
  imgpel* currentFrame     = _fb->getCurrFrame();

/* Same as before, but now consider residual for all three channels.
 * Blocks are processed in channel order: Y,U,V */ 
# if !HARDWARE_OPT
  int* dirList[NCHANS];
  dirList[0] = new int[Y_FSIZE/64];
  dirList[1] = new int[UV_FSIZE/64];
  dirList[2] = new int[UV_FSIZE/64];

  memset(dirList[0], 0, FSIZE/16);
  memset(dirList[1], 0, FSIZE/16);
  memset(dirList[2], 0, FSIZE/16);
# endif

  for (int c = 0; c < NCHANS; c++) { 
    blockCount = 0;
    frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
    frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;
    frameSize = (c == 0) ? Y_FSIZE : UV_FSIZE;
    for (int j = 0; j < frameHeight; j += ResidualBlockSize)
      for (int i = 0; i < frameWidth; i += ResidualBlockSize) {
# if HARDWARE_OPT
        refFrame = backwardRefFrame;
# else // if !HARDWARE_OPT
        int bckDist; // backward distortion
        int fwdDist; // forward distortion

        bckDist = computeSad(CHOFFSET(backwardRefFrame,c) + i + j*frameWidth,
                             CHOFFSET(currentFrame,c) + i + j*frameHeight,
                             frameWidth, frameWidth, 1, 1, ResidualBlockSize);

        fwdDist = computeSad(CHOFFSET(forwardRefFrame,c) + i + j*frameWidth,
                             CHOFFSET(currentFrame,c) + i + j*frameHeight,
                             frameWidth, frameWidth, 1, 1, ResidualBlockSize);

        dirList[c][blockCount] = (bckDist <= fwdDist) ? 0 : 1;

        refFrame = (dirList[c][blockCount] == 0) ? backwardRefFrame : forwardRefFrame;
# endif // HARDWARE_OPT

        for (int y = 0; y < ResidualBlockSize; y++)
          for (int x = 0; x < ResidualBlockSize; x++) {
            int idx = (i+x) + (j+y)*frameWidth;

            CHOFFSET(residue,c)[idx] = CHOFFSET(currentFrame,c)[idx] 
                                       - CHOFFSET(refFrame,c)[idx];
          }

# if !HARDWARE_FLOW
        // Encode motion vector
        _bs[c]->write(dirList[c][blockCount], 1);
# endif
        blockCount++;
      }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int Encoder::computeSad(imgpel* blk1, imgpel* blk2, int width1, int width2, int step1, int step2, int blockSize)
{
  int sad = 0;

  for (int y = 0; y < blockSize; y++)
    for (int x = 0; x < blockSize; x++) {
      imgpel pel1 = *(blk1 + step1*x + step1*y*width1);
      imgpel pel2 = *(blk2 + step2*x + step2*y*width2);

      sad += abs(pel1 - pel2);
    }

  return sad;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::updateMaxValue(int* block, int c)
{
  int frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
  int frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;
  for (int y = 0; y < 4; y++)
    for (int x = 0; x < 4; x++)
      _maxValue[c][y][x] = 0;

  for (int y = 0; y < frameHeight; y += 4)
    for (int x = 0; x < frameWidth; x += 4)
      for (int j = 0; j < 4; j++)
        for (int i = 0; i < 4; i++)
          if (abs(_maxValue[c][j][i]) < abs(block[(x+i) + (y+j)*frameWidth]))
            _maxValue[c][j][i] = block[(x+i) + (y+j)*frameWidth];
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::computeQuantStep(int c)
{
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 4; i++) {
# if RESIDUAL_CODING

#   if !HARDWARE_FLOW
      _bs[c]->write(abs(_maxValue[c][j][i]), 11);
#   endif

      if (QuantMatrix[_qp][j][i] != 0) {
#   if HARDWARE_QUANTIZATION
        _quantStep[c][j][i] = 1 << (MaxBitPlane[j][i]+1-QuantMatrix[_qp][j][i]);
#   else // if !HARDWARE_QUANTIZATION
        int iInterval = 1 << QuantMatrix[_qp][j][i];

        _quantStep[c][j][i] = (int)ceil(double(2*abs(_maxValue[c][j][i]))/double(iInterval-1));
        _quantStep[c][j][i] = Max(_quantStep[c][j][i], MinQStepSize[_qp][j][i]);
#   endif // HARDWARE_QUANTIZATION
      }
      else
        _quantStep[c][j][i] = 1;

# else // if !RESIDUAL_CODING

      if (i != 0 || j != 0) {
        _bs[c]->write(abs(_maxValue[c][j][i]), 11);

#   if HARDWARE_QUANTIZATION
        if (QuantMatrix[_qp][j][i] != 0)
          _quantStep[c][j][i] = 1 << (MaxBitPlane[j][i]+1-QuantMatrix[_qp][j][i]);
        else
          _quantStep[c][j][i] = 0;
#   else // if !HARDWARE_QUANTIZATION
        int iInterval = 1 << (QuantMatrix[_qp][j][i]);

#     if AC_QSTEP
        if (QuantMatrix[_qp][j][i] != 0) {
          _quantStep[c][j][i] = (int)ceil(double(2*abs(_maxValue[c][j][i]))/double(iInterval-1));

          if (_quantStep[c][j][i] < 0)
            _quantStep[c][j][i] = 0;
        }
        else
          _quantStep[c][j][i] = 1;
#     else // if !AC_QSTEP
        if (QuantMatrix[_qp][j][i] != 0)
          _quantStep[c][j][i] = ceil(double(2*abs(_maxValue[c][j][i]))/double(iInterval));
        else
          _quantStep[c][j][i] = 1;
#     endif // AC_QSTEP

#   endif // HARDWARE_QUANTIZATION
      }
      else
        _quantStep[c][j][i] = 1 << (DC_BITDEPTH-QuantMatrix[_qp][j][i]);
# endif // RESIDUAL_CODING
    }
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// TODO: modify to select for each block in YUV
void Encoder::selectCodingMode(int* frame, int c)
{
  int numBands = 0;
  int mode;
  int codingMode;
  int frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
  int frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;
  double th1 = 0.15;
  double th2 = 0.05;
  double th3 = 0.01;

  for (int bandNo = 0; bandNo < 16; bandNo++) {
    int i = ScanOrder[bandNo][0];
    int j = ScanOrder[bandNo][1];
    if (QuantMatrix[_qp][j][i] > 0)
      numBands++;
  }

  memset(_average[c], 0x00, 16*sizeof(double));
  for (int j = 0; j < frameHeight; j++) {
    for (int i = 0; i < frameWidth; i++) {
      int mask, data;
      mask = (0x1 << (QuantMatrix[_qp][j%4][i%4]-1)) - 1;
      data = frame[i + j*frameWidth] & mask;
      _average[c][(i%4) + (j%4)*4] += data;
    }
  }

  for (int i = 0; i < 16; i++)
    if (c == 0)
      _average[c][i] /= Y_BPLEN;
    else
      _average[c][i] /= UV_BPLEN;

  double energy1  = 0.0;
  double energy2  = 0.0;
  double energy3  = 0.0;

  for (int i = 0; i < 16; i++) {
    int x = ScanOrder[i][0];
    int y = ScanOrder[i][1];

    if (i < 3)
      energy1 += _average[c][x + y*4];
    else if (i < 6)
      energy2 += _average[c][x + y*4];
    else
      energy3 += _average[c][x + y*4];
  }

  energy1 /= 3;
  if (numBands > 3) {
    if (numBands >= 6)
      energy2 /= 3;
    else
      energy2 /= (numBands-3);
  }
  else
    energy2 = 0;

  if (numBands > 6)
    energy3 /= (numBands-6);
  else
    energy3 = 0;

  if (energy1 > (th1/(double)(Scale[0][_qp]))) {
    if (energy2 > (th2/(double)(Scale[1][_qp]))) {
      if (energy3 > (th3/(double)(Scale[2][_qp])))
        mode = 0; // channel coding (channel coding for all bands)
      else
        mode = 2; // hybrid mode 2 (channel coding for lower 6 bands
                        //                entropy coding for other bands)
    }
    else
      mode = 1;   // hybrid mode 1 (channel coding for lower 3 bands
                        //                entropy coding for other bands)
  }
  else
    mode = 3;     // entropy coding (entropy coding for all bands)

  _modeCounter[c][mode]++;

  codingMode = mode;

  _bs[c]->write(codingMode, 2);

  if (codingMode == 0) _numChnCodeBands[c] = 16; else
  if (codingMode == 1) _numChnCodeBands[c] =  3; else
  if (codingMode == 2) _numChnCodeBands[c] =  6; else
                       _numChnCodeBands[c] =  0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::generateSkipMask(int c)
{
# if HARDWARE_OPT
  int threshold = 2;
# else // if !HARDWARE_OPT
  int threshold = 5;
# endif // HARDWARE_OPT

# if RESIDUAL_CODING
  int* frame    = _fb->getQuantDctFrame();
# else // if !RESIDUAL_CODING
  int* frame    = new int[FSIZE];
  int* frameDct = new int[FSIZE];

  for (int i = 0; i < FSIZE; i++)
    frame[i] = _fb->getPrevFrame()[i] - _fb->getCurrFrame()[i];

  _trans->dctTransform(CHOFFSET(frame,c), CHOFFSET(frameDct,c), c);
  _trans->quantization(CHOFFSET(frameDct,c), CHOFFSET(frame,c), c);
# endif // !RESIDUAL_CODING

  int frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
  int frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;
  int bplen = (c == 0) ? Y_BPLEN : UV_BPLEN;
  memset(_skipMask[c], 0x0, sizeof(int)*bplen);
  for (int j = 0; j < frameHeight; j += SkipBlockSize) {
    for (int i = 0; i < frameWidth; i += SkipBlockSize) {
      int distortion = 0;
      int blockIndex = i/SkipBlockSize + (j/SkipBlockSize)*(frameWidth/SkipBlockSize);

      for (int y = 0; y < SkipBlockSize; y++) {
        for (int x = 0; x < SkipBlockSize; x++) {
          int mask, data;

          mask = (0x1 << (QuantMatrix[_qp][y][x]-1)) - 1;
          data = CHOFFSET(frame,c)[(i+x) + (j+y)*frameWidth] & mask;

# if HARDWARE_OPT
          distortion += data;
# else // if !HARDWARE_OPT
          distortion += data * data;
# endif // HARDWARE_OPT
        }
      }
      _skipMask[c][blockIndex] = (distortion < threshold) ? 1 : 0;
    }
  }

# if !RESIDUAL_CODING
  delete [] frame;
  delete [] frameDct;
# endif // !RESIDUAL_CODING
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::encodeSkipMask(int c)
{
  // n0 = number of non-skipped blocks
  int code, length;
  int index = 0, n0 = 0, run = 0, bitCount = 0;
  int frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
  int frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;
  int bplen = (c == 0) ? Y_BPLEN : UV_BPLEN;
# if HARDWARE_FLOW
  // Pad zero
  int dummy = 20;
  _bs[c]->write(0, dummy);
# endif // HARDWARE_FLOW

  for (int i = 0; i < bplen; i++)
    if (_skipMask[c][i] == 0)
      n0++;

  int diff = abs(n0 - bplen/2);

# if HARDWARE_OPT
  int type = _prevType[c];
  int nextType = diff / (bplen/6);
  _prevType[c] = nextType;
# else // if !HARDWARE_OPT
  int type = diff / (bplen/6);
# endif // HARDWARE_OPT

  if (type > 2) type = 2;

  _bs[c]->write(type, 2);

  bitCount = 2;

// Directly output first bit to bitstream
  int sign = _skipMask[c][0];

  _bs[c]->write(sign, 1);

  bitCount++;
  run++;
  index++;

# if TESTPATTERN
  File* patternFile;
  FILE* patternFh;

  patternFile = _files->addFile("pattern_rlc", "pattern_rlc.dat");
  patternFile->openFile("w");
  patternFh = patternFile->getFileHandle();

  for (int idx = 1; idx <= 8; idx++) {
    int data = (sign >> (32-idx*4)) & 0xf;
    fprintf(patternFh, "%x", data);
  }
  fprintf(patternFh, "\n");
# endif // TESTPATTERN

  // Huffman code for other bits
  while (index < bplen) {
    if (_skipMask[c][index] == sign) {
      run++;

      if (run == 16) { // reach maximum run length
        bitCount += getHuffmanCode(_qp, type, run-1, code, length);

        _bs[c]->write(code, length);

        run = 1;

# if TESTPATTERN
        for (int idx = 1; idx <= 8; idx++) {
          int data = (code >> (32-idx*4)) & 0xf;
          fprintf(patternFh, "%x", data);
        }
        fprintf(patternFh, "\n");
# endif // TESTPATTERN
      }
    }
    else {
      bitCount += getHuffmanCode(_qp, type, run-1, code, length);

      _bs[c]->write(code, length);

      sign = _skipMask[c][index];
      run = 1;

# if TESTPATTERN
      for (int idx = 1; idx <= 8; idx++) {
        int data = (code >> (32-idx*4)) & 0xf;
        fprintf(patternFh, "%x", data);
      }
      fprintf(patternFh, "\n");
# endif // TESTPATTERN
    }

    index++;
  }

  if (run != 0) {
    bitCount += getHuffmanCode(_qp, type, run-1, code, length);

    _bs[c]->write(code, length);

# if TESTPATTERN
    for (int idx = 1; idx <= 8; idx++) {
      int data = (code >> (32-idx*4)) & 0xf;
      fprintf(patternFh, "%x", data);
    }
    fprintf(patternFh, "\n");
# endif // TESTPATTERN
  }

# if TESTPATTERN
  patternFile->closeFile();
# endif // TESTPATTERN

# if HARDWARE_FLOW
  if (bitCount%32 != 0) { // pad zero
    int dummy = 32 - (bitCount%32);
    _bs[c]->write(0, dummy);
  }
# endif // HARDWARE_FLOW
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int Encoder::getHuffmanCode(int qp, int type, int symbol, int& code, int& length)
{
  int table = qp / 2;

  code   = HuffmanCodeValue [table][type][symbol];
  length = HuffmanCodeLength[table][type][symbol];

  return length;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::encodeFrameLdpca(int* frame, int c)
{
  int   bitPosition;
  bool* accumulatedSyndrome = _parity[c];
  int   frameSize = (c == 0) ? Y_FSIZE : UV_FSIZE;
  int   bplen = (c == 0) ? Y_BPLEN : UV_BPLEN;
  int*  ldpcaSource = new int[bplen + 8];
  for (int band = 0; band < _numChnCodeBands[c]; band++) {
    int i = ScanOrder[band][0];
    int j = ScanOrder[band][1];

# if RESIDUAL_CODING
    for (bitPosition = _rcQuantMatrix[j][i]-1; bitPosition >= 0; bitPosition--)
# else // if !RESIDUAL_CODING
    for (bitPosition = QuantMatrix[_qp][j][i]-1; bitPosition >= 0; bitPosition--)
# endif // RESIDUAL_CODING
    {
      setupLdpcaSource(frame, ldpcaSource, i, j, bitPosition, c);

      if (c == 0) {
# if HARDWARE_LDPC
        for (int n = 0; n < 4; n++) {
          _ldpca->encode(ldpcaSource + n*1584, accumulatedSyndrome);

          computeCRC(ldpcaSource + n*1584, 1584, _crc[c]+n);

          accumulatedSyndrome += bplen/4;
        }

        _crc[c] += 4;
        cout << ".";
# else // if !HARDWARE_LDPC
      _ldpca_cif->encode(ldpcaSource, accumulatedSyndrome);

      cout << ".";

      computeCRC(ldpcaSource, bplen, _crc[c]);

      accumulatedSyndrome += frameSize/16;

      _crc[c]++;
# endif // HARDWARE_LDPC
      } else {
        _ldpca->encode(ldpcaSource, accumulatedSyndrome);

        cout << ".";

        computeCRC(ldpcaSource, bplen, _crc[c]);

        accumulatedSyndrome += frameSize/16;

        _crc[c]++;
      }
    }
  }
  delete [] ldpcaSource;
  cout << endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::setupLdpcaSource(int* frame, int* source,
                               int offsetX, int offsetY,
                               int bitPosition, int c)
{
  int frameWidth = (c == 0) ? Y_WIDTH : UV_WIDTH;
  int frameHeight = (c == 0) ? Y_HEIGHT : UV_HEIGHT;
  for (int y = 0; y < frameHeight; y = y+4)
    for (int x = 0; x < frameWidth; x = x+4) {
      int blockIdx = (x/4) + (y/4)*(frameWidth/4);
      int frameIdx = (x+offsetX) + (y+offsetY)*frameWidth;

# if SKIP_MODE

      if (_skipMask[c][blockIdx] == 1)
        source[blockIdx] = 0;
      else
        source[blockIdx] = (CHOFFSET(frame,c)[frameIdx] >> bitPosition) & 0x1;

# else // if !SKIP_MODE
      source[blockIdx] = (CHOFFSET(frame,c)[frameIdx] >> bitPosition) & 0x1;
# endif // SKIP_MODE

    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::computeCRC(int* data, const int length, unsigned char* crc)
{
  // CRC8 110011011
  const int code[9] = {1, 1, 0, 0, 1, 1, 0, 1, 1};

  int* buffer = new int[length + 8];

  memcpy(buffer, data, length*sizeof(int));

  for (int i = length; i < length+8; i++)
    buffer[i] = 0;

  for (int i = 0; i < length; i++)
    if (buffer[i] == 1)
      for (int j = 0; j < 9; j++)
        buffer[i+j] = code[j] ^ buffer[i+j];

  *crc = 0;

  for (int i = 0; i < 8; i++)
    *crc |= buffer[length+i] << (7-i);

  delete [] buffer;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Encoder::report()
{
# if MODE_DECISION
  char colour_chan[NCHANS] = {'Y', 'U', 'V'};
  cout << "Mode usage: " <<endl;
  for (int c = 0; c < NCHANS; c++) {
  cout << colour_chan[c] << ": ";
    for (int i = 0; i < 4; i++) {
      float usage = (float)_modeCounter[c][i]/75.0 * 100.0;
      cout << usage << " ";
    }
    cout << endl;
  }

  cout << endl;
# endif // MODE_DECISION
}

