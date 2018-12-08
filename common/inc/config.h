
#ifndef COMMON_INC_CONFIG_H
#define COMMON_INC_CONFIG_H

#define PI                      3.14159265
#define BLOCK_SIZE              4
#define MB_BLOCK_SIZE           16
#define DC_BITDEPTH             10

// Macros for encoder/decoder
/////////////////////////////////
#define WHOLEFLOW               1
#define AC_QSTEP                1
#define RESIDUAL_CODING         1
#define SKIP_MODE               1
#define MODE_DECISION           1
#define INTEGER_DCT             1
#define HARDWARE_FLOW           1
#define HARDWARE_QUANTIZATION   1
#define HARDWARE_LDPC           1
#define HARDWARE_CMS            1
#define HARDWARE_OPT            1
#define BIDIRECT_REFINEMENT     0
#define SI_REFINEMENT           1

// frame size/offset macros
#define Y_WIDTH                 352
#define Y_HEIGHT                288
#define Y_FSIZE                 101376
#define Y_BPLEN                 6336

#define UV_WIDTH                176
#define UV_HEIGHT               144
#define UV_FSIZE                25344
#define UV_BPLEN                1584

#define FSIZE                   152064
#define NBANDS                  16
#define NCHANS                  3

const static int channel_start[NCHANS] = {0, 101376, 126720};
#define CHOFFSET(F,X) (F + channel_start[X]) 
// Macros for encoder only
/////////////////////////////////
#ifdef ENCODER
# define TESTPATTERN            1
# define DEBUG                  0
#endif

// Macros for decoder only
/////////////////////////////////
#ifdef DECODER
# define INVERSE_MATRIX         1
#endif

#include "types.h"

#endif // COMMON_INC_CONFIG_H

