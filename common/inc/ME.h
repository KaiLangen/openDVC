#ifndef COMMON_INC_ME_H
#define COMMON_INC_ME_H

#include "types.h"

int calcSAD(imgpel* blk1, imgpel* blk2, int px, int py, 
            int rx,int ry, const int blocksize, int width, int height);

int calcSAD(imgpel* blk1, imgpel* blk2, int width1, int width2, int s1,int s2,int blocksize);

int calcSAD(imgpel* blk1, imgpel* blk2, int width, int blocksize);

void TSS(imgpel* trg, imgpel* ref, mvinfo& mv,
         int step, int center, int width, int height, int blocksize);


#endif // COMMON_INC_ME_H
