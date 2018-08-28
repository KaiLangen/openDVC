#include <cstdlib>
#include <climits>

#include "ME.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int calcSAD(imgpel* blk1, imgpel* blk2, int px, int py, 
            int rx,int ry, const int blocksize, int width, int height)
{
  int sad=0;
  int cx,cy;
  for(int y=0;y<blocksize;y++)
    for(int x=0;x<blocksize;x++)
    {
      cx=px+x;
      cy=py+y;
      if(cx<=0)cx=0;
      if(cx>width-1)cx=width-1;
      if(cy<=0)cy=0;
      if(cy>height-1)cy=height-1;

      imgpel pel1=*(blk1+cx+cy*width);
      imgpel pel2=*(blk2+(rx+x)+(ry+y)*width);
      sad+=abs(pel1-pel2);
    }
  return sad;
}

int calcSAD(imgpel* blk1, imgpel* blk2, int width1, int width2, int s1,int s2,int blocksize)
{
  int sad=0;
  for(int y=0;y<blocksize;y++)
    for(int x=0;x<blocksize;x++)
    {
      imgpel pel1=*(blk1+s1*x+s1*y*width1);
      imgpel pel2=*(blk2+s2*x+s2*y*width2);
      sad+=abs(pel1-pel2);
    }
  return sad;
}

int calcSAD(imgpel* blk1, imgpel* blk2, int width, int blocksize)
{
  int sad=0;
  for(int y=0;y<blocksize;y++)
    for(int x=0;x<blocksize;x++)
    {
      imgpel pel1=*(blk1+x+y*width);
      imgpel pel2=*(blk2+x+y*width);
      sad+=abs(pel1-pel2);
    }
  return sad;
}

static int find_min(unsigned int costs[9])
{
    unsigned int minimum = costs[0];
    int location = 0;
    for(int c = 1; c < 9; ++c)
    {
        if(costs[c] < minimum)
        {
            minimum = costs[c];
            location = c;
        }
    }
    return location;
}

void TSS(imgpel* trg, imgpel* ref, mvinfo& mv,
        int step, int center, int width, int height, int blocksize)
{
    // search start location
    unsigned int costs[9] = {UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,
                   UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX};
    int locations[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
    int loc, og, cx, cy, x, y;
    og = center;
    // calculate the first center
    // avoid recalculating the center within the loop
    costs[4] = calcSAD(&trg[og], &ref[center], width, blocksize);
    locations[4] = center;

    while(step >= 1)
    {
        // center coordinates in the image = (cy, cx)
        cy = center / width;
        cx = center % width;
        // coordinates in the cost matrix = (i,j)
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 3; ++j)
            {
                // the 9 pts formed by stepping away from the center = (y, x)
                //
                // (cy-step, cx-step), (cy-step,      cx), (cy-step, cx+step)
                // (cy     , cx-step), (cy     ,      cx), (cy     , cx+step)
                // (cy+step, cx-step), (cy+step,      cx), (cy+step, cx+step)
                y = cy + (i-1) * step;
                x = cx + (j-1) * step;

                // check if the pt coordinates fall outside of the image
                if(x < 0 || x >= width - blocksize ||
                   y < 0 || y >= height - blocksize ||
                   (i == 1 && j == 1))
                {
//                    costs[i*3 + j] = USHRT_MAX;
                    continue;
                }
                costs[i*3 + j] = calcSAD(&trg[og],
                                         &ref[y*width + x],
                                         width,
                                         blocksize);
                locations[i*3 + j] = y*width + x;
            }
        }
        // re-center the search window on the local minimum
        loc = find_min(costs);
        center = locations[loc];
        step /= 2;

        // set the center and location
        costs[4] = costs[loc];
        locations[4] = center;
    }

    x = og % width;
    y = og / width;
    cx = center % width;
    cy = center / width;
    // MV is new location - original location
    mv.iMvx = cx - x;
    mv.iMvy = cy - y;

}

//int SideInformation::calcSAD(imgpel* blk1, imgpel* blk2,const int blocksize,const int iPadSize){
//  int iWidth;
//  iWidth  = _codec->getFrameWidth();
//  int sad=0;
//  for(int y=0;y<blocksize;y++)
//    for(int x=0;x<blocksize;x++)
//    {
//      imgpel pel1=*(blk1+x+y*(iWidth+2*iPadSize));
//      imgpel pel2=*(blk2+x+y*(iWidth+2*iPadSize));
//      sad+=abs(pel1-pel2);
//    }
//  return sad;
//}
