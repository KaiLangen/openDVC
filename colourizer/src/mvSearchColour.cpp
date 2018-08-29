#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "ME.h"
#include "colourizer.h"

using namespace std;

MvSearchColour::MvSearchColour(map<string, string> configMap)
  : Colourizer(configMap)
{
  _iRange = atoi(configMap["BlockSize"].c_str());
  _nMV    = _width * _height / (_iRange*_iRange);
  _mvs    = new mvinfo[_nMV];
  _param  = atoi(configMap["Param"].c_str());
  _files->addFile("mv", configMap["MVFile"])->openFile("w");
}


/**
 * Description:
 * Recolour the WZ frames using the motion-compensated chroma channels
 * from the Key frames. Motion vectors are generated using the luma channel.
 * ----------------------------------------------------------------------------
 * Param: None
 * Return: None
 */
void MvSearchColour::addColour(imgpel* prevKeyFrame, imgpel* currFrame)
{
  int index, pos, param, n;
  char motionVectorBuffer[100];
  FILE* mvFilePtr = _files->getFile("mv")->getFileHandle();
  // take the nearest power of two to  param to find
  // the step size in pixels
  double L = floor(log2(_param + 1.0));
  double stepMax = pow(2.0, (L-1.0));
  int step = (int)stepMax;
  for (int y = 0; y < _height; y+=_iRange) {
    for (int x = 0; x < _width; x+=_iRange) {
      pos = x+y*_width;
      index = x/_iRange + y*_width/(_iRange*_iRange);
      _mvs[index].iCx = x;
      _mvs[index].iCy = y;
      TSS(currFrame, prevKeyFrame, _mvs[index],
          step, pos, _width, _height, _iRange);
//      cout << "Motion vector at index: " << index << endl;
//      cout << _mvs[index].iCx << " "; 
//      cout << _mvs[index].iCy << " "; 
//      cout << _mvs[index].iMvx << " ";
//      cout << _mvs[index].iMvy << " "; 
//      cout << endl; 
//      n = sprintf(motionVectorBuffer, "%d, %d, %d, %d\n",
//                  _mvs[index].iCx,
//                  _mvs[index].iCy,
//                  _mvs[index].iMvx,
//                  _mvs[index].iMvy);
//      fwrite(motionVectorBuffer, n, 1, mvFilePtr);
    }
  }
  MC(prevKeyFrame, currFrame);
}
