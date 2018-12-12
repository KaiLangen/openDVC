#include <exception>

#include "colourizer.h"

using namespace std;

MvSameColour::MvSameColour(map<string, string> configMap): Colourizer(configMap)
{
  _mvFile = ifstream(configMap["MVFile"]);
  cout << configMap["MVFile"] << endl;
  _iRange=8;
  _nMV = _width * _height / (_iRange*_iRange);
  _mvs = new mvinfo[_nMV];
}

void MvSameColour::readMVFile()
{
  string line, s;
  int quad[4] = {0,0,0,0};
  // zero out the current mv list
  memset(_mvs, 0, _nMV);

  //fill mv list from file
  if (!_mvFile.good()) {
    cerr << "No such mv file" << endl;
    throw invalid_argument("Invalid MV file");
  } else {
    for (int i = 0; i < _nMV; i++) {
      getline(_mvFile, line);
      istringstream is_line(line);
      for(int j = 0; j < 4; j++) {
         if (getline(is_line, s, ','))
           quad[j] = atoi(s.c_str());
         else {
           getline(is_line, s);
           quad[j] = atoi(s.c_str());
         }
      }
      memcpy(&_mvs[i], &quad, sizeof(int)*4);
    }
  }
}


void bilinear(imgpel *source, imgpel *buffer, int width, int height)
{

  int a,b,c,d;
  int x;
  int y;
  for(int j = 0; j < height; j++)
    for(int i = 0; i < width; i++)
    {
      a=source[i + j * width];

      if((i+1)<width) b = source[(i+1)+ j*width];
      else b = a;

      if((j+1)<height) c = source[i + (j+1)*width];
      else c = a;

      if((i+1)<width && (j+1)<height) d = source[(i+1)+(j+1)*width];
      else d = a;

      buffer[2*i + (2*j)*2*width] = a;
      buffer[(2*i+1) + (2*j)*2*width] = (a+b)/2;
      buffer[2*i + (2*j+1)*2*width] = (a+c)/2;
      buffer[(2*i+1) + (2*j+1)*2*width] = (a+b+c+d)/4;
    }
}

void MvSameColour::bidirectional_mc(imgpel* prevKeyFrame, imgpel* nextKeyFrame, imgpel* currFrame)
{
  int px[2],py[2], pos, voffset;
  imgpel pel[4];
  imgpel *imgPrevBuffU, *imgPrevBuffV,*imgNextBuffU, *imgNextBuffV;

  voffset = _width * _height / 4;
  imgpel* imgPrevU = prevKeyFrame + _grayFrameSize;
  imgpel* imgPrevV = imgPrevU + voffset;
  imgpel* imgNextU = nextKeyFrame + _grayFrameSize;
  imgpel* imgNextV = imgNextU + voffset;
  imgpel* imgCurrUV = currFrame + _grayFrameSize;

  imgPrevBuffU = new imgpel[_grayFrameSize];
  imgPrevBuffV = new imgpel[_grayFrameSize];
  imgNextBuffU = new imgpel[_grayFrameSize];
  imgNextBuffV = new imgpel[_grayFrameSize];
  bilinear(imgPrevU, imgPrevBuffU, _width/2, _height/2);
  bilinear(imgPrevV, imgPrevBuffV, _width/2, _height/2);
  bilinear(imgNextU, imgNextBuffU, _width/2, _height/2);
  bilinear(imgNextV, imgNextBuffV, _width/2, _height/2);

  for(int index=0; index < _nMV; index++)
  {
    for(int j = 0; j < _iRange; j+=2)
      for(int i = 0; i < _iRange; i+=2)
      {
        px[0] = _mvs[index].iCx + i + _mvs[index].iMvx/2;
        py[0] = _mvs[index].iCy + j + _mvs[index].iMvy/2;
        px[1] = _mvs[index].iCx + i - _mvs[index].iMvx/2;
        py[1] = _mvs[index].iCy + j - _mvs[index].iMvy/2;

        pel[0] = imgPrevBuffU[px[0] + py[0] * _width];
        pel[1] = imgNextBuffU[px[1] + py[1] * _width];
        pel[2] = imgPrevBuffV[px[0] + py[0] * _width];
        pel[3] = imgNextBuffV[px[1] + py[1] * _width];

        pos = (_mvs[index].iCx+i)/2 + (_mvs[index].iCy+j)/2 * _width/2;
        if (pos >= _grayFrameSize>>2) {
					ostringstream ss;
          ss << "Error! Idx out of bounds: " << pos << endl;
          throw out_of_range(ss.str());
        }
        // U channel pixel
        imgCurrUV[pos] = (pel[0] + pel[1]+1)/2;
        // V channel pixel
        imgCurrUV[pos+voffset] = (pel[2] + pel[3]+1)/2;
      }
  }
  delete [] imgPrevBuffU;
  delete [] imgPrevBuffV;
  delete [] imgNextBuffU;
  delete [] imgNextBuffV;
}

void MvSameColour::addColour(imgpel* prevKeyFrame, imgpel* nextKeyFrame, imgpel* currFrame)
{
  readMVFile();
  bidirectional_mc(prevKeyFrame, nextKeyFrame, currFrame);
}
