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
    throw std::invalid_argument("Invalid MV file");
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

void MvSameColour::addColour(imgpel* prevKeyFrame, imgpel* currFrame)
{
  readMVFile();
  MC(prevKeyFrame + _grayFrameSize,
     currFrame + _grayFrameSize);
}
