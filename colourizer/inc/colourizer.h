#ifndef COULOURIZER_INC_COLOURIZER
#define COULOURIZER_INC_COLOURIZER

#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>

#include "codec.h"
#include "fileManager.h"

using namespace std;

class FileManager;
class SideInformation;
enum Method {SIMPLE=1, MVSEARCH=2, MVSAME=3};

// Base Class
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class Colourizer
{
public:
  Colourizer(map<string, string> configMap);
  virtual ~Colourizer() {};

  void colourize();
  virtual void addColour(imgpel* prevKeyFrame, imgpel* currFrame);
  int getWidth() {return _width;}
  int getHeight() {return _height;}
  void MC(imgpel* imgPrevUV, imgpel* imgCurrUV);

protected:
  FileManager* _files;
  mvinfo* _mvs;
  int _grayFrameSize;
  int _yuvFrameSize;
  int _method;
  int _gop;
  int _width;
  int _height;
  int _nframes;
  int _iRange;
  int _nMV;

};

// Subclasses
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class MvSearchColour:public Colourizer
{
public:
  MvSearchColour(map<string, string> configMap);

  ~MvSearchColour() {delete [] _mvs;}

  void forwardME(imgpel* prevKeyFrame, imgpel* currFrame);

  void addColour(imgpel* prevKeyFrame, imgpel* currFrame);

private:
  int _param;

};

class MvSameColour:public Colourizer
{
public:
  MvSameColour(map<string, string> configMap);

  ~MvSameColour() {_mvFile.close(); delete [] _mvs;}

  void addColour(imgpel* prevKeyFrame, imgpel* currFrame);

private:
  void readMVFile();
private:
  ifstream _mvFile;
};

map<string, string>& readConfig(string filename);

#endif //COULOURIZER_INC_COLOURIZER
