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
enum Method {MVSAME=1, SIMPLE=2, MVSEARCH=3, MVMULTI=4};

// Base Class
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class Colourizer
{
public:
  Colourizer(map<string, string> configMap);
  virtual ~Colourizer() {};

  int getWidth() {return _width;}
  int getHeight() {return _height;}
  virtual void addColour(imgpel* prevKeyFrame, imgpel* currFrame);
  virtual void colourize();

protected:
  virtual void MC(imgpel* prevKeyFrame, imgpel* currFrame);

protected:
  FileManager* _files;
  mvinfo* _mvs;
  int _grayFrameSize;
  int _yuvFrameSize;
  int _method;
  int _gop;
  int _gopLevel;
  int _width;
  int _height;
  int _nframes;
  int _iRange;
  int _nMV;

};

// Subclasses
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------


/**
 * Description: colourizes using block-based motion estimation
 * on the previous frame.
 *
 * Base Class: Colourizer
 */
class MvSearchColour:public Colourizer
{
public:
  MvSearchColour(map<string, string> configMap);

  ~MvSearchColour() {delete [] _mvs;}

  void addColour(imgpel* prevKeyFrame, imgpel* currFrame);

private:
  int _param;
  int _thresh;
};


/**
 * Description: colourizes using the same bilinear motion vectors that were
 * used to generate the Luma side-information. Reads MV's from file.
 *
 * Base Class: Colourizer
 */
class MvSameColour:public Colourizer
{
public:
  MvSameColour(map<string, string> configMap);

  ~MvSameColour() {_mvFile.close(); delete [] _mvs;}

  void addColour(imgpel* prevKeyFrame, imgpel* nextKeyFrame, imgpel* currFrame);

private:
  void readMVFile();

  void bidirectional_mc(imgpel* prevKeyFrame,
                                      imgpel* nextKeyFrame,
                                      imgpel* currFrame);


private:
  ifstream _mvFile;
};


/**
 * Description: colourizes by performing motion estimation using a number of
 * other candidate frames.
 *
 * Base Class: Colourizer
 */
class MvMultiFrame:public Colourizer
{
public:
  MvMultiFrame(map<string, string> configMap);

  ~MvMultiFrame() {delete [] _mvs;}

  void colourize();

  void addColour(imgpel** keyFrameList, imgpel* currFrame);

private:
  void MC(imgpel** keyFrameList, imgpel* currFrame);

private:
  int _param;
  int _nCandidates;
};


// Helper functions
map<string, string>& readConfig(string filename);

void bilinear(imgpel *source, imgpel *buffer, int width, int height);

#endif //COULOURIZER_INC_COLOURIZER
