#ifndef COULOURIZER_INC_COLOURIZER
#define COULOURIZER_INC_COLOURIZER

#include <map>
#include <iostream>

#include "config.h"
#include "codec.h"

using namespace std;


class FileManager;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class Colourizer
{
public:
  enum Method {SIMPLE, MVSEARCH, ML};
  Colourizer(char** argv);
  ~Colourizer() { /* TODO Remember to free memory space */ };

  void colourize();

private:
  map<string, string>& readConfig(string filename);
  void initialize(map<string, string>& configMap);

private:
  FileManager* _files;
  int _grayFrameSize;
  int _colourFrameSize;
  int _method;
  int _gop;
  int _width;
  int _height;
  int _nframes;

};
#endif //COULOURIZER_INC_COLOURIZER
