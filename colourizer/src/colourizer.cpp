#include <sstream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "colourizer.h"
#include "fileManager.h"

using namespace std;

Colourizer::Colourizer(char **argv)
{
  _files = FileManager::getManager();
  map<string, string> configMap = readConfig(argv[1]);
  initialize(configMap);
}

map<string, string>&
Colourizer::readConfig(string filename)
{
  string line;
  ifstream cfile(filename);
  static map<string, string> configMap;
  if (cfile)
  {
    while( getline(cfile, line) )
    {
      string key, value;
      istringstream is_line(line);

      // Check if the first non-whitespace is a #
      if( getline(is_line >> ws, key, '=')  && !key.empty() && key[0] != '#')
      {
        if( getline(is_line >> ws, value) ) 
          configMap[key] = value;
      }
    }

    for (auto const& x : configMap)
    {
      cout << x.first
           << ':'
           << x.second
           << endl;
    }
    cfile.close();
    return configMap;
  }
  else
  {
    cerr << "Invalid config file!" << endl;
    return configMap;
  }
}

void
Colourizer::initialize(map<string, string>& configMap)
{
  // find sequence type in config map
  auto it = configMap.find("SequenceType");
  if (it != configMap.end() && strcmp(it->second.c_str(), "CIF")) {
   _width = 352;
   _height = 288;
  } else if (it != configMap.end() && strcmp(it->second.c_str(), "QCIF")) {
   _width = 352;
   _height = 288;
  } else {
    cerr << "Invalid sequence type!" << endl;
  }

  _gop = atoi(configMap["Gop"].c_str());
  _method = atoi(configMap["Method"].c_str());
  _nframes = atoi(configMap["FramesToBeEncoded"].c_str());
  _files->addFile("wz", configMap["WZFile"]);
  _files->addFile("key", configMap["KeyFile"]);

  _grayFrameSize = _width * _height;
  _colourFrameSize = 3*(_grayFrameSize)>>1;
}

void
Colourizer::colourize()
{
  switch(_method)
  {
    case SIMPLE:
      cout << "Simple" << endl;
      break;
    case MVSEARCH:
      cout << "MVSearch" << endl;
      break;
    case ML:
      cout << "ML" << endl;
      break;
    default:
      cout << "Invalid colourizing method" << endl;
      break;
  }
}
