#include <iostream>

#include "colourizer.h"

using namespace std;

int main(int argc, char** argv)
{
  if (argc != 2) {
    cerr << "Usage: ./colour ";
    cerr << "[config file]";
    cerr << endl;
    return 1;
  } else {
    cout << endl;
    cout << "DVC2.0 - open source distributed video coding tool" << endl;
    cout << endl;

    map<string, string> configMap = readConfig(argv[1]);
    Colourizer* col; 
    int method = atoi(configMap["Method"].c_str());
    switch(method){
      case SIMPLE:
        col = new Colourizer(configMap);
        break;
      case MVSEARCH:
        col = new MvSearchColour(configMap);
        break;
      case MVSAME:
        col = new MvSameColour(configMap);
        break;
      default:
        cerr << "Invalid recolouring method: " << method << endl;
        return 1;
    }
    col->colourize();

    cout << endl;
    cout << "Bye" << endl;
    cout << endl;
    delete col;
  }

}
