
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

    Colourizer* col = new Colourizer(argv);

    col->colourize();

    cout << endl;
    cout << "Bye" << endl;
    cout << endl;
  }

}
