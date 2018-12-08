
#include <iostream>

#include "encoder.h"

using namespace std;

int main(int argc, char** argv)
{
  if (argc != 8) {
    cerr << endl;
    cerr << "Usage: ./enDVC ";
    cerr << "[WZ QP] [key QP] ";
    cerr << "[frame number] [GOP level] ";
    cerr << "[input file] [output file] [recon file]" << endl;
    cerr << endl;
  }
  else {
    cout << endl;
    cout << "DVC2.0 - open source distributed video coding tool" << endl;
    cout << endl;

    Encoder* encoder = new Encoder(argv);

    encoder->encodeKeyFrame();
    encoder->encodeWzFrame();

    cout << endl;
    cout << "Bye" << endl;
    cout << endl;
  }

  //system("pause");
  return 0;
}

