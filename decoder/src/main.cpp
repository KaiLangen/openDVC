
#include <iostream>

#include "decoder.h"

using namespace std;

int main(int argc, char** argv)
{
  if (argc != 6 and argc != 7) {
    cerr << "Usage: ./deDVC ";
    cerr << "[wz varBitstream file] ";
    cerr << "[key frame file] ";
    cerr << "[original video file] ";
    cerr << "[channel] ";
    cerr << "[SI Method] ";
    cerr << "[helper file]" << endl;
    cerr << endl;
    return 1;
  }
  else {
    cout << endl;
    cout << "DVC2.0 - open source distributed video coding tool" << endl;
    cout << endl;

    Decoder* decoder = new Decoder(argv);

    decoder->decodeWZframe();

    cout << endl;
    cout << "Bye" << endl;
    cout << endl;
  }

  //system("pause");
  return 0;
}

