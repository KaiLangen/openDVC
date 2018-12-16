
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>


using namespace std;

int main(int argc, char** argv)
{
  if (argc != 6) {
    cerr << endl;
    cerr << "Usage: ./mux ";
    cerr << "[Y file] [U file] [V file] [output file] [# of frames]" << endl;
    cerr << endl;
  }
  else {
    cout << endl;
    cout << "DVC2.0 - open source distributed video coding tool" << endl;
    cout << endl;

    FILE* yfile = fopen(argv[1], "rb");
    FILE* ufile = fopen(argv[2], "rb");
    FILE* vfile = fopen(argv[3], "rb");
    FILE* outfile = fopen(argv[4], "wb");
    int nframes = atoi(argv[5]);
    const int YSIZE = 101376;
    const int UVSIZE = 25344;
    char yframe[YSIZE];
    char uframe[UVSIZE];
    char vframe[UVSIZE];

    for (int i = 0; i < nframes; i++) {
      fread(yframe, YSIZE, 1, yfile);
      fwrite(yframe, YSIZE, 1, outfile);
      fread(uframe, UVSIZE, 1, ufile);
      fwrite(uframe, UVSIZE, 1, outfile);
      fread(vframe, UVSIZE, 1, vfile);
      fwrite(vframe, UVSIZE, 1, outfile);
    }
    fclose(yfile);
    fclose(ufile);
    fclose(vfile);
    fclose(outfile);
  }
  return 0;
}

