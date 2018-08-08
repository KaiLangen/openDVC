
/*!
 * @file      bitstream.h
 * @brief     Output bitstream class
 * @details   Data are hold in the internal buffer and written to the output
              file when the buffer is full or flushed
 * @author    Chieh-Chuan Chiu (r99943009@ntu.edu.tw)
 * @author    Chester Liu (u931803@gmail.com)
 */

#ifndef ENCODER_INC_BITSTREAM_H
#define ENCODER_INC_BITSTREAM_H

#include <cstdio>

#include "config.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class Bitstream
{
public:
  //! Constructor
  /*! @param[in] bufferCapacity Internal buffer size
      @param[in] fileHandle Output file handle */
  Bitstream(int bufferCapacity, FILE* fileHandle);

  ~Bitstream() { delete [] _streamBuffer; };

  //! Write data to the internal buffer
  /*! @param[in] value   Data value
      @param[in] numBits Data bit length */
  void write(int value, int numBits);

  int read(int numBits);

  bool read1Bit();

  //! Flush the internal buffer to the output file
  void flush();

  int getBitCount() { return _bitCount; };

private:
  //! Write 1-bit to the internal buffer
  /*! @param[in] bitValue Bit value */
  void write1Bit(int bitValue);

  FILE*   _fileHandle;      //!< Output file handle

  int     _numEmptyBits;    //!< Byte buffer empty bit number
  imgpel  _byteBuffer;      //!< Byte buffer

  int     _bufferCapacity;  //!< Internal buffer capacity
  int     _bufferSize;
  int     _byteCount;       //!< Internal buffer usage byte count
  int     _bitCount;
  imgpel* _streamBuffer;    //!< Internal buffer
};

#endif // ENCODER_INC_BITSTREAM_H

