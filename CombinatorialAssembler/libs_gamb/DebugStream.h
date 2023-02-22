#ifndef _DebugStream_h
#define _DebugStream_h

#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include "Timer.h"

/*
CLASS
  DebugStream

  Defines a class for printing logs and debug information. User may determine
  the extent of the debugging information actually printed by setting the debug
  level.

KEYWORDS
  debug

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI> Defining a streamer operator that supports printing endl
16/02/04 Oranit Dror
</LI>
<LI> Using unitbuf flag of Stream for autoflash instead of keeping our own flag
15/02/04 Dina Duhovny and Oranit Dror
</LI>
<LI> Added two functions that allows the usage of percision both regular and
in fixed floating point output
30/8/99 Zipi Fligelman
</LI>
<LI> Added the mode of the DebugStream owning it's ostream. This enables the
 user to construct the DebugStream with a file name and open/close the file
with the DebugStream object liveness.
5/9/99  Ram Nathaniel
<LI> Added the possibility to autoFlush - to make sure that debug messages are
not stuck in the buffers when the program crushes.
5/9/99 Ram Nathaniel
</UL

GOALS
  The DebugStream class was written for the purpose of debugging and generating
  run-time logs. The extent of the information actually printed by the stream
  is controlled by the debug level parameter. A level is attached to each
  message passed to the stream. If this level is smaller then the debug level
  the message is printed.

USAGE
  Using the DebugStream is quite straight-forward. The DebugStream use is
  similar to a regular output stream with a minor difference. A message level
  must be attached to each message. The message level is set by using
  operator() before commencing a sequence of operator<<.

  The following program will output all messages of level 10 and under.
  EXAMPLE
    DebugStream log(cout);
    log.setDebugLevel(10);
    log(10) << "This is a level " << 10 << " message";
    log << " and so is this";
    log(15) << "But this is a level" << 15 << " message and will not be shown";
            << " because the DebugLevel is lower";
    log() << "This is a level" << 0 << " message. You should definitely see it";
  END

  The debug stream can also be initialized by a file name. In this case the
  DebugStream will open the file when created and close it when done.

  The DebugStream class allows the user to control the format of the stream
  output. The user may choose to output the user time, a new-line after every
  message or the message level itself.

  If you wish the DebugStream not to use writing buffers use the autoFlush.
  This option may make the DebugStream work a little slower but will ensure
  that the data was written to the file before the next command was preformed.
*/
class DebugStream
{
public:
  enum{ADJUSTMENT_LEFT, ADJUSTMENT_RIGHT};

  //// Contructor: Intialize the DebugStream class with a true output stream
  // output that the DebugStream decides should be printed will be forwarded
  // to this output stream.
  explicit DebugStream(std::ostream& outStream = std::cout);

  //// Constructor: Initialize the DebugStream class with the name of the
  // output file. This way the class can be initilized when declared as a
  // global varialb. If file is unavailable DebugStream will use cerr and
  // will notify to the problem.
  explicit DebugStream(const char *filename);

  ////destructor: mainly closes the file if ownStream is true.
  ~DebugStream();

  //// Set the debug level controlling the extent of the debugging messages
  // seen. The higher the level the more messages the user will actually see.
  // This parameter may be read from a prameters file controlling the extent
  // of information actually shown during run-time.
  void setDebugLevel(const unsigned int debugLevel);

  //// In line function returning the debug level
  // for logging purposes
  inline unsigned getDebugLevel() const;

  //// flushing mechanism like that of the ostream type
  inline void flush();
  //// Set the new-lines option. If true, starting a new message with
  // operator() causes a new-line to be printed.
  void newLines(const bool print = true);

  //// Set the message levels option. If true, for each mesasage the message
  // level is shown at the beginning. Set to false by default.
  void messageLevels(const bool print = true);

  //// Set the timer option. If true the time of each message is shown at the
  // head of the message. The time show is the system's user time and not the
  // real clocked time.
  void messageTimes(const bool print = true);

  //// Set the indent option. If true log messages are indented according to
  // message levels. 2 spaces per 10 debug levels.
  void indent(const bool print = true);

  //// Resets timer to 0.
  void resetTimer();

  //// Useful Function for working with a specific floating point percision
  void setFloatPrecision(unsigned int prec);

  //// another useful funciton for indenting
  void setWidth(unsigned int prec);

  ////
  // A function for indentation:
  // Params: numOfChars - Sets the width of the filed
  //         adjustment - Defines the adjustment of the charecters
  //         within the field (possible options: left or right)
  void setWidth(unsigned int numOfChars, unsigned int adjustment);

  //// Set the automaticFlush on so that each writing to the stream will be
  // automatially accompanied by a flush.
  void autoFlush(const bool autoflush = true);

  //// Starts a new message specifying the messages level. If no message level
  // is given then a default message level of 0 is assumed and the message is
  // always shown.
  DebugStream& operator()(unsigned int newMsgLevel = 0);

  //// Output operator. If message was started with a low enough message level
  // then data will be passed on to the output stream with which DebugStream
  // was intialized.
  template<class T>
  DebugStream& operator<<(T& data);

  template<class T>
  DebugStream& operator<<(const T& data);

  // Zipi added to previos version since having problems
  //template<class T*>
  //DebugStream& operator<<(const T*& data);

  // Special overload to const char * const
  DebugStream& operator<<(const char* const data);

  DebugStream& operator<<(std::ostream& (*__pf)(std::ostream&));

private:
  // can't be a reference since it can be an ofstream or an ostream
  //according to the constructor used.
  std::ostream *out;
  unsigned int dbgLevel;
  bool withTimer;
  bool withNewLine;
  bool withLevel;
  bool withIndent;
  Timer timer;
  bool active;

  bool ownStream; //for use if constructed with a file name.
};
/**************************************
 *  Inline Function Implementation    *
 **************************************/
unsigned DebugStream::getDebugLevel() const
{
  return dbgLevel;
}

void DebugStream::flush()
{
  out->flush();
}

/**************************************
 *  Template Function Implementation  *
 **************************************/
template<class T>
DebugStream& DebugStream::operator<<(T& data)
{
  if (active)
    *out << data;

  return *this;
}

template<class T>
DebugStream& DebugStream::operator<<(const T& data)
{
  if (active)
    *out << data;

  return *this;
}


#endif
