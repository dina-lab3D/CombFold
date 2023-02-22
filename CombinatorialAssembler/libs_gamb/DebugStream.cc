#include "DebugStream.h"

DebugStream::DebugStream(std::ostream& outStream) :
  out(&outStream), dbgLevel(0), withTimer(false), withNewLine(false),
  withLevel(false), withIndent(false), timer(), active(false),
  ownStream(false) {}

DebugStream::DebugStream(const char *filename):
  out(&std::cerr), /* temporary initilization for a reference can't be blank */
  dbgLevel(0), withTimer(false), withNewLine(false),
  withLevel(false), withIndent(false), timer(), active(false),
  ownStream(false) {
  std::ofstream *outstream=new std::ofstream(filename);
  if (!(*outstream)) std::cerr << "Can't open file : "
                               << filename
                               << " as DebugStream."
                               << std::endl;
  else {
    ownStream = true;
    out=outstream;
  }
}

DebugStream::~DebugStream()
{
    if (ownStream){
        //own the output stream - has to close it.
      ((std::ofstream *)out)->close();
        delete out;
    }
}

void DebugStream::setDebugLevel(unsigned int debugLevel)
{
  dbgLevel = debugLevel;
}

void DebugStream::newLines(const bool print)
{
  withNewLine = print;
}

void DebugStream::messageLevels(const bool print)
{
  withLevel = print;
}

void DebugStream::messageTimes(const bool print)
{
  withTimer = print;
}

void DebugStream::indent(const bool print)
{
  withIndent = print;
}

void DebugStream::resetTimer()
{
  timer.reset();
}

void DebugStream::setFloatPrecision(unsigned int prec)
{
  (*out).setf(std::ios::fixed, std::ios::floatfield);
  (*out).precision(prec);
}

void DebugStream::setWidth(unsigned int perc)
{
    (*out).width(perc);
}


void DebugStream::setWidth(unsigned int numOfChars,
			   unsigned int adjustment) {
  if (adjustment == ADJUSTMENT_LEFT) {
    (*out).setf(std::iostream::left, std::iostream::adjustfield);
  } else {
    (*out).setf(std::iostream::right, std::iostream::adjustfield);
  }

  (*out).width(numOfChars);
}

void DebugStream::autoFlush(const bool autoflush)
{
  if (autoflush)
    out->setf(std::ios::unitbuf);
  else
    out->unsetf(std::ios::unitbuf);
}

DebugStream& DebugStream::operator()(unsigned int messageLevel)
{
  active = (messageLevel <= dbgLevel);

  if (active) {
    if (withNewLine)
      *out << '\n';
    if (withLevel) {
      out->width(3);
      *out << messageLevel;
    }
    if (withTimer) {
      if (withLevel)
        *out << ' ';
      *out << timer << '|' << ' ';
    }
    else if (withLevel)
      *out << ':' << ' ';
    if (withIndent)
      for (unsigned int i=0; i<messageLevel/10; ++i)
        *out << ' ' << ' ';
  }
  return *this;
}

DebugStream& DebugStream::operator<<(const char* const data)
{
  if (active)
    *out << data;

  return *this;
}


DebugStream& DebugStream::operator<<(std::ostream& (*__pf)(std::ostream&)) {
  if (active) {
    *out << *__pf;
  }

  return *this;

}
