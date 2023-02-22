#ifndef _LOGGER_H
#define _LOGGER_H

#include "DebugStream.h"

#include <string>

/*
  CLASS
     Logger

  A singleton utility class that is used to print messages to a log
  file. It is based on the DebugStream class. This means that a level
  is attached to each message passed to the logger and the extent of
  the information actually printed by the logger is controlled by the
  log-level parameter (defined in the program configuratin
  file). Specifically, the message will be printed to the log file
  only if its level is smaller then log-level.

  The user should note that the class is not thread safe.

  KEYWORDS
   Log, DebugStream

  AUTHORS
  Oranit Dror (mailto: oranit@post.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 2004.
CHANGES LOG
<UL>
<LI>9.01.05 Dina: add setLevel method
</LI>
</UL>

  USAGE
  The Logger is unitialized with two parameters logLevel and logFile,
  which are defined in the program configuration file (i.e. defined by
  the "Parameters" singleton class).

  The syntax for printing a message is one of the following:
  Logger::getInstance()(Logger::MSG_LEVEL) << msg;
  Logger::errorMessage() << msg;
  Logger::warningMessage() << msg;
  Logger::infoMessage() << msg;
  Logger::debugMessage() << msg;
*/
class Logger {
public:
  ////
  // The different message levels
  enum MSG_LEVEL {ERROR_LEVEL = 0,
		  WARNING_LEVEL = 1,
		  INFO_LEVEL = 2,
		  DEBUG_LEVEL = 3};


  ////
  static Logger& getInstance();

  ////
  // Gets a reference to the Logger and starts a new error message
  // One should note that using this method and only then streaming
  // the message is not thread safe.
  static Logger& errorMessage() {
    return getInstance()(ERROR_LEVEL);
  }

  ////
  // Gets a reference to the Logger and starts a new warning message
  // One should note that using this method and only then streaming
  // the message is not thread safe.
  static Logger& warningMessage() {
    return getInstance()(WARNING_LEVEL);
  }

  ////
  // Gets a reference to the Logger and starts a new info message
  // One should note that using this method and only then streaming
  // the message is not thread safe.
  static Logger& infoMessage() {
    return getInstance()(INFO_LEVEL);
  }

  ////
  // Gets a reference to the Logger and starts a new debug message
  // One should note that using this method and only then streaming
  // the message is not thread safe.
  static Logger& debugMessage() {
    return getInstance()(DEBUG_LEVEL);
  }

  ////
  // Returns true if messages with the given level will be printed to
  // the log file
  static bool isActiveLevel(MSG_LEVEL level) {
    if (level <= getInstance().getLevel()) {
      return true;
    }

    return false;
  }

  ////
  MSG_LEVEL getLevel() {
    return (MSG_LEVEL) log.getDebugLevel();
  }

  //// set level
  void setLevel(MSG_LEVEL msgLevel) {
    log.setDebugLevel(msgLevel);
  }

  ////
  // Starts a new message and specifies its level.
  // One should note that setting the level outside of the actual
  // message (i.e., first calling this operator and only then
  // streaming the message) is not thread safe.
  Logger& operator()(MSG_LEVEL msgLevel);


  ////
  // If a message was started with a low enough message level then
  // given data will be printed to the log file
  template<class T>
  Logger& operator<<(const T& data);

  ////
  // If a message was started with a low enough message level then
  // given data will be printed to the log file
  template<class T>
  Logger& operator<<(T& data);

  ////
  // The Logger can stream a function that accepts an ostream and
  // returns one.
  Logger& operator<<(std::ostream& (*__pf)(std::ostream&));

  ////
  ~Logger();

private:

  ////
  static void init();


  ////
  static Logger* instance;

  ////
  DebugStream log;


  ////
  Logger(const std::string logFile, unsigned int logLevel);
};


/**************************************
 *  Template Function Implementation  *
 **************************************/

template<class T>
Logger& Logger::operator<<(const T& data) {
  log << data;
  return *this;
}


template<class T>
Logger& Logger::operator<<(T& data) {
  log << data;
  return *this;
}



#endif //_LOGGER_H
