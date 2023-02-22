#include "Logger.h"

#include "macros.h"
#include <iostream>

#include "Parameters.h"

Logger* Logger::instance = NULL;

Logger& Logger::getInstance() {
  if (instance == NULL) {
    init();
  }

  return *instance;
}


void Logger::init() {
  std::string logFile = Parameters::getString("log-file", "log");
  unsigned int logLevel = Parameters::getInt("log-level", WARNING_LEVEL);
  instance = new Logger(logFile, logLevel);
}

Logger::Logger(const std::string logFile, unsigned int logLevel) : log(logFile.c_str()) {
  log.autoFlush(true);
  log.setDebugLevel(logLevel);
  log.newLines(false);
  log.messageTimes(true);
}


Logger& Logger::operator() (MSG_LEVEL msgLevel) {
  std::string msgHeader = "";
  switch(msgLevel) {
  case ERROR_LEVEL:
    msgHeader = "ERROR: ";
    break;
  case WARNING_LEVEL:
    msgHeader = "WARNING: ";
    break;
  case INFO_LEVEL:
    msgHeader = "INFO: ";
    break;
  case DEBUG_LEVEL:
    msgHeader = "DEBUG: ";
    break;
  }

  log(msgLevel) << msgHeader;
  return *this;
}


Logger& Logger::operator<<(std::ostream& (*__pf)(std::ostream&)) {
  log << *__pf;
  return *this;

}

Logger::~Logger() {
  //  delete log;
  delete instance;
  instance = NULL;
}
