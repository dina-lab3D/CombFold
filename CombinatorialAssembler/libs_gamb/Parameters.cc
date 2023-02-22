#include <iostream>
#include <sstream>
#include <vector>
#include "Parameters.h"
#include <boost/algorithm/string.hpp>

typedef Parameters::ParamType ParamType;

class Parameter
{
public:

  Parameter();
  Parameter(const std::string& name, const int val);
  Parameter(const std::string& name, const float val);
  Parameter(const std::string& name, const std::string& val);
  Parameter(const Parameter& param);

  void dump(FILE* outputFile) const;

  ~Parameter() {}
  const std::string& paramLabel() const;
  ParamType paramType() const;

  int intValue() const;
  float floatValue() const;
  const std::string& stringValue() const;

  Parameter& operator=(const Parameter& param);

  friend bool operator<(const Parameter& p, const Parameter& q);
  friend std::ostream& operator<<(std::ostream& out, const Parameter& p);

private:
  std::string paramName;
  ParamType type;
    union {
      int i;
      float f;
    };
  std::string s;
};

typedef std::vector<Parameter> ParamList;

static ParamList paramList;

Parameter::Parameter() : paramName(), type(Parameters::Undef)
{}

Parameter::Parameter(const std::string& name, const int val)
{
  paramName = name;
  i = val;
  type = Parameters::Int;
}

Parameter::Parameter(const std::string& name, const float val)
{
  paramName = name;
  f = val;
  type = Parameters::Float;
}

Parameter::Parameter(const std::string& name, const std::string& val)
{
  paramName = name;
  s = val;
  type = Parameters::Str;
}

Parameter::Parameter(const Parameter& param)
{
  paramName = param.paramName;
  type = param.type;
  if (type == Parameters::Int)
    i = param.i;
  else
    f = param.f;
  s = param.s;
}

const std::string& Parameter::paramLabel() const
{
  return paramName;
}

ParamType Parameter::paramType() const
{
  return type;
}

int Parameter::intValue() const
{
  return i;
}

float Parameter::floatValue() const
{
  return f;
}

const std::string& Parameter::stringValue() const
{
  return s;
}

Parameter& Parameter::operator=(const Parameter& param)
{
  paramName = param.paramName;
  type = param.type;
  if (type == Parameters::Int)
    i = param.i;
  else
    f = param.f;
  s = param.s;
  return *this;
}

bool operator<(const Parameter& p, const Parameter& q)
{
  return (p.paramName < q.paramName);
}

std::ostream& operator<<(std::ostream& out, const Parameter& p) {
  switch(p.type) {
  case Parameters::Int:
    return out << p.paramName << '\t' << "(Int)" << '\t' << p.i;
  case Parameters::Float:
    return out << p.paramName << '\t' << "(Float)" << '\t' << p.f;
  case Parameters::Str:
    return out << p.paramName << '\t' << "(Str)" << '\t' << p.s;
  case Parameters::Undef:
    break;
  }
  return out << '\n';
}


void Parameter::dump(FILE* outputFile) const {
  switch(type) {
  case Parameters::Int:
    fprintf(outputFile, "%s = %d", paramName.c_str(), i);
  case Parameters::Float:
    fprintf(outputFile, "%s = %f", paramName.c_str(), f);
  case Parameters::Str:
    fprintf(outputFile, "%s = %s", paramName.c_str(), s.c_str());
  case Parameters::Undef:
    break;
  }
}


ParamList::iterator findInsertionPoint(ParamList& paramList,
                                       const std::string& paramName)
{
  unsigned short start = 0;
  unsigned short stop = (unsigned short)paramList.size();
  while (stop != start) {
    unsigned short pos = start + (stop-start)/2;
    int comp = paramName.compare(paramList[pos].paramLabel());
    if (comp == 0)
      stop = start = pos;
    else if (comp > 0)
      start = pos + 1;
    else
      stop = pos;
  }

  ParamList::iterator it=paramList.begin();
  it+=stop;
  return it;
}

Parameters::Parameters()
{}

void Parameters::readParameters(std::istream& paramStream)
{
  while (!paramStream.eof()) {
    std::string param;
    getline(paramStream, param);
    boost::trim(param);
    std::string paramName = nextToken(param);

    if (paramName.length() == 0) continue;

    if (*(paramName.data()) == '#') continue;

    updateParameter(paramName, param.c_str());
  }
}


bool Parameters::empty() {
  return paramList.empty();
}

void Parameters::addParameter(const std::string& paramName, const int value)
{
  ParamList::iterator pos = findInsertionPoint(paramList, paramName);
  if (pos != paramList.end() && (*pos).paramLabel() == paramName)
    (*pos) = Parameter(paramName, value);
  else
    paramList.insert(pos, Parameter(paramName, value));
}

void Parameters::addParameter(const std::string& paramName, const double value)
{
  ParamList::iterator pos = findInsertionPoint(paramList, paramName);
  if (pos != paramList.end() && (*pos).paramLabel() == paramName)
    (*pos) = Parameter(paramName, (float)value);
  else
    paramList.insert(pos, Parameter(paramName, (float)value));
}

void Parameters::addParameter(const std::string& paramName, const std::string& value)
{
  ParamList::iterator pos = findInsertionPoint(paramList, paramName);
  if (pos != paramList.end() && (*pos).paramLabel() == paramName)
    (*pos) = Parameter(paramName, value);
  else
    paramList.insert(pos, Parameter(paramName, value));
}

void Parameters::updateParameter(const std::string& paramName,
                                 const char* const value)
{
  ParamList::iterator pos = findInsertionPoint(paramList, paramName);
  if (pos != paramList.end() && (*pos).paramLabel() == paramName)
    switch ((*pos).paramType()) {
    case Int:
      (*pos) = Parameter(paramName, atoi(value));
      break;
    case Float:
      (*pos) = Parameter(paramName, (float)atof(value));
      break;
    case Str:
      (*pos) = Parameter(paramName, std::string(value));
    case Undef:
      (*pos) = Parameter(paramName, std::string(value));
    }
  else
    paramList.insert(pos, Parameter(paramName, std::string(value)));
}

ParamType Parameters::paramType(const std::string& paramName)
{
  ParamList::iterator pos = findInsertionPoint(paramList, paramName);
  if (pos != paramList.end() && (*pos).paramLabel() == paramName)
    return (*pos).paramType();
  else
    return Undef;
}


int Parameters::getInt(const std::string& paramName, const int& defaultValue)
{
  ParamList::iterator pos = findInsertionPoint(paramList, paramName);
  if (pos != paramList.end() && (*pos).paramLabel() == paramName)
    switch ((*pos).paramType()) {
    case Int:
      return (*pos).intValue();
    case Float:
      return (int)(*pos).floatValue();
    case Str:
      return atoi((*pos).stringValue().data());
    case Undef:
      break;
    }
  return defaultValue;
}

float Parameters::getFloat(const std::string& paramName, const float& defaultValue)
{
  ParamList::iterator pos = findInsertionPoint(paramList, paramName);
  if (pos != paramList.end() && (*pos).paramLabel() == paramName)
    switch ((*pos).paramType()) {
    case Float:
      return (*pos).floatValue();
    case Int:
      return (float)(*pos).intValue();
    case Str:
      return (float)atof((*pos).stringValue().data());
    case Undef:
      break;
    }
  return defaultValue;
}

std::string Parameters::getString(const std::string& paramName, const std::string& defaultValue)
{
  ParamList::iterator pos = findInsertionPoint(paramList, paramName);
  if (pos != paramList.end() && (*pos).paramLabel() == paramName)
    switch ((*pos).paramType()) {
    case Str:
      return (*pos).stringValue();
    case Float: {
      return std::to_string((*pos).floatValue());
    }
    case Int: {
      return std::to_string((*pos).intValue());
    }
    case Undef:
      break;
    }
  return defaultValue;
}

void Parameters::dump(std::ostream& out)
{
  for (ParamList::iterator i=paramList.begin(); i != paramList.end(); ++i)
    out << *i << '\n';
}

void Parameters::dump(DebugStream& out, const unsigned int msgLevel)
{
  for (ParamList::iterator i=paramList.begin(); i != paramList.end(); ++i)
    out(msgLevel) << *i;
}

void Parameters::dump(FILE* outputFile) {
  for (ParamList::iterator i = paramList.begin() ; i != paramList.end() ; i++) {
    i->dump(outputFile);
    fprintf(outputFile, "\n");
  }

  fprintf(outputFile, "\n");
}

std::string Parameters::nextToken(std::string& str)
{
  unsigned int start = 0;
  while (start < str.length() &&
         (str[start] == ' ' || str[start] == '\t' || str[start] == '\n'))
    ++start;

  if (start >= str.length()) {
    str = "";
    return "";
  }

  unsigned int stop = start+1;
  while (stop < str.length() &&
         str[stop] != ' ' && str[stop] != '\t' && str[stop] != '\n')
    ++stop;

  unsigned int next = stop;
  while (next < str.length() &&
         (str[next] == ' ' || str[next] == '\t' || str[next] == '\n'))
    ++next;

  std::string result = str.substr((int)start, stop-start);
  str = str.substr((int)next);
  return result;
}
