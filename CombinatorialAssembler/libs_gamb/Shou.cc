#include <stdlib.h>
#include <string>
#include "Shou.h"


bool Shou::isRec(const char* const shouRec)
{
  return (atoi(shouRec) != 0);
}

unsigned short Shou::atomIndex(const char *const shouRec,
			       unsigned int atom)
{
  unsigned short fieldLength=5;
  std::string num;

  num.resize(fieldLength + 1);
  for(int i=0; i<fieldLength; i++)
    num[i] = shouRec[5*atom+i];
  num[fieldLength]='\0';
  return (unsigned short)atoi(num.c_str());
}

float Shou::xCoord(const char* const shouRec)
{
  unsigned short fieldLength=yCoordField-xCoordField;
  std::string num;

  num.resize(fieldLength + 1);
  for(int i=0; i<fieldLength; i++)
    num[i] = shouRec[xCoordField+i];
  num[fieldLength]='\0';
  return (float)atof(num.c_str());
}

float Shou::yCoord(const char* const shouRec)
{
  unsigned short fieldLength=zCoordField-yCoordField;
  std::string num;

  num.resize(fieldLength + 1);
  for(int i=0; i<fieldLength; i++)
    num[i] = shouRec[yCoordField+i];
  num[fieldLength]='\0';
  return (float)atof(num.c_str());
}

float Shou::zCoord(const char* const shouRec)
{
  unsigned short fieldLength=surfaceField-zCoordField;
  std::string num;

  num.resize(fieldLength + 1);
  for(int i=0; i<fieldLength; i++)
    num[i] = shouRec[zCoordField+i];
  num[fieldLength]='\0';
  return (float)atof(num.c_str());
}

float Shou::surface(const char* const shouRec)
{
  unsigned short fieldLength=xNormField-surfaceField;
  std::string num;

  num.resize(fieldLength + 1);
  for(int i=0; i<fieldLength; i++)
    num[i] = shouRec[surfaceField+i];
  num[fieldLength]='\0';
  return (float)atof(num.c_str());
}

float Shou::xNorm(const char* const shouRec)
{
  unsigned short fieldLength=yNormField-xNormField;
  std::string num;

  num.resize(fieldLength + 1);
  for(int i=0; i<fieldLength; i++)
    num[i] = shouRec[xNormField+i];
  num[fieldLength]='\0';
  return (float)atof(num.c_str());
}

float Shou::yNorm(const char* const shouRec)
{
  unsigned short fieldLength=zNormField-yNormField;
  std::string num;

  num.resize(fieldLength + 1);
  for(int i=0; i<fieldLength; i++)
    num[i] = shouRec[yNormField+i];
  num[fieldLength]='\0';
  return (float)atof(num.c_str());
}

float Shou::zNorm(const char* const shouRec)
{
  unsigned short fieldLength=zNormField-yNormField;
  std::string num;

  num.resize(fieldLength + 1);
  for(int i=0; i<fieldLength; i++)
    num[i] = shouRec[zNormField+i];
  num[fieldLength]='\0';
  return (float)atof(num.c_str());
}
