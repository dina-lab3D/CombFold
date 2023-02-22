#include "ChemLib.h"

#include <boost/algorithm/string.hpp>

bool ChemEntry::isSameEntry(const std::string residueName, const std::string atomName) const {
  if(residueName == residueName_ && atomName == atomName_)
    return true;
  return false;
}

ChemLib::ChemLib(const std::string libFileName) {
  loadLibrary(libFileName);
}

void ChemLib::loadLibrary(const std::string libFileName) {
  std::ifstream libFile(libFileName);
  if(!libFile) {
    std::cerr << "Can't find library file: " << libFileName << std::endl;
    exit(1);
  }

  std::string line;
  while (!libFile.eof()) {
    getline(libFile, line);
    boost::trim(line); // remove all spaces
    // skip comments
    if (line[0] == '#' || line[0] == '\0') continue;
    std::vector < std::string > splitResults;
    boost::split(splitResults, line, boost::is_any_of(" "), boost::token_compress_on);
    for(unsigned int i=0; i<splitResults.size(); i++) boost::trim(splitResults[i]);

    ChemEntry entry;

    entry.residueName_ = splitResults[0];
    entry.atomName_  = splitResults[1];
    entry.radius_ = atof(splitResults[2].c_str());
    entry.charge_ = atof(splitResults[3].c_str());
    entry.chemType_ = atoi(splitResults[4].c_str());
    libVector.push_back(entry);

  }
}

const ChemEntry* ChemLib::getLibEntry(const char* residueName, const char* atomPDBType) const {
  std::string resName(residueName);
  std::string atomName(atomPDBType);

  if(resName.length() != 3) {
    std::cerr << "Error in resName: " << resName << " length " << resName.length() << std::endl;
    //    exit(1);
  }

  if(atomName.length() != 4) {
    std::cerr << "Error in atomName: " << atomName << " length " << atomName.length() << std::endl;
    //    exit(1);
  }

  std::string atomTrimmedName = atomName;
  boost::trim(atomTrimmedName);

  const ChemEntry* libEntry = searchLibEntry(resName, atomTrimmedName);
  if(libEntry == NULL) {
    libEntry = searchLibEntry("XXX", atomTrimmedName);
    if(libEntry == NULL) {
      switch(atomTrimmedName.length()) {
      case 1: atomTrimmedName.append("??"); break;
      case 2: atomTrimmedName.append("?"); break;
      case 3: atomTrimmedName[2] = '?'; break;
      default: atomTrimmedName = "???"; break;
      }
      libEntry = searchLibEntry("XXX", atomTrimmedName);

      if(libEntry == NULL) {
	atomTrimmedName[1] = '?';
	libEntry = searchLibEntry("XXX", atomTrimmedName);
	if(libEntry == NULL)
          std::cerr << "Unidentified atom " << atomName << " using default values!!!" << std::endl;
      }
    }
  }
  return libEntry;
}

const ChemEntry* ChemLib::searchLibEntry(const std::string residueName, const std::string atomName) const {
  const ChemEntry* libEntry = NULL;
  for(std::vector<ChemEntry>::const_iterator iter = libVector.begin();
      iter != libVector.end(); iter++)
    if((*iter).isSameEntry(residueName, atomName)) {
      libEntry = &(*iter);
      break;
    }
  return libEntry;
}

float ChemLib::getAtomRadius(const char* residueName, const char* atomPDBType) const {
  const ChemEntry* libEntry = getLibEntry(residueName, atomPDBType);
  if(libEntry!=NULL)
    return libEntry->getEntryRadius();
  return 1.5;
}

float ChemLib::getAtomCharge(const char* residueName, const char* atomPDBType) const {
  const ChemEntry* libEntry = getLibEntry(residueName, atomPDBType);
  if(libEntry!=NULL)
    return libEntry->getEntryCharge();
  return 0.0;
}

int ChemLib::getAtomChemType(const char* residueName, const char* atomPDBType) const {
  const ChemEntry* libEntry = getLibEntry(residueName, atomPDBType);
  if(libEntry!=NULL)
    return libEntry->getEntryChemType();
  return 0;
}

void ChemLib::printLibrary(std::ostream& outFile) const {
  for(std::vector<ChemEntry>::const_iterator iter = libVector.begin();
      iter != libVector.end(); iter++) {
    outFile << *iter << std::endl;
  }
}
