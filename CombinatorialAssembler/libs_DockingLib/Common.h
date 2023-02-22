#ifndef COMMON_H
#define COMMON_H

#include <Surface.h>
#include <Logger.h>

#include "ChemMolecule.h"

#include <fstream>

class Common {
public:

  static void checkFile(const char* fileName, std::ifstream& file);
  static void readChemMolecule(const std::string fileName, ChemMolecule& mol, const ChemLib& protLib);
  static void readChemMolecule(const std::string fileName, ChemMolecule& mol);
  static void readHydrogenMolecule(const std::string fileName, Molecule<Atom>& mol);
  static void readSurface(const std::string fileName, Surface& surface);
};

#endif //COMMON_H
