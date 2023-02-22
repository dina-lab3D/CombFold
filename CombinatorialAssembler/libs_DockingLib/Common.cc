#include "Common.h"

void Common::checkFile(const char* fileName, std::ifstream& file) {
  if(!file) {
    std::cerr << "Can't open file: " << fileName << std::endl;
    Logger::errorMessage() << "Can't open file: " << fileName << std::endl;
    exit(1);
  }
}

void Common::readChemMolecule(const std::string fileName, ChemMolecule& mol, const ChemLib& protLib) {
  std::ifstream molFile(fileName.c_str());
  checkFile(fileName.c_str(), molFile);
  mol.loadMolecule(molFile, protLib);
  molFile.close();
  Logger::infoMessage() << "Chem molecule: " << mol.size() << " atoms were read" << std::endl;
}

void Common::readChemMolecule(const std::string fileName, ChemMolecule& mol) {
  std::ifstream molFile(fileName.c_str());
  checkFile(fileName.c_str(), molFile);
  mol.readPDBfile(molFile, PDB::WaterHydrogenUnSelector());
  molFile.close();
  Logger::infoMessage() << "Chem molecule: " << mol.size() << " atoms were read" << std::endl;
}

void Common::readHydrogenMolecule(const std::string fileName, Molecule<Atom>& mol) {
  std::ifstream molFile(fileName.c_str());
  checkFile(fileName.c_str(), molFile);
  mol.readAllPDBfile(molFile, PDB::WaterUnSelector());
  molFile.close();
  Logger::infoMessage() << "Hydrogen molecule: " << mol.size() << " atoms were read" << std::endl;
}

void Common::readSurface(const std::string fileName, Surface& surface) {
  std::ifstream surfFile(fileName.c_str());
  checkFile(fileName.c_str(), surfFile);
  surface.readShouFile(surfFile);
  surfFile.close();
  Logger::infoMessage() << "Surface: " << surface.size() << " surface points were read" << std::endl;
}
