#ifndef CHEM_LIB
#define CHEM_LIB

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

/*
CLASS
  ChemEntry

  This is a class that keeps protein atom attributes: type, radius and charge

KEYWORD
  ChemEntry, ChemAtom radius, charge, type

AUTHORS
  Dina Schneidman (duhovka@tau.ac.il)

  copyright: GAMBA Group , Tel-Aviv Univ. Israel, 2004.

GOALS
  ChemEntry keeps atom attribures that are not given in PDB file.

CHANGES LOG
<UL>
</UL>

USAGE
*/
class ChemEntry {
 public:
  // GROUP: Inspectors
  ////
  bool isSameEntry(const std::string residueName, const std::string atomName) const;

  //// returns radius
  float getEntryRadius() const { return radius_; }

  //// returns charge
  float getEntryCharge() const { return charge_; }

  //// returns atom type
  int getEntryChemType() const { return chemType_; }

  //// read entry
  friend std::istream& operator>>(std::istream& s, ChemEntry &entry) {
    return s >> entry.residueName_ >> entry.atomName_ >> entry.radius_ >> entry.charge_ >> entry.chemType_;
  }

  //// output entry
  friend std::ostream& operator<<(std::ostream& s, const ChemEntry &entry) {
    s << entry.residueName_ << ' ' << entry.atomName_  << ' '
      << entry.radius_ << ' ' << entry.charge_ << ' ' <<  entry.chemType_;
    return s;
  }

 public:
  std::string residueName_;
  std::string atomName_;
  float radius_;
  float charge_;
  int chemType_;
};


/*
CLASS
  ChemLib

  This is a class that serves as a library for protein attributes:
  type, radius and charge.

KEYWORD
  ChemEntry, ChemAtom, ChemLib radius, charge, type

AUTHORS
  Dina Schneidman (duhovka@tau.ac.il)

  copyright: GAMBA Group , Tel-Aviv Univ. Israel, 2004.

GOALS
  ChemLib keeps attribures for each of the 20 amino acids atoms. These
  attribures are not given in PDB file.

CHANGES LOG
<UL>
</UL>

USAGE
 The ChemLib is initialized with the special library file (mol/lib/gamb++/chem.lib).

*/
class ChemLib {
 public:
  // GROUP: Constructors

  //// Init and load library
  ChemLib(const std::string libFileName);
  void loadLibrary(const std::string libFileName);

  // GROUP: Inspectors

  //// get radius
  float getAtomRadius(const char* residueName, const char* atomPDBType) const;

  //// get charge
  float getAtomCharge(const char* residueName, const char* atomPDBType) const;

  //// get chemical type
  int getAtomChemType(const char* residueName, const char* atomPDBType) const;

  //// print function
  void printLibrary(std::ostream& outFile) const;

  //// get the whole lib entry
  const ChemEntry* getLibEntry(const char* residueName, const char* atomPDBType) const;

 private:
  const ChemEntry* searchLibEntry(const std::string residueName, const std::string atomName) const;

 private:
  std::vector<ChemEntry> libVector;
};

#endif
