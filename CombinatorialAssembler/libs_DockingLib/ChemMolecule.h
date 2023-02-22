#ifndef CHEM_MOL_H
#define CHEM_MOL_H

#include "ChemAtom.h"
#include "ResidueIndexedMolecule.h"
#include "ChemLib.h"
#include "PDB.h"
#include "Logger.h"
#include "EnergyAtom.h"

#include <iostream>
#include <fstream>
#include <set>

/*
CLASS
  ChemMolecule

  This is a class that is a container of ChemAtoms

KEYWORD
  ChemAtom, radius, charge, type

AUTHORS
  Dina Schneidman (duhovka@tau.ac.il)

  copyright: GAMBA Group , Tel-Aviv Univ. Israel, 2004.

GOALS
  ChemMolecule is a molecule of ChemAtoms. ChemAtom extends Atom and has 4
  additional attributes: chemical type, radius, charge and probability of
  the atom of being steady. The type, radius and charge are taken from
  ChemLib. This library is initiated with chem.lib file.
  The chemical types are from:
  Zhang,C., Cornette,J.L. and DeLisi,C. (1997) Determination of atomic
  desolvation energies from the structures of crystallized
  proteins. J. Mol. Biol., 267, 707 726.
  This types are only for protein atoms!
  The radii are taken from CHARMM, however they are slightly bigger
  to account for missing hydrogens. The charges are from CHARMM too.
  The probabilities of being steady were computed based on the maximal
  distance between C-alpha and side-chain atoms of that residue. The
  lower the distance the higher is the prob. of being steady.
CHANGES LOG
<UL>
</UL>

USAGE
  ChemMolecule chemMolecule;
  ifstream pdbFile("mol.pdb");
  ChemLib chemLib("chem.lib");
  chemMolecule.loadMolecule(pdbFile, chemLib);
*/
class ChemMolecule : public ResidueIndexedMolecule<ChemAtom> {
private:
  static const float flexProbs[];
  static const float maxCaDist[];

public:
  //// load molecule from PDB file, using chemLib for attributes.
  void loadMolecule(std::istream &molFile, const ChemLib& chemLib, const PDB::Selector& selector = PDB::WaterHydrogenUnSelector());

  //// reads and assigns default values for radius, charge and chemType
  void readAllPDBfile(std::istream &molFile, const PDB::Selector& selector = PDB::WaterHydrogenUnSelector());

  //// compute probabilities of being steady
  void computeSteadyProbs();

  //// read and save binding site, returns number of residues in site
  unsigned int readBindingSite(const std::string fileName);

  //// read and save blocking site, returns number of residues in site
  unsigned int readBlockingSite(const std::string fileName);

  ////
  const std::set<unsigned int>& getBindingSite(unsigned int siteNumber = 0) const;

  ////
  const std::set<unsigned int>& getBlockingSite(unsigned int siteNumber = 0) const;

  ////
  const ChemAtom& getChemAtom(int atomIndex) const;

  ////
  const Vector3& getCentroid() const { return centroid_; }

  bool isProtein() const { return isProtein_; }

  //// get the volume by summing the volume of all atoms
  double getVolume() const;

  typedef std::pair<int, int> ResidueRange;
  typedef std::pair<char, ResidueRange> Fragment;

  //// get the fragments
  const std::vector<Fragment>& getFragments() const;

  //// add mol2 atom type
  void addMol2Type();

  //// compute accessible surface area
  float computeASAperAtom();
private:
  //// read site file and return number of reasidues that were read
  unsigned int readSiteFile(const std::string fileName, std::set<unsigned int>& site);

  //// checks if the molecule is protein (if there is at least one Ca atom)
  void setIsProtein();

  //// compute fragments of the same chain
  void computeFragments();

protected:
  //// entries of binding site residues according to ResidueIndexedMolecule
  std::vector<std::set<unsigned int> > bindingSites_;

  //// entries of blocking site residues according to ResidueIndexedMolecule
  std::vector<std::set<unsigned int> > blockingSites_;

  ////
  Vector3 centroid_;

  ////
  bool isProtein_;

  //// fragments of the same chain in the molecule
  std::vector<Fragment> fragments_;
};

bool compareFragment(ChemMolecule::Fragment f1, ChemMolecule::Fragment f2);

#endif
