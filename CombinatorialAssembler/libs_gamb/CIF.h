#ifndef _CIF_h
#define _CIF_h

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "PDB.h"

/*
CLASS
  CIF

  The CIF class defines a set of tools used for reading the Protein Data Bank
  file format, used to describe molecular structures, like proteins,
  drugs and Nucleic Acids
KEYWORDS
  CIF, protein, Nucleic Acids (RNA, DNA) , record, atom, residue,
  structure, coordinate, field

AUTHORS
  Dina Schneidman

GOALS
  Class CIF is designed to encapsulate a set of methods used for reading the
  CIF file format. Methods are static requiring no object instantiation
  before their use.

USAGE
  Similar to PDB
  Selectors: uses PDB selectors by converting CIF line to PDB line
*/
class CIF {

public:

  // GROUP: CIF record type determination methods.
  ////
  // fill in the atom_site columns table
  static bool readAtomSiteTable(std::istream& ifile);
  static void writeAtomSiteTable(std::ostream& ofile);

  //
  static bool isAtomSiteLine(const std::string& CIFline);

  ////
  // Returns true if the given record is an ATOM record.
  static bool isATOMrec(const std::string& CIFline);

  ////
  // Returns true if the given record is a HETATM record.
  static bool isHETATMrec(const std::string& CIFline);



  //// Returns an ATOM record's X coordinate.
  static float atomXCoord(const std::string& CIFline);

  //// Returns an ATOM record's Y coordinate.
  static float atomYCoord(const std::string& CIFline);

  //// Returns an ATOM record's Z coordinate.
  static float atomZCoord(const std::string& CIFline);

  //// Returns an ATOM record's atom index number.
  static int atomIndex(const std::string& CIFline);

  //// Returns atomType as a string
  static std::string atomType(const std::string& CIFline);

  //// Returns an ATOM record's chain id.
  static std::string atomChainId(const std::string& CIFline);

  //// Chain id named by author
  static std::string atomChainIdAuthor(const std::string& CIFline);

  //// Returns an ATOM record's residue index.
  static int atomResidueIndex(const std::string& CIFline);

  //// Returns the full string of the Residue Name field
  static std::string const getAtomResidueLongName(const std::string& CIFline);

  //// Return the type of the CIF entry, currently only ATOM and HETATM are identified.
  static AtomEntryType getAtomEntryType(const std::string& CIFline);

  //// Returns the occupancy,  which is a measure of the fraction of molecules in the
  //   crystal in which the current atom actually occupies the specified position
  static float getOccupancy(const std::string& CIFline);

  //// Returns the temperature factor,  which is a measure of how much an atom
  // oscillates or vibrates around the specified position
  static float getTempFactor(const std::string& CIFline);

  static std::string getColumn(const std::string& CIFline, const std::string& key);

  /* TODO:
  // Returns true if the given record is an ATOM record of a C-alpha atom.
  static bool isCaATOMrec(const std::string& CIFline);

  // Returns true if the given record is a record of a nucleic acid residue
  static bool isNucleicRecord(const std::string& CIFline);

  //// returns true if the residue name is of nucleic acid
  static bool isNucleicAcid(const std::string& residueName);

  static bool isCONECTrec(const std::string& CIFline);

  static bool isMODELrec(const std::string& CIFline);

  static bool isENDMDLrec(const std::string& CIFline);

  //// Returns an ATOM record's alternative location indicator.
  static char  atomAltLocIndicator(const std::string& CIFline);

  // code for insertion of residues (looks like not needed for cif)
  static char atomResidueICode(const std::string& CIFline);
  ////
  // Gets a CIF record of an atom and returns a one-letter
  // abbreviation of the type of the residue that the atom belongs
  // to. <br>
  // For a nucleic acid atom returns A for adenine, G for guanine, T
  // for thymine, U for uracil and X for unknown type.
  static char atomResidueType(const std::string& CIFline);

  ////
  // Gets a one-letter abbreviation of the type of the residue that
  // the atom belongs to <br>
  // For a protein atom the method returns the 3-letter abbreviation
  // of its residue. Returns 0 for other types of atoms.
  //static const char* const residueLongName(const char resChar);
  static std::string const residueLongName(const char resChar);
 //// the same for nucleic acid
  static std::string const nucleicAcidLongName(const char resChar);


  ////
  // Given a charachter string with a CIF residue type code, the
  // method returns a one-letter abbreviation of the type of the
  // residue that the atom belongs to: <br>
  // - For a nucleic acid atom returns A for adenine, G for guanine, T
  // for thymine, U for uracil and X for unknown type. <br>
  // - For a protein atom, returns a 3-letter abbreviation of the
  // residue name <br>
  //static char residueShortName(const char *const resName);

  ////
  //
  //static std::vector<unsigned short> connectedAtoms(const std::string& CIFline);
  */
private:
  // Private constructor to prevent object instantiation. Class was meant to
  // be use as a static method pool.
  CIF() {}

  // A CIF dictionary for ATOM records: key = field_name, value = column number
  static std::map<std::string, unsigned int> atomSiteTable_;

};

#endif
