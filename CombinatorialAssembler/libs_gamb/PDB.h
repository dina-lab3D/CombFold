#ifndef _PDB_h
#define _PDB_h

#include <string>
#include <iostream>
#include <vector>

////
enum  AtomEntryType {UNK = 0, ATOM = 1, HETATM = 2};

/*
CLASS
  PDB

  The PDB class defines a set of tools used for reading the Protein Data Bank
  file format, used to describe molecular structures, like proteins,
  drugs and Nucleic Acids
KEYWORDS
  PDB, protein, Nucleic Acids (RNA, DNA) , record, atom, residue,
  structure, coordinate, field

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG

GOALS
  Class PDB is designed to encapsulate a set of methods used for reading the
  PDB file format. Methods are static requiring no object instantiation
  before their use.

USAGE
  A PDB file is comprised of 80-char lines called records. These
  records describe different properties of a molecular structure
  (mainly proteins). A record can describe an atom coordinate, an
  amino-acid sequence, an alpha helix or some other chemical or
  geometrical properties of the molecule, depending on the record's
  type. For more information see the PDB_site.

  Class PDB's methods are all static. Hence the class's constructor was
  declared private so that no objects of type PDB will be declared. PDB's
  methods recieve a PDB record line as an argument and return requested data
  found in the record.
  EXAMPLE
    char rec[81];    pdbFile.getline(line, 81);         // read line from stream.
    if (PDB::isATOMrec(rec))             // if record is of type ATOM
      Vector3 v(PDB::atomXCoord(rec),    // instantiate a vector with
                PDB::atomYCoord(rec),    // atom's coordinates.
                PDB::atomZCoord(rec));
  END
  Note that one must make sure the record is of the correct type before asking
  for data found in the record. For example you cannot ask for an x-coordinate
  of an atom before making sure that you are dealing with an ATOM record.
  Method names are prefixed by the type of record they are meant to handle
  (e.g. atomXCoord is meant to extract the X coordinate from an ATOM record).
  Methods beginning with the is prefix (e.g. isATOMrec) check for record type.
*/
class PDB {

public:

  // GROUP: PDB record type determination methods.

  ////
  // Returns true if the given record is an ATOM record.
  static bool isATOMrec(const std::string& pdb_line);

  ////
  // Returns true if the given record is a HETATM record.
  static bool isHETATMrec(const std::string& pdb_line);

  ////
  // Returns true if the given record is an ATOM record of a C-alpha atom.
  static bool isCaATOMrec(const std::string& pdb_line);

  ////
  // Returns true if the given record is a record of a nucleic acid residue
  // <br> Author: Oranit Dror (oranit@tau.ac.il)
  static bool isNucleicRecord(const std::string& pdb_line);

  //// returns true if the residue name is of nucleic acid
  static bool isNucleicAcid(const std::string& residueName);

  ////
  // Returns true if the given record is a record of a connection between atoms
  // <br> Author: Bluvshtein Nelly (bluvsh@post.tau.ac.il)
  static bool isCONECTrec(const std::string& pdb_line);

  static bool isMODELrec(const std::string& pdb_line);

  static bool isENDMDLrec(const std::string& pdb_line);


  // GROUP: ATOM record costants and methods.
  // ATOM Record Format (from PDB)
  // COLUMNS        DATA TYPE       CONTENTS
  //  --------------------------------------------------------------------------------
  // 1 -  6         Record name     "ATOM  "
  // 7 - 11         Integer         Atom serial number.
  // 13 - 16        Atom            Atom name.
  // 17             Character       Alternate location indicator.
  // 18 - 20        Residue name    Residue name.
  // 22             Character       Chain identifier.
  // 23 - 26        Integer         Residue sequence number.
  // 27             AChar           Code for insertion of residues.
  // 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
  // 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
  // 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
  // 55 - 60        Real(6.2)       Occupancy.
  // 61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
  // 73 - 76        LString(4)      Segment identifier, left-justified.
  // 77 - 78        LString(2)      Element symbol, right-justified.
  // 79 - 80        LString(2)      Charge on the atom.

  static const unsigned short atomEntryType     =   0;
  static const unsigned short atomIndexField    =   6;
  static const unsigned short atomTypeField     =  12;
  static const unsigned short atomAltLocField   =  16;
  static const unsigned short atomResTypeField  =  17;
  static const unsigned short atomChainIdField  =  21;
  static const unsigned short atomResIndexField =  22;
  static const unsigned short atomResIndexCharField =  26;
  static const unsigned short atomXCoordField   =  30;
  static const unsigned short atomYCoordField   =  38;
  static const unsigned short atomZCoordField   =  46;
  static const unsigned short atomOccupancy     =  54;
  static const unsigned short atomTempFactor    =  60;

  //// Returns an ATOM record's X coordinate.
  static float atomXCoord(const std::string& pdb_line);

  //// Returns an ATOM record's Y coordinate.
  static float atomYCoord(const std::string& pdb_line);

  //// Returns an ATOM record's Z coordinate.
  static float atomZCoord(const std::string& pdb_line);

  //// Returns an ATOM record's atom index number.
  static int atomIndex(const std::string& pdb_line);

  //// Returns a pointer to an ATOM's record atom type string.
  // The atom type is a 4 charachter long field. A pointer is returned to
  // the correct position withing the given string PDBrec. Hence the pointer
  // should not be freed or deleted.
  static std::string atomType(const std::string& pdb_line);

  //// Returns an ATOM record's alternative location indicator.
  static char  atomAltLocIndicator(const std::string& pdb_line);

  //// Returns an ATOM record's chain id.
  static char  atomChainId(const std::string& pdb_line);

  //// Returns an ATOM record's residue index.
  static short atomResidueIndex(const std::string& pdb_line);

  ////
  // code for insertion of residues
  // Returns the value stored in the atomResIndexCharField of the given PDB
  // record (column 26).
  //<br>
  // ATOM    279  CA  SER    60       8.961   9.167  37.481  1.00 20.10      1SGT 477 <br>
  // ATOM    285  CA  GLY    60A      5.640  10.901  38.407  1.00 21.41      1SGT 483
  // <br> In this case for atom no. 285 the atomResidueICode method
  // will return 'A' and for the atom no. 279 it will return ' '
  static char atomResidueICode(const std::string& pdb_line);

  ////
  // Gets a PDB record of an atom and returns a one-letter
  // abbreviation of the type of the residue that the atom belongs
  // to. <br>
  // For a nucleic acid atom returns A for adenine, G for guanine, T
  // for thymine, U for uracil and X for unknown type.
  static char atomResidueType(const std::string& pdb_line);

  ////
  // Gets a one-letter abbreviation of the type of the residue that
  // the atom belongs to <br>
  // For a protein atom the method returns the 3-letter abbreviation
  // of its residue. Returns 0 for other types of atoms.
  //static const char* const residueLongName(const char resChar);
  static std::string const residueLongName(const char resChar);

  //// the same for nucleic acid
  static std::string const nucleicAcidLongName(const char resChar);

  //// Returns the full string of the Residue Name field
  static std::string const getAtomResidueLongName(const std::string& pdb_line);

  //// Return the type of the PDB entry, currently only ATOM and HETATM are identified.
  static AtomEntryType getAtomEntryType(const std::string& pdb_line);

  //// Returns the occupancy,  which is a measure of the fraction of molecules in the
  //   crystal in which the current atom actually occupies the specified position
  static float getOccupancy(const std::string& pdb_line);

  //// Returns the temperature factor,  which is a measure of how much an atom
  // oscillates or vibrates around the specified position
  static float getTempFactor(const std::string& pdb_line);

  ////
  // Given a charachter string with a PDB residue type code, the
  // method returns a one-letter abbreviation of the type of the
  // residue that the atom belongs to: <br>
  // - For a nucleic acid atom returns A for adenine, G for guanine, T
  // for thymine, U for uracil and X for unknown type. <br>
  // - For a protein atom, returns a 3-letter abbreviation of the
  // residue name <br>
  static char residueShortName(const std::string& resName);

  ////
  //
  static std::vector<unsigned short> connectedAtoms(const std::string& pdb_line);


  //// Selector is a general purpose class used to select records from a PDB
  // file. Using descendants of this class one may implement arbitrary
  // selection functions with operator() and pass them to PDB reading functions
  // for object selection. (See Molecule for examples).
  class Selector {
  public:
    virtual bool operator()(const char*) const {
      return true;
    }
    virtual ~Selector() {}
  };

  //// Defines a selector that will pick only C-alpha atoms.
  class CAlphaSelector : public Selector{
   public:
    bool operator()(const char *PDBrec)const{
      std::string type = atomType(PDBrec);
      char alt = atomAltLocIndicator(PDBrec);
      if(type[1] == 'C' && type[2] == 'A' && type[3] == ' ' &&
         (alt == ' ' || alt == 'A')) {
	// && (*(type+4) == ' ' || *(type+4) == 'A')){
	return true;
      }else
	return false;
    }
  };

  //// Defines a selector that will pick only C-beta atoms.
  class CBetaSelector: public Selector{
   public:
    bool operator()(const char *PDBrec) const {
      const std::string residueName = getAtomResidueLongName(PDBrec);
      if (residueName[0] == 'G' && residueName[1] == 'L' && residueName[2] == 'Y') {
	return CAlphaSelector()(PDBrec);
      }

      std::string type = atomType(PDBrec);
      char alt = atomAltLocIndicator(PDBrec);
      if(type[1] == 'C' && type[2] == 'B' && type[3] == ' ' &&
         (alt == ' ' || alt == 'A'))
	return true;
      else
	return false;
    }
  };

  //// Defines a selector that will pick only C atoms. (not Ca or Cb)
  class CSelector: public Selector {
     public:
    bool operator()(const char *PDBrec) const {
      std::string type = atomType(PDBrec);
      return (type[1] == 'C' && type[2] == ' ');
    }
  };

  //// Defines a selector that will pick only N atoms.
  class NSelector: public Selector {
   public:
    bool operator()(const char *PDBrec) const {
      std::string type = atomType(PDBrec);
      return (type[1] == 'N' && type[2] == ' ');
    }
  };

  //// Defines a selector that picks backbone atoms
  class BBSelector : public Selector {
   public:
    bool operator()(const char *PDBrec) const {
      std::string type = atomType(PDBrec);
      if (type[1] == 'N' && type[2] == ' ')
        return true;
      else if (type[1] == 'C' && type[2] == 'A' && type[3] == ' ')
        return true;
      else if (type[1] == 'C' && type[2] == ' ')
        return true;
      else if (type[1] == 'O' && type[2] == ' ')
        return true;
      return false;
    }
  };

  //// Defines a selector that will pick every atom.
  class AllSelector : public Selector {
   public:
    bool operator()(const char* ) const{
      return true;
    }
  };

  //// Selector that picks atoms of a given chain.
  class ChainSelector : public Selector {
   public:
    // constructor
    ChainSelector(const std::string &chainID): chains(chainID) {}
    virtual ~ChainSelector() {}
    bool operator()(const char* PDBrec) const {
      for(int i=0; i< (int)chains.length(); i++)
	if(atomChainId(PDBrec) == chains[i])
	  return true;
      return false;
    }
   private:
    std::string chains;
  };


  //// Selector that check if the line is water record
  class WaterSelector : public Selector {
   public:
    bool operator()(const char* PDBrec) const {
      return ((PDBrec[17]=='H' && PDBrec[18]=='O' && PDBrec[19]=='H') ||
	      (PDBrec[17]=='D' && PDBrec[18]=='O' && PDBrec[19]=='D'));
    }
  };


  //// Selector that check if the line is hydrogen record
  class HydrogenSelector : public Selector {
   public:
    bool operator()(const char* PDBrec) const {
      return (PDBrec[13] == 'H' || PDBrec[13] == 'D');
    }
  };

  //// Selector that picks non water and non hydrogen atoms
  class WaterHydrogenUnSelector : public Selector {
   public:
    bool operator()(const char* PDBrec) const {
      return( !((PDBrec[17]=='H' && PDBrec[18]=='O' && PDBrec[19]=='H') ||
		(PDBrec[17]=='D' && PDBrec[18]=='O' && PDBrec[19]=='D')) &&
	      ! (PDBrec[13] == 'H' || PDBrec[13] == 'D' || PDBrec[12] == 'H' || PDBrec[12] == 'D') &&
              (PDB::atomAltLocIndicator(PDBrec) == ' ' || PDB::atomAltLocIndicator(PDBrec) == 'A') ); // skip alternatives
    }
  };

  //// Selector that picks non water atoms
  class WaterUnSelector : public Selector {
   public:
    bool operator()(const char* PDBrec) const {
      return( !((PDBrec[17]=='H' && PDBrec[18]=='O' && PDBrec[19]=='H') ||
		(PDBrec[17]=='D' && PDBrec[18]=='O' && PDBrec[19]=='D')));
	}
  };

  ////
  // A PDB Selector that picks only Phosphate atoms. <br>
  class PSelector : public Selector {
  public:
    bool operator()(const char *PDBrec) const {
      std::string type = PDB::atomType(PDBrec);
      return ((type[1] == 'P' && type[2] == ' ') || (type[0] == 'P' && type[1] == ' '));
    }
  };

  ////
  // A PDB Selector that picks only nucleic acid atoms. <br>
  class NucleicSelector : public Selector {
  public:
    bool operator()(const char *PDBrec) const {
      return PDB::isNucleicRecord(PDBrec);
    }
  };

  ////
  // A PDB Selector that ignores all alternative location atoms. <br>
  class IgnoreAlternativesSelector : public Selector {
  public:
    bool operator()(const char *PDBrec) const {
      return ((PDB::atomAltLocIndicator(PDBrec) == ' ') ||
	      (PDB::atomAltLocIndicator(PDBrec) == 'A'));
    }
  };

private:
  // Private constructor to prevent object instantiation. Class was meant to
  // be use as a static method pool.
  PDB() {}


};

#endif
