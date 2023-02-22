#ifndef _Atom_h
#define _Atom_h

#include "Vector3.h"
#include "RigidTrans3.h"
#include "PDB.h"
#include "AminoAcid.h"
#include <string>

/*
CLASS
  Atom

  This is a general class to hold an atom as it appears in a PDB file.

KEYWORD
    atom, PDB, particle

AUTHORS
    Meir Fuchs (meirfux@math.tau.ac.il), Zipi Fligelman (zipo@math.tau.ac.il)

    copyright: SAMBA Group , Tel-Aviv Univ. Israel, 1997, 1999.

CHANGES LOG
<UL>

<LI>
13.11.07 Elad Donsky
Adding alternate location field
<\LI>
<LI>
2.08.07 Efrat
Adding isCB method
<\LI>
<LI>
10.05.05 Nelly Bluvshtein
1.
String data member for residueSeqID ( = "45A") was replaced by residueSeqICode ( = 'A').
The residue number is stored in residueId (as it was).
Atom constructors were changed to get resICode (with default = ' ').
residueSequenceID() returns concatenation of residueId and residueSeqICode.
2.
isH(),isWater() added for Molecule<Atom>::WaterHydrogenUnSelector
<\LI>
<LI>
10.10.04 Dina Schneidman
Residue name is stored in a 3-characters string instead of char[3]
<\LI>
<LI>
16.07.04 Oranit Dror
A small update in the << operator
<\LI>
<LI>
17.3.04
pair(char,short) sseInfo added.
<LI>
13.01.04 Dina and Oranit
Atom name is stored in a 4-characters string instead of char[4]. If
the name length is shorter than 4 characters then blanks are appended.
The "type" method is now returns a terminated char[] ('\0' at the end).
<LI>
15.01.03 Max
<p>
class Atom : public Particle -> class Atom : public Vector3
<p>
All functions are changed accordingly.
<\LI>
<LI>
Occupancy and temperature factor were added. 8.01.03 Alex
<\LI>
<LI>
getAtomEntryType() - function that tells if an atom is
ATOM or HETATM record in PDB file.  8.01.03 Alex
<\LI>
<LI>
isBackbone() - function that returns true for backbone atoms
 8.01.03 Dina
<\LI>
<LI>
Added: operator RESIDUE_INDEX()const (14.12.03 Max)
<\LI>
</UL>

USAGE
*/
class Atom : public Vector3
{

public:
  // GROUP: Constructors

  //// An empty constructor that gives the atom members null values
  Atom();

  //// An explicit constructor
  Atom(const Vector3& p, const char ch, const unsigned int aid,
       const unsigned int resid,
       const char* const t,
       const char resTp,
       const char resICode = ' ');

  ////
  Atom(const Vector3& p, const AtomEntryType atp, const char ch,
       const unsigned int aid, const unsigned int resid,
       const char* const t, const char resTp, const char* const res, const char resICode = ' ');

  ////
  Atom(const Vector3& p, const AtomEntryType atp, const char ch,
       const unsigned int aid, const unsigned int resid,
       const char* const atomName_, const char resTp, const char* const res,
       const float occup, const float tmp, const char resICode = ' ');

  ////
  Atom(const Vector3& p, unsigned int aid);

  //// A constructor that intiate an atom type from a line of a PDB or CIF record
  explicit Atom(const std::string& line, bool cif = false);

  static const unsigned short resNameLen;
  static const unsigned short atomNameLen;
  static const float defaultOccupancy;
  static const float defaultTempFactor;

  // GROUP: Inspectors

  //// Returns the chainId - important for working with complexes
  char chainId() const;

  //// Returns the alternative location indicator.
  char altLocIndicator() const;

  //// returns the atom index in the PDB file
  unsigned int atomIndex() const;

  //// returns the residue (thas the atom belongs to) index in the PDB file
  int residueIndex() const;

  // Return the value stored in the ResSeq field of the given PDB
  // record (columns 23-26).
  // The returned value is the same as the value returned by the
  // residueIndex method, except for cases like:
  // ATOM    279  CA  SER    60       8.961   9.167  37.481  1.00 20.10      1SGT 477
  // ATOM    285  CA  GLY    60A      5.640  10.901  38.407  1.00 21.41      1SGT 483
  std::string residueSequenceID() const;

  //// returns space if there are no alternatives for this atom, column 26 in PDB format
  char residueICode() const;

  //// Returns the atom type
  const char *const type() const;

  //// Returns the residue type in one letter code.
  char residueType() const;

  //// Returns the residue type in 3-letter mnemonic. If no matching mnemonic
  // is found a NULL pointer is returned.
  const char *const residueName() const;

  //// Return the type of the PDB entry, currently only ATOM and HETATM are identified.
  AtomEntryType getAtomEntryType() const;

  //// Returns the occupancy,  which is a measure of the fraction of molecules in the
  //   crystal in which the current atom actually occupies the specified position
  float getOccupancy() const;

  //// Returns the temperature factor,  which is a measure of how much an atom
  // oscillates or vibrates around the specified position
  float getTempFactor() const;

  //// Sets new temperature factor
  void setTempFactor(float ftemp);

  //// Returns true if the atom is backbone atom: N, Ca, C, O
  bool isBackbone() const;

  //// Returns true if the atom is CA atom
  bool isCA() const;

  //// Returns true if the atom is CB atom
  bool isCB() const;

  //// Returns true if the atom is N atom
  bool isN() const;

  //// Returns true if the atom is C atom
  bool isC() const;

  //// Returns true if the atom is O atom
  bool isO() const;

  //// Returns true if the atom is CG atom
  bool isCG() const;

  //// Returns true if the atom is CD atom
  bool isCD() const;

  //// Returns true if the atom is CE atom
  bool isCE() const;

  //// Returns true if the atom is CZ atom
  bool isCZ() const;

  //// Returns true if the atom is CD1 atom
  bool isCD1() const;

  //// Returns true if the atom is CD2 atom
  bool isCD2() const;

  //// Returns true if the atom is CG1 atom
  bool isCG1() const;

  //// Returns true if the atom is CG2 atom
  bool isCG2() const;

  //// Returns true if the atom is CE1 atom
  bool isCE1() const;

  //// Returns true if the atom is CE2 atom
  bool isCE2() const;

  //// Returns true if the atom is CE3 atom
  bool isCE3() const;

  //// Returns true if the atom is CZ2 atom
  bool isCZ2() const;

  //// Returns true if the atom is CZ3 atom
  bool isCZ3() const;

  //// Returns true if the atom is CH2 atom
  bool isCH2() const;


  //// Returns true if the atom is OD1 atom
  bool isOD1() const;

  //// Returns true if the atom is OD2 atom
  bool isOD2() const;

  //// Returns true if the atom is OE1 atom
  bool isOE1() const;

  //// Returns true if the atom is OE2 atom
  bool isOE2() const;

  //// Returns true if the atom is OH atom
  bool isOH() const;

  //// Returns true if the atom is OG atom
  bool isOG() const;

  //// Returns true if the atom is OG1 atom
  bool isOG1() const;


  //// Returns true if the atom is NE atom
  bool isNE() const;

  //// Returns true if the atom is ND1 atom
  bool isND1() const;

  //// Returns true if the atom is ND2 atom
  bool isND2() const;

  //// Returns true if the atom is ND1 atom
  bool isNE1() const;

  //// Returns true if the atom is ND2 atom
  bool isNE2() const;

  //// Returns true if the atom is ND1 atom
  bool isNH1() const;

  //// Returns true if the atom is ND2 atom
  bool isNH2() const;

  //// Returns true if the atom is NZ atom
  bool isNZ() const;

  //// Returns true if the atom is SD atom
  bool isSD() const;

  //// Returns true if the atom is SG atom
  bool isSG() const;

  //// Returns true if the atom is H atom
  bool isH() const;

  //// Returns true if the atom is water's atom
  bool isWater() const;

  //// Returns true if the atom is phosphate
  bool isP() const;

  //// Returns true if the atom is a backbone atom
  bool isRNABackbone() const;

  //// Returns true if the atom belongs to the sugar group of a nucleotide
  bool isSugarAtom() const;

  //// Returns true if the atom is a sugar Carbon atom
  bool isSugarCarbon() const;

  //// Returns true if the atom is a sugar Oxygen atom
  bool isSugarOxygen() const;

  //// Returns true if the atom is a phosphate Oxygen atom
  bool isPhosphateOxygen() const;

  //// Returns true if the atom belongs to a nitrogenous base
  bool isNitrogenousBaseAtom() const;

  // GROUP: Adaptors

  //// Sets the chain ID for various implementation like the VMD coloring
  // scheme, or RasMol web assignment. Please make sure you use it with
  // accordance to the appropriate software util.
  void setChainId(const char newChainId);

  //// Set atomId
  void setAtomId(unsigned int newAtomId);

  //// set residue index
  void setResidueIndex(int residueIndex) { residueId = residueIndex; }

  ////applies transformation on Atom coordinates
  Atom& operator*=(const RigidTrans3& rt) {
    (Vector3 &)(*this) *= rt;
    return *this;
  }

  ////Converts Atom into RESIDUE_INDEX
  //(atoms of the same residue type will have the same RESIDUE_INDEX)
  operator RESIDUE_INDEX() const {
    return AminoAcid<>::getResidueTypeIndx(resType);
  }

  ////Updates SSE information
  void setSSE(const std::pair<char,short>& isseInfo){sseInfo=isseInfo;}

  ////Updates SSE information
  std::pair<char,short> getSSE() const{ return sseInfo; }

  ////
  static std::string parseAtomName(const char* const atomName_);

  //// Writing back the PBD line - without the endline
  friend std::ostream& operator<<(std::ostream& s, const Atom& at);

  ////
  void output2cif(std::ostream& s, int modelNum=1) const;

protected:
  void initPDB(const std::string& PDBrec);
  void initCIF(const std::string& PDBrec);

protected:

  AtomEntryType  atomEntryType; // The type of the entity (ATOM/HETATM)
  char  chId;                 // chain.
  char  altLoc = ' ';          // alternative location indicator.
  unsigned int atomId;        // atom index.
  int residueId;            // residue index.
  char residueSeqICode = ' ';       // code for insertion of residues (for example for residue "45A"
                              // residueId = 45 and residueSeqICode = 'A')
  std::string atomName;       // Four-characters string. If the name
			      // is less than four characters, spaces
			      // are appended (this is the standard
			      // length according to the PDB).
  char resType;               // One letter Amino Acid Code
  std::string resName;        // Three-characters string. The full name of the entity:
                              // Residue name for ATOM record or ligand name for HETATM
  float occupancy;             // The occupancy n(j) of atom j is a measure of the fraction
                              // of molecules in the crystal in which atom j actually
                              // occupies the position specified in the model
  float tempFactor;            // The temperature factor or B-factor can be thought of as a
                              // measure of how much an atom oscillates or vibrates around
                              // the position specified in the model
  std::pair<char,short> sseInfo;     // secondary structure information: char sse_type - according to SecondStruct class,
                              // int sse_id - index of SSE in molecule

};

#endif
