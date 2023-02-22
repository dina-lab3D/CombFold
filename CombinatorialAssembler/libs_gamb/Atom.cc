#include "Atom.h"
#include "PDB.h"
#include "CIF.h"

#include <string>
#include <string.h>
#include <math.h>
#include "macros.h"

const unsigned short Atom::resNameLen            =   3;
const unsigned short Atom::atomNameLen           =   4;
const float Atom::defaultOccupancy               =   1.0;
const float Atom::defaultTempFactor              =   0.0;

std::string Atom::parseAtomName(const char* const atomName_) {
  std::string atomname;
  switch (strlen(atomName_)) {
  case 1:
    atomname = " " + std::string(atomName_) + "  ";
    break;
  case 2:
    atomname = " " + std::string(atomName_) + " ";
    break;
  case 3:
    atomname = " " + std::string(atomName_ );
    break;
  default:
    atomname = std::string(atomName_).substr(0,4);
  }

  return atomname;
}


Atom::Atom()
  : Vector3(), atomEntryType(ATOM), chId(' '), altLoc(' '), atomId(0),
    residueId(0), residueSeqICode(' '), atomName("    "), resType(' '), resName("   "),
    occupancy(defaultOccupancy), tempFactor(defaultTempFactor),sseInfo(std::pair<char,short>(0,0))
{}

Atom::Atom(const Vector3& p, const char ch, const unsigned int aid,
	   const unsigned int resid,
	   const char* const atomName_, const char resTp, const char resICode) :
  Vector3(p), atomEntryType(ATOM), chId(ch), altLoc(' '), atomId(aid),
  residueId(resid), residueSeqICode(resICode), resType(resTp),
  occupancy(defaultOccupancy), tempFactor(defaultTempFactor),sseInfo(std::pair<char,short>(0,0))
{
  atomName = parseAtomName(atomName_);
  resName = std::string(PDB::residueLongName(resTp)).substr(0,3);
}

Atom::Atom(const Vector3& p, const AtomEntryType atp, const char ch,
    const unsigned int aid, const unsigned int resid,
    const char* const atomName_, const char resTp, const char* const res, const char resICode) :
    Vector3(p), atomEntryType(atp), chId(ch), altLoc(' '), atomId(aid),
    residueId(resid), residueSeqICode(resICode), resType(resTp), resName("   "),
    occupancy(defaultOccupancy), tempFactor(defaultTempFactor),sseInfo(std::pair<char,short>(0,0))
{
  atomName = parseAtomName(atomName_);
  resName = std::string(res).substr(0,3);
}

Atom::Atom(const Vector3& p, const AtomEntryType atp, const char ch,
	   const unsigned int aid, const unsigned int resid,
	   const char* const atomName_, const char resTp, const char* const res,
	   const float occup, const float tmp, const char resICode) :
  Vector3(p), atomEntryType(atp), chId(ch), altLoc(' '), atomId(aid),
  residueId(resid), residueSeqICode(resICode), resType(resTp), resName("   "),
  occupancy(occup), tempFactor(tmp), sseInfo(std::pair<char,short>(0,0))
{
  atomName = parseAtomName(atomName_);
  resName = std::string(res).substr(0,3);
}

Atom::Atom(const Vector3& p, unsigned int aid) :
  Vector3(p), atomEntryType(ATOM), chId(' '), altLoc(' '), atomId(aid), residueId(0),
  residueSeqICode(' '), atomName("    "), resType(' '), resName("   "),
  occupancy(defaultOccupancy), tempFactor(defaultTempFactor), sseInfo(std::pair<char,short>(0,0))
{}

Atom::Atom(const std::string& line, bool cif) {
  if(cif) initCIF(line);
  else initPDB(line);
}

void Atom::initPDB(const std::string& PDBrec) {
  (*this)[0] = PDB::atomXCoord(PDBrec);
  (*this)[1] = PDB::atomYCoord(PDBrec);
  (*this)[2] = PDB::atomZCoord(PDBrec);
  atomEntryType = PDB::getAtomEntryType(PDBrec);
  chId = PDB::atomChainId(PDBrec);
  altLoc = PDB::atomAltLocIndicator(PDBrec);
  atomId = PDB::atomIndex(PDBrec);
  residueId = PDB::atomResidueIndex(PDBrec);
  residueSeqICode = PDB::atomResidueICode(PDBrec);
  resType = PDB::atomResidueType(PDBrec);
  occupancy = PDB::getOccupancy(PDBrec);
  tempFactor = PDB::getTempFactor(PDBrec);
  sseInfo = std::pair<char,short>(0,0);
  atomName = PDB::atomType(PDBrec);
  resName = PDB::getAtomResidueLongName(PDBrec);
}

void Atom::initCIF(const std::string& CIFrec) {
  (*this)[0] = CIF::atomXCoord(CIFrec);
  (*this)[1] = CIF::atomYCoord(CIFrec);
  (*this)[2] = CIF::atomZCoord(CIFrec);
  atomEntryType = CIF::getAtomEntryType(CIFrec);
  chId = CIF::atomChainIdAuthor(CIFrec)[0];
  //altLoc = CIF::atomAltLocIndicator(CIFrec); //TODO
  atomId = CIF::atomIndex(CIFrec);
  residueId = CIF::atomResidueIndex(CIFrec);
  //residueSeqICode = CIF::atomResidueICode(CIFrec);
  occupancy = CIF::getOccupancy(CIFrec);
  tempFactor = CIF::getTempFactor(CIFrec);
  sseInfo = std::pair<char,short>(0,0);
  atomName = parseAtomName(CIF::atomType(CIFrec).c_str());
  resName = CIF::getAtomResidueLongName(CIFrec);
  resType = PDB::residueShortName(resName);
}

char Atom::chainId() const {
  return chId;
}

char Atom::altLocIndicator() const {
  return altLoc;
}

unsigned int Atom::atomIndex() const {
  return atomId;
}

int Atom::residueIndex() const {
  return residueId;
}

std::string Atom::residueSequenceID() const {
  static const size_t MAX_CHARS_LENGTH = 5;
  char cResId[MAX_CHARS_LENGTH];
  sprintf(cResId, "%d", residueId);
  std::string sResId(cResId);
  sResId += residueSeqICode;
  return sResId;
}

char Atom::residueICode() const {
  return residueSeqICode;
}

const char *const Atom::type() const {
  return atomName.c_str();
}

char Atom::residueType() const {
  return resType;
}

const char *const Atom::residueName() const
{
  return resName.c_str();
}

void Atom::setChainId(const char newChainId) {
  chId = newChainId;
}

void Atom::setAtomId(unsigned int newAtomId) {
  atomId = newAtomId;
}

AtomEntryType Atom::getAtomEntryType() const {
  return atomEntryType;
}

float Atom::getOccupancy() const {
  return occupancy;
}

float Atom::getTempFactor() const {
  return tempFactor;
}

void Atom::setTempFactor(float ftemp) {
  tempFactor=ftemp;
}

bool Atom::isBackbone() const
{
  const char* atomType = type();
  if (atomType[1] == 'N' && atomType[2] == ' ')
    return true;
  else if (atomType[1] == 'C' && atomType[2] == 'A' && atomType[3] == ' ')
    return true;
  else if (atomType[1] == 'C' && atomType[2] == ' ')
    return true;
  else
    return atomType[1] == 'O' && atomType[2] == ' ';
}

bool Atom::isCA() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'A' && atomType[3] == ' ';
}

bool Atom::isCB() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'B' && atomType[3] == ' ';
}

bool Atom::isN() const
{
  const char* atomType = type();
  return atomType[1] == 'N' && atomType[2] == ' ';
}

bool Atom::isC() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == ' ';
}

bool Atom::isO() const
{
  const char* atomType = type();
  return atomType[1] == 'O' && atomType[2] == ' ';
}

bool Atom::isCG() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'G' && atomType[3] == ' ';
}

bool Atom::isCD() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'D' && atomType[3] == ' ';
}

bool Atom::isCE() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'E' && atomType[3] == ' ';
}

bool Atom::isCD1() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'D' && atomType[3] == '1';
}

bool Atom::isCD2() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'D' && atomType[3] == '2';
}

bool Atom::isCG1() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'G' && atomType[3] == '1';
}

bool Atom::isCG2() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'G' && atomType[3] == '2';
}

bool Atom::isCE1() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'E' && atomType[3] == '1';
}

bool Atom::isCE2() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'E' && atomType[3] == '2';
}

bool Atom::isCE3() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'E' && atomType[3] == '3';
}

bool Atom::isCZ() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'Z' && atomType[3] == ' ';
}

bool Atom::isCZ2() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'Z' && atomType[3] == '2';
}

bool Atom::isCZ3() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'Z' && atomType[3] == '3';
}

bool Atom::isCH2() const
{
  const char* atomType = type();
  return atomType[1] == 'C' && atomType[2] == 'H' && atomType[3] == '2';
}

bool Atom::isOD1() const
{
  const char* atomType = type();
  return atomType[1] == 'O' && atomType[2] == 'D' && atomType[3] == '1';
}

bool Atom::isOD2() const
{
  const char* atomType = type();
  return atomType[1] == 'O' && atomType[2] == 'D' && atomType[3] == '2';
}

bool Atom::isOE1() const
{
  const char* atomType = type();
  return atomType[1] == 'O' && atomType[2] == 'E' && atomType[3] == '1';
}

bool Atom::isOE2() const
{
  const char* atomType = type();
  return atomType[1] == 'O' && atomType[2] == 'E' && atomType[3] == '2';
}

bool Atom::isOH() const
{
  const char* atomType = type();
  return atomType[1] == 'O' && atomType[2] == 'H' && atomType[3] == ' ';
}

bool Atom::isOG() const
{
  const char* atomType = type();
  return atomType[1] == 'O' && atomType[2] == 'G' && atomType[3] == ' ';
}

bool Atom::isOG1() const
{
  const char* atomType = type();
  return atomType[1] == 'O' && atomType[2] == 'G' && atomType[3] == '1';
}

bool Atom::isNE() const
{
  const char* atomType = type();
  return atomType[1] == 'N' && atomType[2] == 'E' && atomType[3] == ' ';
}

bool Atom::isND1() const
{
  const char* atomType = type();
  return atomType[1] == 'N' && atomType[2] == 'D' && atomType[3] == '1';
}

bool Atom::isND2() const
{
  const char* atomType = type();
  return atomType[1] == 'N' && atomType[2] == 'D' && atomType[3] == '2';
}

bool Atom::isNE1() const
{
  const char* atomType = type();
  return atomType[1] == 'N' && atomType[2] == 'E' && atomType[3] == '1';
}

bool Atom::isNE2() const
{
  const char* atomType = type();
  return atomType[1] == 'N' && atomType[2] == 'E' && atomType[3] == '2';
}

bool Atom::isNH1() const
{
  const char* atomType = type();
  return atomType[1] == 'N' && atomType[2] == 'H' && atomType[3] == '1';
}

bool Atom::isNH2() const
{
  const char* atomType = type();
  return atomType[1] == 'N' && atomType[2] == 'H' && atomType[3] == '2';
}

bool Atom::isNZ() const
{
  const char* atomType = type();
  return atomType[1] == 'N' && atomType[2] == 'Z' && atomType[3] == ' ';
}

bool Atom::isSD() const
{
  const char* atomType = type();
  return atomType[1] == 'S' && atomType[2] == 'D' && atomType[3] == ' ';
}

bool Atom::isSG() const
{
  const char* atomType = type();
  return atomType[1] == 'S' && atomType[2] == 'G' && atomType[3] == ' ';
}

bool Atom::isH() const
{
  const char* atomType = type();
  return atomType[1] == 'H' || atomType[1] == 'D' || atomType[0] == 'H' || atomType[0] == 'D';
}

bool Atom::isWater() const
{
  const char* resName = residueName();
  return strncmp(resName, "HOH", 3) || strncmp(resName,"DOD", 3);
}

/* RNA backbone
//        OP1    H5' H4'  O4'   residue
//         |      |    \ /   \  /
//        -P-O5'-C5'---C4'    C1'
//         |      |     \     / \
//        OP2    H5''   C3'--C2' H1'
//                     / \   / \
//                  O3' H3' O2' H2''
//                   |       |
//                          H2'
*/
bool Atom::isRNABackbone() const {
  const char* atomType = type();
  return ((atomType[1] == 'O' && atomType[2] == 'P') || isSugarAtom() || isP());
}

bool Atom::isP() const {
  const char* atomType = type();
  return atomType[1] == 'P' && atomType[2] == ' ';
}

bool Atom::isSugarAtom() const {
  const char* atomType = type();
  // '*' for old convention
  return atomType[3] == '*' || atomType[3] == '\'';
}

bool Atom::isSugarCarbon() const
{
  const char *atomType = type();
  return atomType[1] == 'C' && atomType[3] == '\'';
}

bool Atom::isSugarOxygen() const
{
  const char *atomType = type();
  return atomType[1] == 'O' && atomType[3] == '\'';
}

bool Atom::isPhosphateOxygen() const
{
  const char *atomType = type();
  return atomType[1] == 'O' && atomType[2] == 'P';
}

bool Atom::isNitrogenousBaseAtom() const
{
  const char *atomType = type();
  return (atomType[1] == 'N' || atomType[1] == 'O' || atomType[1] == 'C')
	  && atomType[2] >= '0' && atomType[2] <= '9'
	  && (atomType[3] == ' ')
	  ;
}


std::ostream& operator<<(std::ostream& s, const Atom& at)
{
  // first 6 chars (1-6) stand for field type
  s.setf(std::ios::left, std::ios::adjustfield);
  s.width(6);
  if (at.atomEntryType == ATOM)
    s << "ATOM";
  else if (at.atomEntryType == HETATM)
    s << "HETATM";
  else
    s << "UNK";

  // next 5 chars (7-11) stand for atom id
  s.setf(std::ios::right, std::ios::adjustfield);
  s.width(5);
  s << at.atomId;

  // skip an undefined position (12)
  s.width(1);
  s << " ";

  // 4 positions (13-16) for atom name
  s.width(4); //worked well with gcc-2.95
  s << at.atomName;

  // skip the alternate indication position (17)
  s.width(1);
  s << at.altLoc;

  // 3 positions for residue name (18-20)
  //s.width(3); worked well with gcc-2.95
  if(at.resName.length()==3) {
    s << at.resName;
  } else {
    s << "   ";
  }

  // 1 position for chain identifier (22) + 1 skipped (21)
  s.width(1);
  s << " " << at.chId; // skip the chain identifier

  // 4 position for residue number (23-26)
  s.width(4);
  s << at.residueIndex(); // residue sequence number
  s.width(1);
  s << at.residueICode(); // residue insertion code (27)
  s.setf(std::ios::fixed, std::ios::floatfield);
  s << "   "; // skip 3 undefined positions (28-30)

  // coordinates (31-38,39-46,47-54)
  s.width(8);
  s.precision(3);
  s << at[0];
  s.width(8);
  s.precision(3);
  s << at[1];
  s.width(8);
  s.precision(3);
  s << at[2];
  // occupancy (55-60)
  s.width(6);
  s.precision(2);
  s << at.occupancy;
  // temp. factor (61-66)
  s.width(6);
  s.precision(2);
  s << at.tempFactor;

  //  s << endl;
  return s;
}

void Atom::output2cif(std::ostream& s, int modelNum) const
{
  if (atomEntryType == ATOM)
    s << "ATOM ";
  else if (atomEntryType == HETATM)
    s << "HETATM ";

  s << atomId; //_atom_site.id
  s << " ? "; // we don't save the element (_atom_site.type_symbol)
  s << atomName; // _atom_site.label_atom_id
  s << " . "; // _atom_site.label_alt_id
  s << resName;//_atom_site.label_comp_id
  s << " ";
  s << chId; //_atom_site.label_asym_id
  s << " ? "; //_atom_site.label_entity_id  - chain number
  s << residueId; //_atom_site.label_seq_id
  s << " ? "; //_atom_site.pdbx_PDB_ins_code

  s << (*this)[0] << " " << (*this)[1] << " " << (*this)[2] << " "; // coordinates
  s << occupancy << " " << tempFactor << " ? ";
  s << residueId; //_atom_site.auth_seq_id
  s << " ";
  s << resName;//_atom_site.auth_comp_id
  s << " ";
  s << chId; //_atom_site.auth_asym_id
  s << " ";
  s << atomName; // _atom_site.auth_atom_id
  s << " ";
  s << modelNum; // _atom_site.pdbx_PDB_model_num
}
