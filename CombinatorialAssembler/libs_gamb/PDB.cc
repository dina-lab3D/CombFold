#include <stdlib.h>
#include <stdio.h>
#include <cctype>

#include "PDB.h"
#include <string.h>

static const unsigned short NUM_OF_PROTEIN_RESIDUES = 23;
static const std::string PROTEIN_RESIDUES =
"ALAA ARGR ASPD ASNN CYSC GLNQ GLUE GLYG HISH ILEI LEUL LYSK METM PHEF PROP SERS THRT TYRY TRPW VALV AS B GL Z UNKX";

static const unsigned short NUM_OF_NUCLEIC_ACIDS = 10;
static const std::string NUCLEIC_ACIDS = "  AA   GG   TT   UU   CC ADEA CYTC GUAG THYT URAU  DAA  DCC  DGG  DTT  DUU UNKX";

static const unsigned short NUM_OF_NUCLEIC_RESIDUES = 27;
static const char* NUCLEIC_RESIDUES[NUM_OF_NUCLEIC_RESIDUES] = {"  A", " +A", "  C", " +C", "  G", " +G", "  I", " +I", "  T", " +T", "  U", " +U", "1MA", "5MC", "OMC", "1MG", "2MG", "M2G", "7MG", "OMG", " YG", "H2U", "5MU", "PSU", "4SU", "MIA", "UNK"};

bool PDB::isATOMrec(const std::string& PDBline)
{
  return(PDBline[0]=='A' && PDBline[1]=='T' &&
	 PDBline[2]=='O' && PDBline[3]=='M');
}

bool PDB::isHETATMrec(const std::string& PDBline)
{
  return(PDBline[0]=='H' && PDBline[1]=='E' && PDBline[2]=='T' &&
	 PDBline[3]=='A' && PDBline[4]=='T' && PDBline[5]=='M');
}

bool PDB::isCaATOMrec(const std::string& PDBline)
{
  return(isATOMrec(PDBline) && PDBline[13] == 'C' && PDBline[14] == 'A');
}

bool PDB::isNucleicRecord(const std::string& PDBline) {
  const std::string residueName = getAtomResidueLongName(PDBline);
  return isNucleicAcid(residueName);
}

bool PDB::isNucleicAcid(const std::string& residueName)
{
  if (residueName[0] == 'U' && residueName[1] == 'N' && residueName[2] == 'K') {
    std::cerr << "PDB::isNucleicRecord - The given residue name is unknown :"
              << residueName<< ":" << std::endl;
    return false;
  }

  for (unsigned int i = 0 ; i < NUM_OF_NUCLEIC_RESIDUES ; i++) {
    if (residueName[0] == NUCLEIC_RESIDUES[i][0] && residueName[1] == NUCLEIC_RESIDUES[i][1] &&
	residueName[2] == NUCLEIC_RESIDUES[i][2]) {
      return true;
    }
  }

  return false;
}

bool PDB::isCONECTrec(const std::string& PDBline)
{
  return(PDBline[0]=='C' && PDBline[1]=='O' &&
	 PDBline[2]=='N' && PDBline[3]=='E' && PDBline[4]=='C' && PDBline[5]=='T');
}

bool PDB::isMODELrec(const std::string& PDBline)
{
  return (PDBline[0] == 'M' && PDBline[1] == 'O' && PDBline[2] == 'D' &&
          PDBline[3] == 'E' && PDBline[4] == 'L');
}

bool PDB::isENDMDLrec(const std::string& PDBline)
{
  return (PDBline[0] == 'E' && PDBline[1] == 'N' && PDBline[2] == 'D' &&
          PDBline[3] == 'M' && PDBline[4] == 'D' && PDBline[5] == 'L');
}

char PDB::residueShortName(const std::string& resName)
{
  for (unsigned int i = 0 ; i < NUM_OF_PROTEIN_RESIDUES*5; i += 5) {
    if (resName[0] == PROTEIN_RESIDUES[i] && resName[1] == PROTEIN_RESIDUES[i+1] &&
	resName[2] == PROTEIN_RESIDUES[i+2]) {
      return PROTEIN_RESIDUES[i+3];
    }
  }

  for (unsigned int i = 0 ; i < NUM_OF_NUCLEIC_RESIDUES ; i++) {
	if (resName[0] == NUCLEIC_RESIDUES[i][0] &&
		resName[1] == NUCLEIC_RESIDUES[i][1] &&
		resName[2] == NUCLEIC_RESIDUES[i][2]) {
	  return NUCLEIC_RESIDUES[i][2];
	}
  }
  // for (unsigned int i = 0 ; i < NUM_OF_NUCLEIC_ACIDS*5 ; i += 5) {
    // if (resName[0] == NUCLEIC_ACIDS[i] &&
        // resName[1] == NUCLEIC_ACIDS[i+1] &&
	// resName[2] == NUCLEIC_ACIDS[i+2]) {
      // return NUCLEIC_ACIDS[i+3];
    // }
  // }

  return 'X';
}

std::string const PDB::residueLongName(const char resChar)
{
  std::string retName;
  for (int i = 3 ; i < NUM_OF_PROTEIN_RESIDUES*5 ; i += 5) {
    if (resChar == PROTEIN_RESIDUES[i]) {
      retName = PROTEIN_RESIDUES.substr(i-3,3);
      return retName;
    }
  }

  std::cerr << "PDB::residueLongName - A unique long name does not exist for non-protein residue :"
            << resChar << ":" << std::endl;
  return retName;
}

std::string const PDB::nucleicAcidLongName(const char resChar)
{
  std::string retName;
  for (int i = 3 ; i < NUM_OF_NUCLEIC_ACIDS*5 ; i += 5) {
    if (resChar == NUCLEIC_ACIDS[i]) {
      retName = NUCLEIC_ACIDS.substr(i-3,3);
      return retName;
    }
  }

  std::cerr << "PDB::nucleicAcidLongName - A unique long name does not exist for nucleic acid:"
            << resChar << ":" << std::endl;
  return retName;
}

float PDB::atomXCoord(const std::string& PDBline)
{
  return (float)atof(PDBline.substr(atomXCoordField, 8).c_str());
}

float PDB::atomYCoord(const std::string& PDBline)
{
  return (float)atof(PDBline.substr(atomYCoordField, 8).c_str());
}

float PDB::atomZCoord(const std::string& PDBline)
{
  return (float)atof(PDBline.substr(atomZCoordField, 8).c_str());
}

int PDB::atomIndex(const std::string& PDBline)
{
  return atoi(PDBline.substr(atomIndexField).c_str());
}

char PDB::atomAltLocIndicator(const std::string& PDBline)
{
  return PDBline[atomAltLocField];
}

char  PDB::atomChainId(const std::string& PDBline)
{
  return PDBline[atomChainIdField];
}

short PDB::atomResidueIndex(const std::string& PDBline)
{
  return (short)atoi(PDBline.substr(atomResIndexField, 4).c_str());
}

char PDB::atomResidueICode(const std::string& PDBline)
{
  return PDBline[atomResIndexCharField];
}

char PDB::atomResidueType(const std::string& PDBline)
{
  return residueShortName(PDBline.substr(atomResTypeField, 3).c_str());
}

std::string PDB::atomType(const std::string& PDBline)
{
  return PDBline.substr(atomTypeField, 4);
}

const std::string PDB::getAtomResidueLongName(const std::string& PDBline)
{
  return PDBline.substr(atomResTypeField, 3);
}

AtomEntryType PDB::getAtomEntryType(const std::string& PDBline)
{
  int field_len = atomIndexField - atomEntryType;
  std::string entryType = PDBline.substr(atomEntryType, field_len);
  if (strcmp (entryType.c_str(), "HETATM") == 0)
    return(HETATM);
  else if (strncmp (entryType.c_str(), "ATOM", field_len - 2) == 0)
    return(ATOM);
  else
    return(UNK);
}

float PDB::getOccupancy(const std::string& PDBline)
{
  if(PDBline.length() >= atomOccupancy + 6)
    return (float)atof(PDBline.substr(atomOccupancy, 6).c_str());
  return 0.0;
}

float PDB::getTempFactor(const std::string& PDBline)
{
  if(PDBline.length() >= atomTempFactor + 6)
    return (float)atof(PDBline.substr(atomTempFactor, 6).c_str());
  return 0.0;
}

std::vector<unsigned short> PDB::connectedAtoms(const std::string& PDBline)
{
  std::vector<unsigned short> connAtoms;
  if (isCONECTrec(PDBline)) {
    // (1-6) - "CONECT"
    // (7-11), (12-16), ...,(57-61) - serial numbers of connected atoms
    const unsigned short atomIndexLen = 5;
    unsigned short atomIndexStart = 6;
    unsigned short lineLen = (unsigned short)PDBline.length();
    if (lineLen > 61) lineLen = 61;

    while ((atomIndexStart+atomIndexLen) <= lineLen) {
      //char atomIndex[atomIndexLen+1];
      //strncpy(atomIndex, PDBline + atomIndexStart, atomIndexLen);
      //atomIndex[atomIndexLen] = '\0';
      short sh_atomIndex = (short)atoi(PDBline.substr(atomIndexStart, atomIndexLen).c_str());
      if (sh_atomIndex > 0) connAtoms.push_back(sh_atomIndex);
      atomIndexStart += atomIndexLen;
    }
  }
  return connAtoms;
}
