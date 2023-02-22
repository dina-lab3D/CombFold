#ifndef AminoAcid_hpp
#define AminoAcid_hpp

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "PDB.h"
#include "RigidTrans3.h"

class Atom;

enum RESIDUE_INDEX {
    ARG,LYS,GLU,GLN,ASP,ASN,PRO,GLY,SER,THR,
    HIS,ALA,ILE,LEU,VAL,PHE,CYS,MET,TYR,TRP,OTHER,GAP,MAX_AA
};

template<class AtomT = Atom>
/*
CLASS
  AminoAcid

  Container that stores Atom's of amino acid.
  In addition the class provides a general information about amino acid residues.

KEYWORD
  AminoAcid, SeqScore

AUTHORS
    Maxim Shatsky (maxshats@math.tau.ac.il)

    copyright: SAMBA Group , Tel-Aviv Univ. Israel, 1997, 1999, 2002

CHANGES LOG
   19.12.02 Changed to inherit from vector<Atom>. Added function 'parsePDB'.
   16.01.03 Added new residue index OTHER, to support PCA, PTR ...
   10.05.05 parsePDB was changed to differentiate between residues with the same residueId but different residueSeqId
            (for example, 45 and 45A)

GOALS
   Container that stores Atom's of amino acid.
   In addition, AminoAcid class should help to construct user-defined classes
   that deal with amino acid residues. For example it provides one-to-one mapping
   function from character that identifies the amino acid residue to
   index which defined in range of [0,19] (see USAGE for more details),for
   this purpose class provides global enum :
   RESIDUE_INDEX {ARG,LYS,GLU,GLN,ASP,ASN,PRO,GLY,SER,
   THR,HIS,ALA,ILE,LEU,VAL,PHE,CYS,MET,TYR,TRP,OTHER,MAX_AA};


   AminoAcid can provide statistical information about amino acid residues,
   like "Frequency of Occurrence in Proteins".
USAGE
 See SeqScore.
 Complexity of all functions is O(1) besides:
 parsePDB - O(n)
*/
class AminoAcid
  : public std::vector<AtomT> {
private:
  static const float m_fFreqOfOccurrence[];
  static const RESIDUE_INDEX m_indexes[];
  static const char m_revIndexes[];
public:

  ////returns CA atom position (-1 if not found)
  int getCA() const;

  ////returns N atom position (-1 if not found)
  int getN() const;

  ////returns C atom position (-1 if not found)
  int getC() const;

  ////returns O atom position (-1 if not found)
  int getO() const;

  ////returns CA atom position (-1 if not found)
  int getCB() const;

  //// Applies a rigid transformation on all atom coordinates.
  AminoAcid& operator*=(const RigidTrans3 &rt);

  ////Returns frequency of occurrence of each amino acid residue in proteins
  //taken from M.H.Klapper,Biochem.Biophys.Res.Commun.78:1018-1024,1977
  static float getFreqOfOccurrence(RESIDUE_INDEX index);

  ////Receives one-letter abbreviation of amino acid residue.
  //Returns index as appears in enum RESIDUE_INDEX.
  inline static RESIDUE_INDEX getResidueTypeIndx(const char cRes){
    if ( cRes == '-' || cRes == '.' ) return GAP;

    int indx=(int)cRes-(int)'A';

    if (indx<0 || indx>25) { //RESIDUE_INDEX AminoAcid::m_indexes[26]
      std::cerr<<"SequenceScore: unknown residue type:"<< cRes <<std::endl;
      return(RESIDUE_INDEX) -1;
    } else
      return m_indexes[indx];
  }

  ////Receives RESIDUE_INDEX  and returns one-letter abbreviation of amino acid residue.
  inline static char getResidueType(const RESIDUE_INDEX indx){
    return m_revIndexes[indx];
  }

  ////Parses PDB molecule container and extracts AminoAcid objects
  static std::vector<AminoAcid<AtomT> > parsePDB(const std::vector<AtomT>& i_contAtom);

  ////Parses PDB file and extracts AminoAcid objects
  //possible to extract even 'alt' residues (alternative residue entries having the same id)
  static std::vector<AminoAcid<AtomT> > parseFilePDB(const std::string& pdbFile);

  template<class T>
  friend std::ostream& operator<<(std::ostream& s, const AminoAcid& pm);
};

template<class AtomT>
std::vector<AminoAcid<AtomT> > AminoAcid<AtomT>::parsePDB(const std::vector<AtomT>& i_contAtom) {

  std::vector<AminoAcid<AtomT> > o_vAA;

  //residue index in order to differentiate between
  //different residues
  int iCurrRes;
  std::string sCurrRes;
  char chainCurrRes;

  if (i_contAtom.size()==0)
    return o_vAA;

  iCurrRes=i_contAtom[0].residueIndex();
  sCurrRes=i_contAtom[0].residueSequenceID();
  chainCurrRes=i_contAtom[0].chainId();

  AminoAcid vRes;

  for (int i=0;i<(int)i_contAtom.size();i++) {

    if (i_contAtom[i].residueIndex()!=iCurrRes || i_contAtom[i].residueSequenceID()!=sCurrRes ||
	i_contAtom[i].chainId()!=chainCurrRes) {

      o_vAA.push_back(vRes);
      vRes.clear();

      iCurrRes=i_contAtom[i].residueIndex();
      sCurrRes=i_contAtom[i].residueSequenceID();
      chainCurrRes=i_contAtom[i].chainId();
    }

    vRes.push_back(i_contAtom[i]);
  }

  //last residue
  if (vRes.size() > 0) {
    o_vAA.push_back(vRes);
  }

  return o_vAA;
}

template<class AtomT>
std::vector<AminoAcid<AtomT> > AminoAcid<AtomT>::parseFilePDB(const std::string& pdbFile) {
  std::ifstream pdbStream(pdbFile.c_str());

  if (!pdbStream) {
    std::cerr << "File does not exist: " << pdbFile << std::endl;
    exit(1);
  }

  std::vector<AminoAcid<AtomT> > o_vAA;

  //residue index in order to differentiate between
  //different residues
  int iCurrRes=-99999;
  char cCurrAlt='-';
  char cCurrChain='9';

  AminoAcid vRes;

  while (!pdbStream.eof()) {
    std::string line;
    getline(pdbStream, line);
    if (! PDB::isATOMrec(line.c_str()))
      continue;

    char alt=line.c_str()[26];
    AtomT atom(line.c_str());

    if (cCurrAlt=='-')
      cCurrAlt=alt;
    if (iCurrRes == -99999)
      iCurrRes=atom.residueIndex();
    if (cCurrChain == '9')
      cCurrChain = atom.chainId();

    if (atom.residueIndex()!=iCurrRes ||  cCurrAlt!=alt || cCurrChain!=atom.chainId()) {
      o_vAA.push_back(vRes);
      iCurrRes=atom.residueIndex();
      cCurrAlt=alt;
      cCurrChain=atom.chainId();
      vRes.clear();
    }

    vRes.push_back(atom);
  }

  if (vRes.size()>0) {
    o_vAA.push_back(vRes);
  }

  pdbStream.close();

  return o_vAA;
}

//ALA,RES_TMP0,CYS,ASP,GLU,PHE,GLY,HIS,ILE,RES_TMP1,LYS,LEU,MET,ASN,RES_TMP2,PRO,GLN,ARG,SER,THR,RES_TMP3,VAL,TRP,RES_TMP4,TYR
template<class AtomT>
const RESIDUE_INDEX AminoAcid<AtomT>::m_indexes[26]={ALA,OTHER,CYS,ASP,GLU,PHE,GLY,HIS,ILE,OTHER,LYS,LEU,MET,ASN,OTHER,PRO,GLN,ARG,SER,THR,OTHER,VAL,TRP,OTHER,TYR,OTHER};


//ARG,LYS,GLU,GLN,ASP,ASN,PRO,GLY,SER,THR,HIS,ALA,ILE,LEU,VAL,PHE,CYS,MET,TYR,TRP,OTHER,GAP,MAX_AA
template<class AtomT>
const char AminoAcid<AtomT>::m_revIndexes[MAX_AA]={'R','K','E','Q','D','N','P','G','S','T','H','A','I','L','V',
					    'F','C','M','Y','W','.','-'};

//ARG,LYS,GLU,GLN,ASP,ASN,PRO,GLY,SER,THR,HIS,ALA,ILE,LEU,VAL,PHE,CYS,MET,TYR,TRP
template<class AtomT>
const float AminoAcid<AtomT>::m_fFreqOfOccurrence[MAX_AA]={0.047,
                                        0.07,
                                        0.062,
                                        0.039,
                                        0.055,
                                        0.044,
                                        0.046,
                                        0.075,
                                        0.071,
                                        0.06,
                                        0.021,
                                        0.09,
                                        0.046,
                                        0.075,
                                        0.069,
                                        0.035,
                                        0.028,
                                        0.017,
                                        0.035,
                                        0.011,
					0.0,//OTHER
					0.0//GAP
                                        };

template<class AtomT>
int AminoAcid<AtomT>::getCA() const{
  for(int i=0;i<(int)(this->size());++i){
    if((*this)[i].isCA())
    	return i;
  }
  return -1;
}

template<class AtomT>
int AminoAcid<AtomT>::getC() const{
  for(int i=0;i<(int)(this->size());++i){
    if((*this)[i].isC())
    	return i;
  }
  return -1;
}

template<class AtomT>
int AminoAcid<AtomT>::getN() const{
  for(int i=0;i<(int)(this->size());++i){
    if((*this)[i].isN())
    	return i;
  }
  return -1;
}

template<class AtomT>
int AminoAcid<AtomT>::getO() const{
  for(int i=0;i<(int)(this->size());++i){
    if((*this)[i].isO())
    	return i;
  }
  return -1;
}

template<class AtomT>
int AminoAcid<AtomT>::getCB() const{
  for(int i=0;i<(int)(this->size());++i){
    if((*this)[i].isCB())
    	return i;
  }
  return -1;
}

template<class AtomT>
float AminoAcid<AtomT>::getFreqOfOccurrence(RESIDUE_INDEX index){
  return m_fFreqOfOccurrence[index];
}

template<class AtomT>
std::ostream& operator<<(std::ostream& s, const AminoAcid<AtomT>& i_aa){
  for(int i=0;i<(int)i_aa.size();++i)
    s<< i_aa[i] << std::endl;

  return s;
}

template<class AtomT>
AminoAcid<AtomT>& AminoAcid<AtomT>::operator*=(const RigidTrans3 &rt){
 for(int i=0;i<(int)(this->size());++i) (*this)[i]*=rt;
 return *this;
}

#endif
