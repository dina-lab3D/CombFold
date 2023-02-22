#ifndef INDEX_MOL_H
#define INDEX_MOL_H

#include <Molecule.h>
#include <AminoAcid.h>
#include <Logger.h>

#include <boost/algorithm/string.hpp>

#include <map>
#include <set>
#include <string>

template<class AtomT> // AtomT must implement: residueSequenceID() and chainId()
class ResidueIndexedMolecule : public Molecule<AtomT> {
 public:
  ResidueIndexedMolecule() : Molecule<AtomT>() {isInitialized_ = false;}

  //// check if the mapping is initialized
  bool isInitialized() const { return isInitialized_; }

  //// clear and init
  void resetResEntries();

  ////
  int residueEntry(char chainId, std::string resSequenceID ) const {
    if (!isInitialized_) {
      ResidueIndexedMolecule<AtomT>* nonConstThis = const_cast<ResidueIndexedMolecule<AtomT>*>(this);
      nonConstThis->resetResEntries();
    }
    auto aPair = resEntries_.find(indexOf(chainId, resSequenceID));
    if (aPair == resEntries_.end())
      return -1;
    return (*aPair).second;
  }

  ////
  int getFirstAtomEntryForResEntry(unsigned int entry) const {
    if (!isInitialized_) {
      ResidueIndexedMolecule<AtomT>* nonConstThis = const_cast<ResidueIndexedMolecule<AtomT>*>(this);
      nonConstThis->resetResEntries();
    }
    auto it = entries2atoms_.find(entry);
    if(it == entries2atoms_.end())
      return -1;
    return it->second;
  }

  ////
  int getFirstAtomEntryForResIndex(char chainId, std::string resSequenceID) const {
    return getFirstAtomEntryForResEntry(residueEntry(chainId, resSequenceID));
  }

  ////
  int getClosestAtomEntryForResIndex(char chainId, std::string resSequenceID,
                                     int maxOffset, int& offset) const;

  //// mapping from atom entry (=index in array) to residue entry (=index in array)
  int atomEntryToResEntry(unsigned int entry) const {
    if (!isInitialized_) {
      ResidueIndexedMolecule<AtomT>* nonConstThis = const_cast<ResidueIndexedMolecule<AtomT>*>(this);
      nonConstThis->resetResEntries();
    }
    const AtomT &atom = (*this)[entry];
    return residueEntry(atom.chainId(), atom.residueSequenceID());
  }

  ////
  unsigned int numberOfRes() const {
    if (!isInitialized_) {
      ResidueIndexedMolecule<AtomT>* nonConstThis = const_cast<ResidueIndexedMolecule<AtomT>*>(this);
      nonConstThis->resetResEntries();
    }
    return resEntries_.size();
  }

private:
  bool isInitialized_;

  // mapping data structure
  std::map<std::string, unsigned int> resEntries_;
  std::map<unsigned int, unsigned int> entries2atoms_;
  std::set<char> chainIDs_; // chain ids of the molecules

private:
  // compute hash map key
  std::string indexOf(char chainId, std::string resSequenceID) const;

  //// initialization
  void initResEntries();
};


template<class AtomT>
std::string ResidueIndexedMolecule<AtomT>::indexOf(char chainId, std::string resSequenceID) const {
  std::string indexStr = resSequenceID;
  boost::trim(indexStr);
  indexStr += "_";
  indexStr += chainId;
  return indexStr;
}

template<class AtomT>
void ResidueIndexedMolecule<AtomT>::resetResEntries() {
  resEntries_.clear();
  entries2atoms_.clear();
  chainIDs_.clear();
  initResEntries();
}

template<class AtomT>
void ResidueIndexedMolecule<AtomT>::initResEntries() {
  std::vector<AminoAcid<AtomT> > aaMol = AminoAcid<AtomT>::parsePDB(*this);
  //  cerr << "aaMol.size() = " << aaMol.size() << endl;
  for(unsigned int residueCounter=0; residueCounter<aaMol.size(); residueCounter++) {
    int CAindex = aaMol[residueCounter].getCA();
    if(CAindex == -1)
      CAindex=0;
    char chainId = aaMol[residueCounter][CAindex].chainId();
    std::string residueSequenceID = aaMol[residueCounter][CAindex].residueSequenceID();
    resEntries_[indexOf(chainId, residueSequenceID)] = residueCounter;
    chainIDs_.insert(chainId);
  }

  unsigned int atomCounter=0;

  for(unsigned int residueCounter=0; residueCounter<aaMol.size(); residueCounter++) {
    int CAindex = aaMol[residueCounter].getCA();
    if(CAindex == -1)
      CAindex=0;
    char chainId = aaMol[residueCounter][CAindex].chainId();
    std::string resSequenceID = aaMol[residueCounter][CAindex].residueSequenceID();
    std::map<std::string, unsigned int>::const_iterator aPair = resEntries_.find(indexOf(chainId, resSequenceID));
    if (aPair != resEntries_.end()) {
      int currEntry = (*aPair).second;
      entries2atoms_[currEntry] = atomCounter+CAindex;
    }
    atomCounter+=aaMol[residueCounter].size();
  }

  isInitialized_ = true;
}

template<class AtomT>
int ResidueIndexedMolecule<AtomT>::getClosestAtomEntryForResIndex(char chainId,
                                                                  std::string resSequenceId,
                                                                  int maxOffset,
                                                                  int& offset) const {
  if (!isInitialized_) {
    ResidueIndexedMolecule<AtomT>* nonConstThis = const_cast<ResidueIndexedMolecule<AtomT>*>(this);
    nonConstThis->resetResEntries();
  }

  if(chainIDs_.find(chainId) == chainIDs_.end())  return -1;

  // 1. look for CA of the given residue
  int atomIndex = getFirstAtomEntryForResIndex(chainId, resSequenceId);
  if(atomIndex != -1)
    return atomIndex;

  // 2. look for +/- offset
  int resIndex = std::stoi(resSequenceId);
  for(int i=1; i<=maxOffset; i++) {
    // try +/- i
    std::string resString = std::to_string(resIndex + i);
    int entry = residueEntry(chainId, resString);
    if(entry != -1) { // match found
      offset = i;
      Logger::warningMessage() << "Can't find residue " << resSequenceId << " chain "
                               << chainId << " using closest " << resString << std::endl;
      return getFirstAtomEntryForResEntry(entry);
    }
    // -i
    resString = std::to_string(resIndex - i);
    entry = residueEntry(chainId, resString);
    if(entry != -1) { // match found
      offset = -i;
      Logger::warningMessage() << "Can't find residue " << resSequenceId << " chain "
                               << chainId << " using closest " << resString << std::endl;
      return getFirstAtomEntryForResEntry(entry);
    }
  }
  Logger::warningMessage() << "Can't find residue " << resSequenceId
                           << " chain " << chainId << std::endl;
  return -1;
}

#endif
