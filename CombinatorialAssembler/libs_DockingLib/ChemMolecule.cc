#include "ChemMolecule.h"
#include "SASurface.h"

/* ARG,LYS,GLU,GLN,ASP,
   ASN,PRO,GLY,SER,THR,
   HIS,ALA,ILE,LEU,VAL,
   PHE,CYS,MET,TYR,TRP */
const float ChemMolecule::flexProbs[MAX_AA]={0.24, 0.38, 0.14, 0.23, 0.06,
					     0.11, 0.05, 0.05, 0.06, 0.07,
					     0.05, 0.00, 0.13, 0.13, 0.08,
					     0.015,0.02, 0.24, 0.07, 0.02,
					     0.0//OTHER
};

const float ChemMolecule::maxCaDist[MAX_AA]={8.0, 6.3, 5.5, 5.3, 4.2,
					     4.2, 3.0, 2.0, 3.0, 3.0,
					     5.4, 3.0, 4.5, 4.4, 3.1,
					     5.7, 3.5, 5.7, 7.0, 6.6,
					     0.0//OTHER
};

bool compareFragment(ChemMolecule::Fragment f1, ChemMolecule::Fragment f2) {
  if(f1.first != f2.first) { //different chains
    return f1.first < f2.first;
  }
  // same chain
  return f1.second.first < f2.second.first;
}

void ChemMolecule::loadMolecule(std::istream &molFile, const ChemLib& chemLib, const PDB::Selector& selector) {
  Molecule<ChemAtom>::readAllPDBfile(molFile, selector);
  for(Molecule<ChemAtom>::iterator molIter = begin(); molIter != end(); molIter++) {
    if(!molIter->isHydrogen()) {
      const ChemEntry* entry=chemLib.getLibEntry(molIter->residueName(), molIter->type());
      if(entry != NULL) {
	molIter->setRadius(entry->getEntryRadius());
	molIter->setCharge(entry->getEntryCharge());
	molIter->setChemType(entry->getEntryChemType());
      }
    }
  }

  centroid_ = centroid();
  setIsProtein();
}

void ChemMolecule::readAllPDBfile(std::istream &molFile, const PDB::Selector& selector) {
  Molecule<ChemAtom>::readAllPDBfile(molFile, selector);
  computeASAperAtom();
  centroid_ = centroid();
  setIsProtein();
}

float ChemMolecule::computeASAperAtom() {
  SASurface sas(*this, 1.8, 10); //probe radius = 1.8, density = 10
  float asa = sas.computeASAForAtoms();
  Logger::infoMessage() << "ASA of the molecule = " << asa << std::endl;
  return asa;
}

void ChemMolecule::computeSteadyProbs() {
  Vector3 caAtom;
  for(Molecule<ChemAtom>::iterator molIter = begin(); molIter != end(); molIter++) {
    RESIDUE_INDEX resIndex = RESIDUE_INDEX(*molIter);
    if(molIter->isBackbone()) {
      // glycine backbone can be flexible, since it has no side-chain
      if(resIndex == GLY)
	molIter->setSteadyProb(1-flexProbs[resIndex]);
      else {
	molIter->setSteadyProb(1-flexProbs[OTHER]);
	if(molIter->isCA())
	  caAtom = molIter->position();
      }
    } else {
      float distanceFromCa = caAtom.dist(molIter->position());
      float normFactor= distanceFromCa/maxCaDist[resIndex];
      molIter->setSteadyProb(1-flexProbs[resIndex]*normFactor);
    }
  }
}

void ChemMolecule::setIsProtein() {
  for(Molecule<ChemAtom>::iterator molIter = begin(); molIter != end(); molIter++) {
    if(molIter->isCA()) {
      isProtein_ = true;
      return;
    }
  }
  isProtein_ = false;
}

unsigned int ChemMolecule::readBindingSite(const std::string fileName) {
  bindingSites_.push_back(std::set<unsigned int>());
  return readSiteFile(fileName, bindingSites_.back());
}

unsigned int ChemMolecule::readBlockingSite(const std::string fileName) {
  blockingSites_.push_back(std::set<unsigned int>());
  return readSiteFile(fileName, blockingSites_.back());
}

unsigned int ChemMolecule::readSiteFile(const std::string fileName, std::set<unsigned int>& site) {
  std::ifstream file(fileName.c_str());
  if(!file) {
    std::cerr << "Can't open file " << fileName << std::endl;
    Logger::errorMessage() << "Can't open file " << fileName << std::endl;
    exit(1);
  }

  char chainId;
  std::string residueSequenceID;

  while(!file.eof()) {
    std::string line;
    getline(file, line);
    boost::trim(line);
    if(line.length() == 0) continue;
    int entry = -1;

    char record[line.size()+1];
    strcpy(record, line.c_str());

    char* splittedRecord = strtok(record, " ");
    residueSequenceID = std::string(splittedRecord);

    splittedRecord = strtok (NULL, " ");
    if(splittedRecord != NULL) {
      chainId = splittedRecord[0];
    } else {
      chainId = ' ';
    }
    //cerr << residueSequenceID << ";" << chainId << ";" << endl;
    entry = residueEntry(chainId, residueSequenceID);

    if(entry != -1) {
      site.insert(entry);
    } else {
      Logger::errorMessage() << "Invalid line in site file: " << line << std::endl;
      std::cerr << "Invalid line in site file: " << line << std::endl;
      exit(1);
    }
  }

  Logger::infoMessage() << site.size() << " residues were read from site file: " << fileName << std::endl;
  file.close();
  return site.size();
}


const std::set<unsigned int>& ChemMolecule::getBindingSite(unsigned int siteNumber) const {
  const static std::set<unsigned int> emptySet;
  if(bindingSites_.size() > siteNumber) {
    return bindingSites_[siteNumber];
  }
  return emptySet;
}

const std::set<unsigned int>& ChemMolecule::getBlockingSite(unsigned int siteNumber) const {
  const static std::set<unsigned int> emptySet;
  if(blockingSites_.size() > siteNumber) {
    return blockingSites_[siteNumber];
  }
  return emptySet;
}

const ChemAtom& ChemMolecule::getChemAtom(int atomIndex) const {
  for(Molecule<ChemAtom>::const_iterator molIter = begin(); molIter != end(); molIter++) {
    if((int)molIter->atomIndex() == atomIndex) {
      return *molIter;
    }
  }
  Logger::errorMessage() << "Input Error: Can't find atom with atom index " << atomIndex << "!" << std::endl;
  std::cerr << "Input Error: Can't find atom with atom index " << atomIndex << "!" << std::endl;
  exit(1);
}

double ChemMolecule::getVolume() const {
  double volume = 0.0;
  double c = (4.0/3.0)*M_PI;
  for(Molecule<ChemAtom>::const_iterator molIter = begin(); molIter != end(); molIter++) {
    double r = molIter->getRadius();
    volume += c * r * r * r;
  }
  return volume;
}

const std::vector<ChemMolecule::Fragment>& ChemMolecule::getFragments() const {
  if(fragments_.size() == 0) {
    ChemMolecule* nonConstThis = const_cast<ChemMolecule*>(this);
    nonConstThis->computeFragments();
  }
  return fragments_;
}

void ChemMolecule::computeFragments() {
  // calculate endpoints
  char currChain;
  int firstResIndex, prevResIndex;
  bool currChainSet = false;
  for(auto i = begin(); i != end(); i++) {
    if(!i->isCA() || i->getAtomEntryType() == HETATM) continue;
    char chain = i->chainId();
    int resIndex = i->residueIndex();
    //if(i->getAtomEntryType() == HETATM) continue;
    // one more residue of the same chain - advance
    if(currChainSet && currChain == chain &&
       resIndex > prevResIndex) {
      prevResIndex = resIndex;
    } else { // new chain
      if(currChainSet) { // save currChain
        ResidueRange range(firstResIndex, prevResIndex);
        Fragment fragment(std::make_pair(currChain, range));
        fragments_.push_back(fragment);
      }
      // update
      currChain = chain;
      firstResIndex = prevResIndex = resIndex;
      currChainSet = true;
    }
  }
  // save last fragment
  if(currChainSet) { // save currChain
    ResidueRange range(firstResIndex, prevResIndex);
    Fragment fragment(std::make_pair(currChain, range));
    fragments_.push_back(fragment);
  }

  std::sort(fragments_.begin(), fragments_.end(), compareFragment);

  for(int i=0; i<(int)fragments_.size(); i++) {
    std::cout << "Fragment " << i << " chainId "
              << fragments_[i].first << " range "
              << fragments_[i].second.first << ":"
              << fragments_[i].second.second << std::endl;
  }
}

void ChemMolecule::addMol2Type() {


  for(auto i = begin(); i != end(); i++) {
    const char* type = i->type();
    char res = i->residueType();
    i->setMol2Type(UNK_MOL2_TYPE);
    if(i->isN()) i->setMol2Type(N2);
    if(i->isC()) i->setMol2Type(C2);
    if(i->isCA()) i->setMol2Type(C3);
    if(i->isO()) i->setMol2Type(O2); // TODO: the last residue should be O.co2
    if(i->isCB()) i->setMol2Type(C3);

    // carbon
    if(i->isCD() && (res == 'K' || res == 'P' || res == 'R')) i->setMol2Type(C3);
    if(i->isCD() && (res == 'Q' || res == 'E' )) i->setMol2Type(C2);

    if(i->isCG() && (res == 'K' || res == 'P' || res == 'R' || res == 'M' || res == 'L' || res == 'Q' || res == 'E' )) i->setMol2Type(C3);
    if(i->isCG() && (res == 'D' || res == 'N')) i->setMol2Type(C2);
    if(i->isCG() && (res == 'H' || res == 'F' || res == 'W' || res == 'Y')) i->setMol2Type(Car);
    if(i->isCZ() && (res == 'F' || res == 'Y')) i->setMol2Type(Car);
    if(i->isCZ() && res == 'R') i->setMol2Type(Ccat);
    if((i->isCD1() || i->isCD2()) && (res == 'F' || res == 'W' || res == 'Y' || res == 'H')) i->setMol2Type(Car);
    if((i->isCD1() || i->isCD2()) && (res == 'L'  || res == 'I')) i->setMol2Type(C3);
    if(i->isCE() || i->isCG1() || i->isCG2()) i->setMol2Type(C3);
    if(i->isCE1() || i->isCE2() || i->isCE3() || i->isCZ2() || i->isCZ3() || i->isCH2() ) i->setMol2Type(Car);

    // for RNA
    if (i->isSugarCarbon())
      i->setMol2Type(C3);

    // oxygen
    if((i->isOD1() || i->isOE1()) && ((res == 'N' || res == 'Q'))) i->setMol2Type(O2);
    if((i->isOD1() || i->isOE1()|| i->isOD2()|| i->isOE2()) && (res == 'D' || res == 'E')) i->setMol2Type(Oco2);
    if(i->isOH() || i->isOG()|| i->isOG1()) i->setMol2Type(O3);

    // for RNA
    if (i->isSugarOxygen())
      i->setMol2Type(O3);
    if (i->isPhosphateOxygen())
      i->setMol2Type(Oco2);

    // nitrogen
    if(i->isNH1() || i->isNH2() || i->isNE()) i->setMol2Type(Npl3);
    if((i->isND2() && res == 'N') || (i->isNE2() && res == 'Q')) i->setMol2Type(Nam);
    if(((i->isND1() || i->isNE2()) && res == 'H') ||  (i->isNE1() && res == 'W')) i->setMol2Type(Nar);
    if(i->isNZ()) i->setMol2Type(N4);

    // for RNA
    if (i->isNitrogenousBaseAtom()) {
      if (type[1] == 'N') {
        if ((type[2] == '4' && res == 'C') || (type[2] == '6' && res == 'A') || (type[2] == '2' && res == 'G'))
          i->setMol2Type(Npl3);
        else
          i->setMol2Type(Nar);
      }
      else if (type[1] == 'O')
        i->setMol2Type(O2);
      else if (type[1] == 'C')
        i->setMol2Type(Car);
    }

    // sulfur
    if(i->isSD() || i->isSG()) i->setMol2Type(S3);

    // phosphate
    if(i->isP()) i->setMol2Type(P3);

    if(i->getMol2Type() == UNK_MOL2_TYPE)
      std::cout << "Can't assign mol2 type:" << type[0] << type[1] << type[2] << type[3] << ":" << std::endl;
  }
}
