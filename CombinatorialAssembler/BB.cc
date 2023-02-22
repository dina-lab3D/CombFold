#include "BB.h"

#include <Common.h>
#include <connolly_surface.h>
// #include <ShuoSurface.h>

// static declarations
float BB::penetrationThreshold_ = 3.5;
ScoreParams BB::scoreParams_;
float BB::gridResolution_(0.0);
float BB::gridMargins_(0.0);

namespace {
bool isNonPolr(const char res, const char* atom) {
  //cout << " isNonPolr " << res << " " << atom << endl;
  if (atom[1] == 'S') return true;
  if (atom[1] != 'C') return false;
  if ((atom[2] == 'A' || atom[2] == 'B' ) && atom[3] == ' ') return true;
  if ( atom[2] == 'H') return true;
  if ( atom[2] == 'G') {
    if (res == 'N') return false;
    return true;
  }
  if ( atom[2] == 'D') {
    if (res == 'Q') return false;
    return true;
  }
  if ( atom[2] == 'E') {
    if (res == 'H') return false;
    return true;
  }
  if ( atom[2] == 'Z') {
    if (res == 'R') return false;
    return true;
  }
  if ( atom[2] == 'Z') {
    if (res == 'R') return false;
    return true;
  }
  return false;
}
};


BB::BB(int id, const std::string pdbFileName, int groupID, int groupSize, const ChemLib& lib) :
  id_(id), groupId_(groupID), groupSize_(groupSize), pdbFileName_(pdbFileName) {
  // read atoms
  Common::readChemMolecule(pdbFileName_, allAtoms_, lib);
  std::cout << "Done reading ChemMolecule " << allAtoms_.size() << std::endl;
  // read backbone atoms
  std::ifstream pdb2(pdbFileName_);
  backBone_.readPDBfile(pdb2, PDB::BBSelector());
  pdb2.close();
  cm_ = backBone_.centroid();
  numOfAtoms_ = allAtoms_.size();

  // read CA atoms
  std::ifstream pdb3(pdbFileName_);
  caAtoms_.readPDBfile(pdb3, PDB::CAlphaSelector());
  pdb3.close();


  // compute ms surface
  msSurface_ = get_connolly_surface(allAtoms_, 10, 1.8);
  std::cout << "Surface size " << msSurface_.size() << std::endl;

  // compute grid
  grid_ = new BBGrid(msSurface_, gridResolution_, gridMargins_, 1.5);
  grid_->computeDistFromSurface(msSurface_);
  grid_->markTheInside(allAtoms_);
  grid_->markResidues(backBone_);
  std::cout << "Done compute grid " << pdbFileName_ << std::endl;

  // compute shuo surface from Connolly MS
  /*ShuoSurface shuo(msSurface_, allAtoms_, *grid_);
  shuo.computeShuoSurface(surface_);
  std::cout << "Shuo surface done " << surface_.size() << std::endl;
  */
  std::vector<Atom*> atomsMap;
  int j = 0;
  atomsMap.push_back(&(*allAtoms_.begin()));
  for(ChemMolecule::iterator i= allAtoms_.begin(); i != allAtoms_.end(); i++) {
    atomsMap.push_back(&(*i));
  }
  for(Surface::iterator it = surface_.begin(); it != surface_.end(); it++) {
    bool nonPolar = isNonPolr(atomsMap[it->atomIndex(0)]->residueType(),
                              atomsMap[it->atomIndex(0)]->type());
    if (nonPolar) {
      //cerr << " nonPolar " << j << endl;
      j++;
    }
    ipSet_.push_back(InterestPoint(it->position(),
                                   it->position() +
                                   grid_->calcVolumeNormal(it->position(), 4.0)*PROB_RADIUS,
                                   nonPolar));

  }
  std::cout << "# nonPolar points" << j << std::endl;

  dockScore_ = new GeomScore(msSurface_, scoreParams_.weights, scoreParams_.max_penetration, 0.5);
  dockScore_->buildTree();

  computeFragments(); // get the endpoints

  for (Molecule<Atom>::const_iterator it = caAtoms_.begin(); it != caAtoms_.end(); it++) {
      resIndexToCAAtom[it->residueIndex()] = *it;
  }

  std::cout << " done reading BB " << pdbFileName_.c_str() << std::endl;
}


const Atom& BB::getAtomByResId(unsigned int resId) const {
    return resIndexToCAAtom.at(resId);
}

void BB::computeFragments() {
  // calculate endpoints
  char currChain;
  int firstResIndex, prevResIndex;
  bool currChainSet = false;
  for(auto i = allAtoms_.begin(); i != allAtoms_.end(); i++) {
    // only CA atoms are considered
    if(!i->isCA()) continue;
    char chain = i->chainId();
    int resIndex = i->residueIndex();
    // one more residue of the same chain - advance
    if(currChainSet && currChain == chain) {
      prevResIndex = resIndex;
    } else { // new chain
      if(currChainSet) { // save currChain
        ResidueRange range(firstResIndex, prevResIndex);
        fragmentEndpoints_.push_back(std::make_pair(currChain, range));
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
    fragmentEndpoints_.push_back(std::make_pair(currChain, range));
  }

  for(int i=0; i<(int)fragmentEndpoints_.size(); i++) {
    std::cout << "Fragment " << i << " chainId "
              << fragmentEndpoints_[i].first << " range "
              << fragmentEndpoints_[i].second.first << ":"
              << fragmentEndpoints_[i].second.second << std::endl;
  }
}

void BB::getChainConnectivityConstraints(const BB& otherBB,
                                         std::vector<std::pair<char, std::pair<int, int>>>& constraints) const {
  for(int i=0; i<(int)fragmentEndpoints_.size(); i++) {
    char chainId1 = fragmentEndpoints_[i].first;
    int resIndex1N = fragmentEndpoints_[i].second.first;
    int resIndex1C = fragmentEndpoints_[i].second.second;

    for(int j=0; j<(int)otherBB.fragmentEndpoints_.size(); j++) {
      char chainId2 = otherBB.fragmentEndpoints_[j].first;
      if(chainId1 != chainId2) continue;
      int resIndex2N = otherBB.fragmentEndpoints_[j].second.first;
      int resIndex2C = otherBB.fragmentEndpoints_[j].second.second;
      if(resIndex1N < resIndex2N) { // add constraint on 1C and 2N
        constraints.push_back(std::make_pair(chainId1, std::make_pair(resIndex1C, resIndex2N)));
      } else { // add constraint on 2C and 1N
        constraints.push_back(std::make_pair(chainId1, std::make_pair(resIndex1N, resIndex2C)));
      }
    }
  }
}


//  float BB::scoreForTrans(Score &score, RigidTrans3 trans, BB &other) {
//    float RADIUS = 1.5;
//    float RADIUS2 = RADIUS * RADIUS;
//    GeomHash <Triangle, TransformationAndScore*> *gh = trans_[ENTRY(id_, other.id_)];
//    //cerr << "GH " << gh << endl;
//    Triangle t(t_);
//    t *= trans;
//    HashResult<TransformationAndScore*> result;
//    gh->query(t, RADIUS, result); //!!!!!
//    float lscore = 0;
//    TransformationAndScore *bestTrans(NULL);
//    for (HashResult<TransformationAndScore*>::iterator it = result.begin(); it != result.end(); it++) {
//      float dist2 = t.dist2((*it)->t_);
//      //cerr << "trans " << trans << " t " << t << " (*it)->refFrame_ " << (*it)->refFrame_ << " (*it)->t_ " << (*it)->t_ << " dist2 " << dist2 << " score " << (*it)->score_.totalScore_ << endl;
//      if (dist2 < RADIUS2 &&
//          lscore < (*it)->score_.totalScore_) {
//          lscore = (*it)->score_.totalScore_;
//          bestTrans = *it;
//      }
//    }
//    if (bestTrans != NULL) {
//      for (int i = 0; i < NO_OF_RANGES; i++) {
//        score.ps_[i].count_ += bestTrans->score_.ps_[i].count_;
//        score.ps_[i].surface_ += bestTrans->score_.ps_[i].surface_;
//      }
//      score.totalScore_ += bestTrans->score_.totalScore_;
//      lscore = bestTrans->score_.totalScore_;
//    } else {
//      //cerr << " for trans " << trans  << " no entry " << " bb1 " << id_ << " bb2 " << other.id_ << endl;
//    }
//    return lscore;
//  }

bool BB::isPenetrating(const RigidTrans3 &trans, const BB &other, float threshold) const {
  //  for (vector<SurfaceAtom>::iterator it = other.sas_.begin(); it != other.sas_.end(); it++) {
    //float penetration = getDistFromSurface(trans * it->pos_) - it-> radius_;
  //float max = 1000;
  for (Surface::const_iterator it = other.surface_.begin(); it != other.surface_.end(); it++) {
    float penetration = getDistFromSurface(trans * it->position());
    if (threshold > penetration) {
      //cerr << "trans " << trans << " penetrates " << penetration << endl;
      return true;
    }
    //
    //if (max > penetration) {
    //  max = penetration;
    //}
  }
  //cerr << id_ << " " << other.id_ << " maxPen " << max << endl;
  return false;
}

float BB::maxPenetration(const RigidTrans3 &trans, const BB &other) const {
  //  for (vector<SurfaceAtom>::iterator it = other.sas_.begin(); it != other.sas_.end(); it++) {
    //float penetration = getDistFromSurface(trans * it->pos_) - it-> radius_;
  float max = 1000;
  for (Surface::const_iterator it = other.surface_.begin(); it != other.surface_.end(); it++) {
    float penetration = getDistFromSurface(trans * it->position());
    if (max > penetration) {
      max = penetration;
    }
  }
  return max;
}

void BB::computeRanges(const std::vector<std::shared_ptr<const BB> >&bbs, int margin) const {
  for (size_t i = 0; i < bbs.size(); i++) {
    std::pair<int, int> range;// = computeRange(*bbs[i], margin);
    if(groupId_ == bbs[i]->groupId_) // same group
      range = std::make_pair(1, groupSize_ + margin);
    else
      range = std::make_pair(std::max(groupSize_,
                                      bbs[i]->groupSize_) - margin,
                             trans_.size());
    transformationRanges_.push_back(range);
  }
}

bool BB::isIdent(const BB& otherBB) const {
    if(getNumOfAtoms() != otherBB.getNumOfAtoms()) {
        std::cout << "different num of atoms" << std::endl;
        return false;
    }
    for (unsigned int i=0; i < allAtoms_.size(); i++) {
        if(!(allAtoms_[i].position() - otherBB.allAtoms_[i].position()).isZero()){
            std::cout << "different atom position" << i << std::endl;
            return false;
        }

    }
    std::cout << "checking transforms" << std::endl;
    if(trans_.size() != otherBB.trans_.size())
        return false;
    for (unsigned int i=0; i < trans_.size(); i++) {
        if (i == id_ || i == otherBB.id_)
            continue;

        if(trans_[i].size() != otherBB.trans_[i].size()){
            std::cout << "different number of trans " << i << " " << trans_[i].size() << " != " << otherBB.trans_[i].size() << std::endl;
            return false;
        }

        for (unsigned int j=0; j < trans_[i].size(); j++) {
            if (trans_[i][j]->score() != otherBB.trans_[i][j]->score()){
                std::cout << "different trans score" << i << std::endl;
                return false;
            }

            // TODO: better way to compare transformations (this has problem with more than 2 copies of chain
//            if(!(trans_[i][j]->trans().translation() - otherBB.trans_[i][j]->trans().translation()).isZero()){
//                std::cout << "different trans translation" << i << std::endl;
//                return false;
//            }
        }
    }
    return true;
}
