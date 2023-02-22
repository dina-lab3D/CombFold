#include "SuperBB.h"

#include "HierarchicalFold.h"

unsigned int SuperBB::counters_[100][100] = {0};
//unsigned int SuperBB::failsCounter_ = 0;
//unsigned int SuperBB::failsCounters_[20][20] = {0};

SuperBB::SuperBB(std::shared_ptr<const BB> bb, unsigned int numberOfBBs) :
        numberOfBBs_(numberOfBBs), backBonePen_(0), newScore_(0.0), transScore_(0), weightedTransScore_(100),
        contact_(numberOfBBs * numberOfBBs, false), rmsd_(0.0), yetAnother_(0) {

    id_ = bb->id_;
    bbs_.push_back(bb);
    Vector3 v(0, 0, 0);
    Matrix3 M(1);
    RigidTrans3 T(M, v);
    trans_.push_back(T);
    size_ = 1;
    bitIDS_ = bb->bitId();
    numOfAtoms_ = bb->numOfAtoms_;
    sps_ = bb->ipSet_.size();
}

float SuperBB::compactness(int numOfAtoms, int sps) {
    float f_numOfAtoms = numOfAtoms * 3.0 / 100.0;
    float f_sps = ((float) sps) / 1000.0;
    //cerr << " numOfAtoms " << numOfAtoms << " sps " << sps << endl;
    //cerr << " numOfAtoms " << f_numOfAtoms << " sps " << f_sps << endl;
    return (f_numOfAtoms * f_numOfAtoms) / (f_sps * f_sps * f_sps);

}

void SuperBB::join(const RigidTrans3 &trans, const SuperBB &other, unsigned int burriedSps) {
    // now joining is simple
    for (unsigned int j = 0; j < other.size_; j++) {
        bbs_.push_back(other.bbs_[j]);
        trans_.push_back(trans * other.trans_[j]);
    }
    size_ += other.size_;
    child1_ = bitIDS_;
    child2_ = other.bitIDS_;
    bitIDS_ |= other.bitIDS_;
    id_ = std::min(other.id_, id_);
    sps_ += other.sps_ - burriedSps;
    numOfAtoms_ += other.numOfAtoms_;
    counters_[size_][id_]++;
}

/*
void SuperBB::reportCounters() {
  int total = 0;
  int totalFails = 0;
  for (unsigned int size = 2; size <= numberOfBBs_; size++) {
    int level = 0;
    int failLevel = 0;
    for (unsigned int id = 0; id <= numberOfBBs_ - size; id++) {
      std::cerr << " counters_["<<size<<"]["<<id<<"] = " << counters_[size][id];
      std::cerr << " failsCounters_["<<size<<"]["<<id<<"] = " << failsCounters_[size][id];
      level += counters_[size][id];
      failLevel += failsCounters_[size][id];
    }
    std::cerr << " total of size " << size << " is: " << level << " fails " << failLevel << std::endl;
    total += level;
    totalFails += failLevel;
  }
  std::cerr << " total " << total << std::endl;
  std::cerr << " totalFails " << totalFails << std::endl;
}
*/

float getWeightedTransScore(std::vector<FoldStep> steps, std::vector <std::shared_ptr<const BB>> bbs) {
    unsigned int totalCa = 0;
    float totalScore = 0;
    std::map<unsigned int, unsigned int> bbIdToCaLength;
    for(std::shared_ptr<const BB> bb: bbs) {
        bbIdToCaLength[bb->getID()] = bb->getNumOfAtoms();
    }

    for(FoldStep step : steps) {
        std::vector<unsigned int> right;
        std::vector<unsigned int> left;

        right.push_back(step.i_);
        left.push_back(step.j_);

        std::vector<FoldStep> notUsedSteps(steps);

        while(!notUsedSteps.empty()) {
            FoldStep step2 = notUsedSteps.back();
            notUsedSteps.pop_back();
            if(step2.i_ == step.i_ && step2.j_ == step.j_)
                continue;
            if(std::find(right.begin(), right.end(), step2.i_) != right.end() && std::find(right.begin(), right.end(), step2.j_) == right.end()) {
                right.push_back(step2.j_);
            } else if(std::find(right.begin(), right.end(), step2.j_) != right.end() && std::find(right.begin(), right.end(), step2.i_) == right.end()) {
                right.push_back(step2.i_);
            } else if(std::find(left.begin(), left.end(), step2.i_) != left.end() && std::find(left.begin(), left.end(), step2.j_) == left.end()) {
                left.push_back(step2.j_);
            } else if(std::find(left.begin(), left.end(), step2.j_) != left.end() && std::find(left.begin(), left.end(), step2.i_) == left.end()) {
                left.push_back(step2.i_);
            } else {
                notUsedSteps.insert(notUsedSteps.begin(), step2);
            }
        }

        unsigned int rightSize = 0;
        unsigned int leftSize = 0;

        for (unsigned int elem : right)
            rightSize += bbIdToCaLength[elem];

        for (unsigned int elem : left)
            leftSize += bbIdToCaLength[elem];

        unsigned int usedSize = rightSize;
        if (rightSize > leftSize)
            usedSize = leftSize;

        totalScore += usedSize * step.tScore_;
        totalCa += usedSize;
    }
    return totalScore / totalCa;

}

void SuperBB::join(const RigidTrans3 &trans, const SuperBB &other, unsigned int burriedSps,
                   int bbPen, float newScore, FoldStep &step, float transScore) {
    // now joining is simple
    for (unsigned int j = 0; j < other.size_; j++) {
        bbs_.push_back(other.bbs_[j]);
        trans_.push_back(trans * other.trans_[j]);
    }
    size1_ = size_;
    size2_ = other.size_;
    size_ += other.size_;
    child1_ = bitIDS_;
    child2_ = other.bitIDS_;
    bitIDS_ |= other.bitIDS_;
    id_ = std::min(other.id_, id_);
    sps_ += other.sps_ - burriedSps;
    backBonePen_ = bbPen;
    newScore_ += other.newScore_ + newScore;

    int totalAtomNum = numOfAtoms_ + other.numOfAtoms_;

    numOfAtoms_ = totalAtomNum;
    transScore_ += other.transScore_ + transScore;

    if (!other.foldSteps_.empty()) {
        foldSteps_.insert(foldSteps_.end(), other.foldSteps_.begin(), other.foldSteps_.end());
    }
    foldSteps_.push_back(step);
    counters_[size_][id_]++;

    weightedTransScore_ = getWeightedTransScore(foldSteps_, bbs_);
}


bool SuperBB::isPenetrating(const RigidTrans3 &trans, const SuperBB &other, float threshold) const {
    for (unsigned int i = 0; i < size_; i++) {
        RigidTrans3 t = (!trans_[i]) * trans;
        for (unsigned int j = 0; j < other.size_; j++) {
            if (bbs_[i]->isPenetrating(t * other.trans_[j], *other.bbs_[j], threshold)) return true;
        }
    }
    return false;
}

float SuperBB::calcCompact() {
    int sps = 0;
    int numOfAtoms = 0;
    for (unsigned int i = 0; i < size_; i++) {
        numOfAtoms += bbs_[i]->numOfAtoms_;
        for (Surface::const_iterator it = bbs_[i]->surface_.begin(); it != bbs_[i]->surface_.end(); it++) {
            Vector3 p = trans_[i] * (it->position() + (it->normal() * 1.8));
            bool isSP = true;
            for (unsigned int j = 0; j < size_; j++) {
                if (j == i) continue;
                float penetration = bbs_[j]->getDistFromSurface((!trans_[j]) * p);
                if (penetration <= 0) {
                    isSP = false;
                    break;
                }
            }
            if (isSP) sps++;
        }
    }
    float f_numOfAtoms = numOfAtoms * 3.0 / 100.0;
    float f_sps = ((float) sps) / 1000.0;
    std::cerr << " numOfAtoms " << numOfAtoms << " sps " << sps << std::endl;
    std::cerr << " numOfAtoms " << f_numOfAtoms << " sps " << f_sps << std::endl;

    return (f_numOfAtoms * f_numOfAtoms) / (f_sps * f_sps * f_sps);
}

void SuperBB::calcFinalScore() {
    backBonePen_ = 0;
    maxPen_ = 0.0;
    numOfBuried_ = 0;
    numOfXX_ = 0;
    numOfHydBuried_ = 0;
    newScore_ = 0;
    multPen_ = 0;
    singlePen_ = 0;
    for (unsigned int i = 0; i < size_; i++) {
        const BB &bb1 = *bbs_[i];
        std::vector<RigidTrans3> tr;
        for (unsigned int j = 0; j < size_; j++) {
            tr.push_back((!trans_[j]) * trans_[i]);
        }
        for (unsigned int sp = 0; sp < bb1.ipSet_.size(); sp++) { // loop over ips
            const InterestPoint &ip = bb1.ipSet_[sp];
            bool unburried = true;
            bool xx = false;
            int countPen = 0;

            for (unsigned int j = 0; j < size_; j++) {
                if (i == j) continue;
                const BB &bb2 = *bbs_[j];
                float dist = bb2.getDistFromSurface(tr[j] * ip.position_); // check penetration
                if (dist < maxPen_) {
                    maxPen_ = dist;
                }
                if (dist <= 0) countPen++;
                float dist1 = bb2.getDistFromSurface(tr[j] * ip.probePosition_); // check is burried
                if (dist1 > dist && dist < 0) {
                    numOfXX_++;
                    xx = true;
                }
                if (unburried && !xx && dist <= 1.0) { // if burried
                    unburried = false;
                }
            }
            if (!unburried && !xx) {
                numOfBuried_++;
                if (ip.hydrophobic_) numOfHydBuried_++;
            }

            if (countPen > 1) multPen_++;
            if (countPen == 1) singlePen_++;

        }
        //std::cout << " multPen_ " << multPen_ << " singlePen_ " << singlePen_ << std::endl;
        for (unsigned int j = 0; j < size_; j++) {
            if (i == j) continue;
            const BB &bb2 = *bbs_[j];
            //newScore_ += bb1.dockScore_->score(tr[j], *bb2.grid_);
            bb1.dockScore_->setGrid(bb2.grid_);
            newScore_ += bb1.dockScore_->score(tr[j]);
            for (Molecule<ChemAtom>::const_iterator it = bb1.backBone_.begin(); it != bb1.backBone_.end(); it++) {
                Vector3 v = tr[j] * it->position();
//                if (bb2.getDistFromSurface(v) < 0) {
                  if (bb2.getDistFromSurface(v) < -1.0) {
                    //cerr << "a_penetration: " << bb2.getDistFromSurface(v) << endl;
                    if (bb2.grid_->getResidueEntry(v) < 0) {
                        //cerr << "b_penetration: " << backBonePen_ << endl;
                        backBonePen_++;
                    }
                }
            }
        }
    }
}

// RMSD relative to input position
double SuperBB::calcRmsd(RigidTrans3 &trans, float &compact2) const {
  Match match;
  Molecule<Atom> A, B;
  int j = 0;
  Vector3 cm(0, 0, 0);
  for (unsigned int i = 0; i < size_; i++) {
    const Molecule<ChemAtom> &backBone = bbs_[i]->backBone_;
    const RigidTrans3 &tr = trans_[i];
    for (Molecule<ChemAtom>::const_iterator it = backBone.begin(); it != backBone.end(); it++) {
      if (it->type()[1] == 'C' && it->type()[2] == 'A') {
        Vector3 pos(tr * it->position());
        Atom a(*it);//; SET_ATOM_POS(it->position());
        Atom b(a);
        b *= tr;
        A.add(a);
        B.add(b);
        match.add(j, j, 1);
        j++;
        cm += pos;
      }
    }
  }
  cm /= j;
  //cerr << "cm " << cm<< endl;
  match.calculateBestFit(A, B);
  trans = match.rigidTrans();

  float ms = 0;
  for (unsigned int i = 0; i < size_; i++) {
    const Molecule<ChemAtom> &backBone = bbs_[i]->backBone_;
    const RigidTrans3 &tr = trans_[i];
    for (Molecule<ChemAtom>::const_iterator it = backBone.begin(); it != backBone.end(); it++) {
      if (it->type()[1] == 'C' && it->type()[2] == 'A') {
        Vector3 pos(tr * it->position() - cm);
        float dist2 = pos.norm2();
        ms += dist2;
      }
    }
  }
  ms /= j;
  compact2 = sqrt(ms);
  return match.rmsd();
}


// RMSD between two SBBs (assuming same BBs)
double SuperBB::calcRmsd(const SuperBB &other, std::vector<std::vector<unsigned int>>& identGroups) const {
    std::vector<std::vector<unsigned int>> presentIdentGroups;

    for (std::vector<unsigned int> identGroup: identGroups){
        std::vector<unsigned int> possiblyRelevantIdentGroup;
        for (unsigned int i: identGroup){
            if((bitIDS_ & (1 << i)) != 0)
                possiblyRelevantIdentGroup.push_back(i);
        }
        if(possiblyRelevantIdentGroup.size() >= 2)
            presentIdentGroups.push_back(possiblyRelevantIdentGroup);
    }

    std::map<unsigned int, unsigned int> bbIdToThisBBIndex;
    std::map<unsigned int, unsigned int> bbIdToOtherBBIndex;
    std::vector<unsigned int> bbIdsInThis;
    for(unsigned int i =0; i < bbs_.size(); i++) {
        bbIdsInThis.push_back(bbs_[i]->getID());
        bbIdToThisBBIndex[bbs_[i]->getID()] = i;
        for(unsigned int j = 0; j < size_; j++) {
            if(bbs_[i]->getID() == other.bbs_[j]->getID()) { // same BB
                bbIdToOtherBBIndex[bbs_[i]->getID()] = j;
                break;
            }
        }
    }

    // base rmsd find
    Match cmMatch;
    Molecule<Vector3> cmA, cmB;

    for(unsigned int bitId: bbIdsInThis) {
        Vector3 a = trans_[bbIdToThisBBIndex[bitId]]*bbs_[bbIdToThisBBIndex[bitId]]->cm_;
        cmA.add(a);
        Vector3 b = other.trans_[bbIdToOtherBBIndex[bitId]] * other.bbs_[bbIdToOtherBBIndex[bitId]]->cm_;
        cmB.add(b);
    }
    for(unsigned int i = 0; i < size_; i++) {
        cmMatch.add(i,i);
    }
    cmMatch.calculateBestFit(cmA, cmB);
    float foundRMSD = cmMatch.rmsd();

    // std::cout << "trying to cluster " << numOfBuried_ << " " << other.numOfBuried_ << " " << foundRMSD << std::endl;

    int MAX_ITER = 10;
    // try replacing chains with each other
    for(int i = 0; i < MAX_ITER; i++){
	    //std::cout << "iter " << i << std::endl;
	    bool somethingChanged = false;
        for (std::vector<unsigned int> identGroup: presentIdentGroups){
            for(unsigned int i = 0; i < identGroup.size(); i++) {
                for(unsigned int j = i + 1; j < identGroup.size(); j++){
                    Match cmMatch;
                    Molecule<Vector3> cmA, cmB;

                    for(unsigned int bitId: bbIdsInThis) {
                        Vector3 a;
                        if(bitId == identGroup[i]) {
                            unsigned int replacedBitId = identGroup[j];
                            a = trans_[bbIdToThisBBIndex[replacedBitId]]*bbs_[bbIdToThisBBIndex[replacedBitId]]->cm_;
                        } else if(bitId == identGroup[j]) {
                            unsigned int replacedBitId = identGroup[i];
                            a = trans_[bbIdToThisBBIndex[replacedBitId]]*bbs_[bbIdToThisBBIndex[replacedBitId]]->cm_;
                        } else {
                            a = trans_[bbIdToThisBBIndex[bitId]]*bbs_[bbIdToThisBBIndex[bitId]]->cm_;
                        }
                        cmA.add(a);
                        Vector3 b = other.trans_[bbIdToOtherBBIndex[bitId]] * other.bbs_[bbIdToOtherBBIndex[bitId]]->cm_;
                        cmB.add(b);
                    }
                    for(unsigned int i = 0; i < size_; i++) {
                        cmMatch.add(i,i);
                    }
                    cmMatch.calculateBestFit(cmA, cmB);
                    float newRMSD = cmMatch.rmsd();
                    if(newRMSD < foundRMSD) {
                        // std::cout << "replacing chains in calcRMSD " << newRMSD << " " << foundRMSD << " " << identGroup[i] << " " << identGroup[j] << " " << bbIdToThisBBIndex[identGroup[i]] << " " << bbIdToThisBBIndex[identGroup[j]] << std::endl;
                        unsigned int tmp = bbIdToThisBBIndex[identGroup[i]];
                        bbIdToThisBBIndex[identGroup[i]] = bbIdToThisBBIndex[identGroup[j]];
                        bbIdToThisBBIndex[identGroup[j]] = tmp;
                        foundRMSD = newRMSD;
                        somethingChanged = true;
                    } else {
                        // std::cout << "not replacing chains in calcRMSD " << newRMSD << " " << foundRMSD << " " << identGroup[i] << " " << identGroup[j] << std::endl;
                    }
                }
            }
        }
        if(!somethingChanged)
            break;
        if(i == MAX_ITER - 1)
            std::cout << "stops clustering because max iteration " << bitIDS_ << std::endl;
    }

    if(foundRMSD > 1.5) return foundRMSD;

    Match match;
    Molecule<Atom> A, B;
    for(unsigned int bitId: bbIdsInThis) {
        Molecule<Atom> molA = bbs_[bbIdToThisBBIndex[bitId]]->caAtoms_;
        molA.rigidTrans(trans_[bbIdToThisBBIndex[bitId]]);
        A.concat(molA);

        Molecule<Atom> molB = other.bbs_[bbIdToOtherBBIndex[bitId]]->caAtoms_;
        molB.rigidTrans(other.trans_[bbIdToOtherBBIndex[bitId]]);
        B.concat(molB);
    }

    for(unsigned int i = 0; i < A.size(); i++) {
        match.add(i,i);
    }
    match.calculateBestFit(A, B);
    //std::cerr << A.size() << " " << B.size() << " rmsd " <<  match.rmsd() << std::endl;
    return match.rmsd();
}



double SuperBB::calcRmsd(const SuperBB &other) const {
  // calculate on CMs (centroids) first
  Match cmMatch;
  Molecule<Vector3> cmA, cmB;

  // create transformed cmA and cmB mols
  for(unsigned int i = 0; i < size_; i++) {
    Vector3 a = trans_[i]*bbs_[i]->cm_;
    cmA.add(a);

    // look for the same BB in other
    for(unsigned int j = 0; j < size_; j++) {
      if(bbs_[i]->getID() == other.bbs_[j]->getID()) { // same BB
        Vector3 b = other.trans_[j] * other.bbs_[j]->cm_;
        cmB.add(b);
        break;
      }
    }
  }

  for(unsigned int i = 0; i < size_; i++) {
    cmMatch.add(i,i);
  }
  cmMatch.calculateBestFit(cmA, cmB);
  if(cmMatch.rmsd() > 1.5) return cmMatch.rmsd();

  Match match;
  Molecule<Atom> A, B;
  // create transformed A and B mols
  for(unsigned int i = 0; i < size_; i++) {
    Molecule<Atom> mol = bbs_[i]->caAtoms_;
    mol.rigidTrans(trans_[i]);
    A.concat(mol);

    // look for the same BB in other
    for(unsigned int j = 0; j < size_; j++) {
      if(bbs_[i]->getID() == other.bbs_[j]->getID()) { // same BB
        Molecule<Atom> mol = other.bbs_[j]->caAtoms_;
        mol.rigidTrans(other.trans_[j]);
        B.concat(mol);
        break;
      }
    }
  }

  for(unsigned int i = 0; i < A.size(); i++) {
    match.add(i,i);
  }
  match.calculateBestFit(A, B);
  //std::cerr << A.size() << " " << B.size() << " rmsd " <<  match.rmsd() << std::endl;
  return match.rmsd();
}


void SuperBB::fullReport(std::ostream &s) {
  //HierarchalFold::timer_.stop();

  HierarchicalFold::countResults_++;
  RigidTrans3 tr;
  float compact2(0);
  if (yetAnother_ == 0)
    yetAnother_ = newScore_ * ((float) singlePen_ / (singlePen_ - multPen_));
  if (rmsd_ == 0)
    rmsd_ = calcRmsd(tr, compact2);
  //cerr << " compact2 " << compact2 << endl;
  s << "size_ " << size_ << " rmsd " << rmsd_ << " numOfBuried_ " << numOfBuried_ << " numOfHydBuried_ "
    << numOfHydBuried_ << " buried_diff " << ((float) numOfHydBuried_) / (numOfBuried_ - numOfHydBuried_)
    << " newScore_ " << newScore_ << " yetAnother_ " << yetAnother_ << " transScore_ " << transScore_ << " XX_ "
    << numOfXX_ << " compact2 " << compact2 << " multPen_ " << multPen_ << " singlePen_ " << singlePen_ << " diffPen "
    << singlePen_ - multPen_ << " backBonePen_ " << backBonePen_ << " maxPen_ " << maxPen_
    << " restraintsRatio_ " << restraintsRatio_ << " weightedTransScore " << weightedTransScore_;
  s << " contact_ ";
  for (std::vector<bool>::iterator i = contact_.begin(); i < contact_.end(); ++i)
    s << (*i ? '1' : '0');
  s << " [";
  int i = 0;
  for (auto it = bbs_.begin(); it != bbs_.end(); it++, i++) {
    s << (*it)->id_ << "(" << trans_[i] << ")";
    if (it + 1 != bbs_.end()) s << ",";
  }
  s << "]  " << tr;
  FoldStep::outputFoldSteps(s, foldSteps_);
  s << std::endl;

  //HierarchalFold::timer_.start();
}

std::ostream &operator<<(std::ostream &s, const SuperBB &sbb) {
  // scores
  s << "size_ " << sbb.size_ << " newScore_ " << sbb.newScore_ << " backBonePen_ " << sbb.backBonePen_
    << " restraintsRatio_ " << sbb.restraintsRatio_;

  // transformations
  s << " [";
  int i = 0;
  for (auto it = sbb.bbs_.begin(); it != sbb.bbs_.end(); it++, i++) {
    s << (*it)->getID() << "(" << sbb.trans_[i] << ")";
    if (it + 1 != sbb.bbs_.end()) s << ",";
  }
  s << "]  ";

  FoldStep::outputFoldSteps(s, sbb.foldSteps_);
  s << std::endl;
  return s;
}

void SuperBB::preClustering() {
    for (unsigned int i = 0; i < size_; i++) {
        const BB &bb1 = *bbs_[i];
        std::vector<RigidTrans3> tr;
        for (unsigned int j = 0; j < size_; j++) {
            tr.push_back((!trans_[j]) * trans_[i]);
        }
        for (unsigned int sp = 0; sp < bb1.ipSet_.size(); sp++) { // loop over ips
            for (unsigned int j = 0; j < size_; j++) {
                if (i == j) continue;
                const BB &bb2 = *bbs_[j];
                const InterestPoint &ip = bb1.ipSet_[sp];
                float dist = bb2.getDistFromSurface(tr[j] * ip.position_); // check penetration
                if (dist < 1) {
                    contact_[bb1.id_ * numberOfBBs_ + bb2.id_] = contact_[bb2.id_ * numberOfBBs_ + bb1.id_] = true;
                }
            }
        }
    }
}


TransIterator2::TransIterator2(const SuperBB &bb1, const SuperBB &bb2, int pbb1, int pbb2) :
        bb1_(bb1), bb2_(bb2), index1_(0), index2_(0) {


    std::shared_ptr<const BB> b1 = bb1_.bbs_[index1_];
    while (b1->getID() != pbb1) {
        b1 = bb1_.bbs_[++index1_];
    }
    std::shared_ptr<const BB> b2 = bb2_.bbs_[index2_];
    while (b2->getID() != pbb2) {
        b2 = bb2_.bbs_[++index2_];
    }

    mediatorTrans_ = !bb2_.trans_[index2_];
    trans_ = &(b1->getTransformations(b2->getID())); //BB::trans(*b1, *b2);
    it_ = trans_->begin();
    if (!isAtEnd()) {
      // check if bb1.size+bb2.size is in the range of the transformations validity
      //std::cerr << "TransIterator2: bb1 " << b1->getID() << " bb2 " << b2->getID() << " trans size " << trans_->size() << std::endl;
      int newSize = bb1_.size() + bb2_.size();
      const std::pair<int, int>& range = b1->getTransformsRange(b2->getID());
      //std::cerr << "newSize " << newSize << " range [" << range.first << " : " << range.second << "]" << std::endl;
      if( newSize >= range.first && newSize <= range.second)
       generateTransformation();
       else
       it_ = trans_->end();
      //generateTransformation();
    }
}
