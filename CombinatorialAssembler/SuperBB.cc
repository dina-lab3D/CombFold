#include "SuperBB.h"
#include "HierarchicalFold.h"


SuperBB::SuperBB(std::shared_ptr<const BB> bb)
    : backBonePen_(0), transScore_(0), weightedTransScore_(100) {
    bbs_.push_back(bb);
    Vector3 v(0, 0, 0);
    Matrix3 M(1);
    RigidTrans3 T(M, v);
    trans_.push_back(T);
    size_ = 1;
    bitIDS_ = bb->bitId();
}

float getWeightedTransScore(std::vector<FoldStep> steps, std::vector<std::shared_ptr<const BB>> bbs) {
    unsigned int totalCa = 0;
    float totalScore = 0;
    std::map<unsigned int, unsigned int> bbIdToCaLength;
    for (std::shared_ptr<const BB> bb : bbs) {
        bbIdToCaLength[bb->getID()] = bb->getNumOfAtoms();
    }

    for (FoldStep step : steps) {
        std::vector<unsigned int> right;
        std::vector<unsigned int> left;

        right.push_back(step.i_);
        left.push_back(step.j_);

        std::vector<FoldStep> notUsedSteps(steps);

        while (!notUsedSteps.empty()) {
            FoldStep step2 = notUsedSteps.back();
            notUsedSteps.pop_back();
            if (step2.i_ == step.i_ && step2.j_ == step.j_)
                continue;
            if (std::find(right.begin(), right.end(), step2.i_) != right.end() &&
                std::find(right.begin(), right.end(), step2.j_) == right.end()) {
                right.push_back(step2.j_);
            } else if (std::find(right.begin(), right.end(), step2.j_) != right.end() &&
                       std::find(right.begin(), right.end(), step2.i_) == right.end()) {
                right.push_back(step2.i_);
            } else if (std::find(left.begin(), left.end(), step2.i_) != left.end() &&
                       std::find(left.begin(), left.end(), step2.j_) == left.end()) {
                left.push_back(step2.j_);
            } else if (std::find(left.begin(), left.end(), step2.j_) != left.end() &&
                       std::find(left.begin(), left.end(), step2.i_) == left.end()) {
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

void SuperBB::join(const RigidTrans3 &trans, const SuperBB &other, int bbPen, FoldStep &step, float transScore) {
    // now joining is simple
    for (unsigned int j = 0; j < other.size_; j++) {
        bbs_.push_back(other.bbs_[j]);
        trans_.push_back(trans * other.trans_[j]);
    }
    size_ += other.size_;
    bitIDS_ |= other.bitIDS_;
    backBonePen_ = bbPen;

    transScore_ += other.transScore_ + transScore;

    if (!other.foldSteps_.empty()) {
        foldSteps_.insert(foldSteps_.end(), other.foldSteps_.begin(), other.foldSteps_.end());
    }
    foldSteps_.push_back(step);

    weightedTransScore_ = getWeightedTransScore(foldSteps_, bbs_);
}

void SuperBB::replaceIdentBB(unsigned long oldBBBitId, std::shared_ptr<const BB> bb) {
    if((bb->bitId() & bitIDS_) != 0)
        throw std::runtime_error("Tried to replace BB with BB that already exists in SuperBB");

    std::shared_ptr<const BB> oldBB;
    unsigned int bbIndex = -1;
    for (unsigned int i = 0; i < size_; i++) {
        if (bbs_[i]->bitId() == oldBBBitId) {
            oldBB = bbs_[i];
            bbIndex = i;
            break;
        }
    }
    if(bbIndex == -1)
        throw std::runtime_error("SuperBB::replaceIdentBB: oldBBBitId not found in SuperBB");

    bbs_[bbIndex] = bb;
    bitIDS_ -= oldBB->bitId();
    bitIDS_ += bb->bitId();

    for(FoldStep &step : foldSteps_) {
        if(step.i_ == oldBB->getID())
            step.i_ = bb->getID();
        if(step.j_ == oldBB->getID())
            step.j_ = bb->getID();
    }
}

bool SuperBB::isPenetrating(const RigidTrans3 &trans, const SuperBB &other, float threshold) const {
    for (unsigned int i = 0; i < size_; i++) {
        RigidTrans3 t = (!trans_[i]) * trans;
        for (unsigned int j = 0; j < other.size_; j++) {
            if (bbs_[i]->isPenetrating(t * other.trans_[j], *other.bbs_[j], threshold))
                return true;
        }
    }
    return false;
}

double SuperBB::calcRmsd(const SuperBB &other, std::vector<std::vector<unsigned int>> &identGroups) const {
    return calcRmsd(other, identGroups, false);
}

// RMSD between two SBBs (assuming same BBs in each SBB)
// double SuperBB::calcRmsd(const SuperBB &other, std::vector<std::vector<unsigned int>> &identGroups, bool verbose) const {
//     std::vector<std::vector<unsigned int>> presentIdentGroups;

//     for (std::vector<unsigned int> identGroup : identGroups) {
//         std::vector<unsigned int> possiblyRelevantIdentGroup;
//         for (unsigned int i : identGroup) {
//             if ((bitIDS_ & (1 << i)) != 0)
//                 possiblyRelevantIdentGroup.push_back(i);
//         }
//         if (possiblyRelevantIdentGroup.size() >= 2)
//             presentIdentGroups.push_back(possiblyRelevantIdentGroup);
//     }

//     std::map<unsigned int, unsigned int> bbIdToThisBBIndex;
//     std::map<unsigned int, unsigned int> bbIdToOtherBBIndex;
//     std::vector<unsigned int> bbIdsInThis;
//     for (unsigned int i = 0; i < bbs_.size(); i++) {
//         bbIdsInThis.push_back(bbs_[i]->getID());
//         bbIdToThisBBIndex[bbs_[i]->getID()] = i;
//         for (unsigned int j = 0; j < size_; j++) {
//             if (bbs_[i]->getID() == other.bbs_[j]->getID()) { // same BB
//                 bbIdToOtherBBIndex[bbs_[i]->getID()] = j;
//                 break;
//             }
//         }
//     }
//     if(verbose) {
//         std::cout << "start map:" << std::endl;
//         std::cout << "bbIdToThisBBIndex: ";
//         for (const auto &[a, b] : bbIdToThisBBIndex)
//             std::cout << a << "->" << bbs_[b]->getID() << ", ";
        
//         std::cout << std::endl << "bbIdToOtherBBIndex: ";
//         for (const auto &[a, b] : bbIdToOtherBBIndex)
//             std::cout << a << "->" << other.bbs_[b]->getID() << ", ";
//         std::cout << std::endl;
//     }   

//     // base rmsd find
//     Match cmMatch;
//     Molecule<Vector3> cmA, cmB;

//     for (unsigned int bitId : bbIdsInThis) {
//         Vector3 a = trans_[bbIdToThisBBIndex[bitId]] * bbs_[bbIdToThisBBIndex[bitId]]->cm_;
//         cmA.add(a);
//         Vector3 b = other.trans_[bbIdToOtherBBIndex[bitId]] * other.bbs_[bbIdToOtherBBIndex[bitId]]->cm_;
//         cmB.add(b);
//     }
//     for (unsigned int i = 0; i < size_; i++) {
//         cmMatch.add(i, i);
//     }
//     cmMatch.calculateBestFit(cmA, cmB);
//     float foundRMSD = cmMatch.rmsd();

//     int MAX_ITER = 10;
//     // try replacing chains with each other
//     for (int i = 0; i < MAX_ITER; i++) {

//         if (foundRMSD < 0.1)
//             break;

//         // std::cout << "iter " << i << std::endl;
//         bool somethingChanged = false;
//         for (std::vector<unsigned int> identGroup : presentIdentGroups) {
//             int best_found_i = -1;
//             int best_found_j = -1;
//             float best_found_rmsd = foundRMSD;
//             for (unsigned int i = 0; i < identGroup.size(); i++) {
//                 for (unsigned int j = i + 1; j < identGroup.size(); j++) {
//                     Match cmMatch;
//                     Molecule<Vector3> cmA, cmB;

//                     for (unsigned int bitId : bbIdsInThis) {
//                         Vector3 a;
//                         if (bitId == identGroup[i]) {
//                             unsigned int replacedBitId = identGroup[j];
//                             a = trans_[bbIdToThisBBIndex[replacedBitId]] * bbs_[bbIdToThisBBIndex[replacedBitId]]->cm_;
//                         } else if (bitId == identGroup[j]) {
//                             unsigned int replacedBitId = identGroup[i];
//                             a = trans_[bbIdToThisBBIndex[replacedBitId]] * bbs_[bbIdToThisBBIndex[replacedBitId]]->cm_;
//                         } else {
//                             a = trans_[bbIdToThisBBIndex[bitId]] * bbs_[bbIdToThisBBIndex[bitId]]->cm_;
//                         }
//                         cmA.add(a);
//                         Vector3 b =
//                             other.trans_[bbIdToOtherBBIndex[bitId]] * other.bbs_[bbIdToOtherBBIndex[bitId]]->cm_;
//                         cmB.add(b);
//                     }
//                     for (unsigned int i = 0; i < size_; i++) {
//                         cmMatch.add(i, i);
//                     }
//                     cmMatch.calculateBestFit(cmA, cmB);
//                     float newRMSD = cmMatch.rmsd();
//                     if (newRMSD <= foundRMSD) {
//                         // unsigned int tmp = bbIdToThisBBIndex[identGroup[i]];
//                         // bbIdToThisBBIndex[identGroup[i]] = bbIdToThisBBIndex[identGroup[j]];
//                         // bbIdToThisBBIndex[identGroup[j]] = tmp;
//                         // foundRMSD = newRMSD;
//                         // somethingChanged = true;


//                         // if(verbose) {
//                         //     std::cout << "replacing chains in calcRMSD " << newRMSD << " " << foundRMSD << " " <<
//                         //     identGroup[i] << " " << identGroup[j] << " " << bbIdToThisBBIndex[identGroup[i]] << " " <<
//                         //     bbIdToThisBBIndex[identGroup[j]] << std::endl;
//                         //     std::cout << "bbIdToThisBBIndex: ";
//                         //     for (const auto &[a, b] : bbIdToThisBBIndex)
//                         //         std::cout << a << "->" << bbs_[b]->getID() << ", ";
//                         //     std::cout << std::endl;
//                         // }
//                         if (newRMSD < best_found_rmsd) {
//                             best_found_rmsd = newRMSD;
//                             best_found_i = i;
//                             best_found_j = j;
//                         }
//                     } else if (verbose){
//                         std::cout << "not replacing chains in calcRMSD " << newRMSD << " " << foundRMSD << " " <<
//                         identGroup[i] << " " << identGroup[j] << std::endl;
//                     }
//                 }
//             }
//             if (best_found_i != -1) {
//                 unsigned int tmp = bbIdToThisBBIndex[identGroup[best_found_i]];
//                 bbIdToThisBBIndex[identGroup[best_found_i]] = bbIdToThisBBIndex[identGroup[best_found_j]];
//                 bbIdToThisBBIndex[identGroup[best_found_j]] = tmp;
//                 if(verbose) {
//                     std::cout << "replacing chains in calcRMSD " << best_found_rmsd << " " << foundRMSD << " " <<
//                     identGroup[best_found_i] << " " << identGroup[best_found_j] << " " << bbIdToThisBBIndex[identGroup[best_found_i]] << " " <<
//                     bbIdToThisBBIndex[identGroup[best_found_j]] << std::endl;
//                     std::cout << "bbIdToThisBBIndex: ";
//                     for (const auto &[a, b] : bbIdToThisBBIndex)
//                         std::cout << a << "->" << bbs_[b]->getID() << ", ";
//                     std::cout << std::endl;
//                 }
//                 foundRMSD = best_found_rmsd;
//                 somethingChanged = true;
//             }
//         }
//         if (!somethingChanged)
//             break;
        
//         if (i == MAX_ITER - 1)
//             std::cout << "stops clustering because max iteration " << bitIDS_ << std::endl;
//     }

//     if(verbose) {
//         std::cout << "rmsd: " << foundRMSD << ", final map:" << std::endl;
//         std::cout << "bbIdToThisBBIndex: ";
//         for (const auto &[a, b] : bbIdToThisBBIndex)
//             std::cout << a << "->" << bbs_[b]->getID() << ", ";
        
//         std::cout << std::endl << "bbIdToOtherBBIndex: ";
//         for (const auto &[a, b] : bbIdToOtherBBIndex)
//             std::cout << a << "->" << other.bbs_[b]->getID() << ", ";
//         std::cout << std::endl;
//     }

//     if (foundRMSD > 1.5)
//         return foundRMSD;

//     Match match;
//     Molecule<Atom> A, B;
//     for (unsigned int bitId : bbIdsInThis) {
//         Molecule<Atom> molA = bbs_[bbIdToThisBBIndex[bitId]]->caAtoms_;
//         molA.rigidTrans(trans_[bbIdToThisBBIndex[bitId]]);
//         A.concat(molA);

//         Molecule<Atom> molB = other.bbs_[bbIdToOtherBBIndex[bitId]]->caAtoms_;
//         molB.rigidTrans(other.trans_[bbIdToOtherBBIndex[bitId]]);
//         B.concat(molB);
//     }

//     for (unsigned int i = 0; i < A.size(); i++) {
//         match.add(i, i);
//     }
//     match.calculateBestFit(A, B);
//     return match.rmsd();
// }

double SuperBB::calcRmsdWithMap(const SuperBB &other, std::map<int, int> &bbIdMap) const {
    double centroidsRMSD = calcRmsdCentroidsWithMap(other, bbIdMap);
    if (centroidsRMSD > 1.5)
        return centroidsRMSD;
    Match match;
    Molecule<Atom> A, B;

    // create transformed A and B mols
    for (unsigned int i = 0; i < size_; i++) {
        Molecule<Atom> mol = bbs_[i]->caAtoms_;
        mol.rigidTrans(trans_[i]);
        A.concat(mol);

        // look for the same BB in other
        int otherBBId = (bbIdMap.count(bbs_[i]->getID()) > 0) ? bbIdMap[bbs_[i]->getID()] : bbs_[i]->getID();
        for (unsigned int j = 0; j < size_; j++) {
            if (otherBBId == other.bbs_[j]->getID()) { // same BB
                Molecule<Atom> mol = other.bbs_[j]->caAtoms_;
                mol.rigidTrans(other.trans_[j]);
                B.concat(mol);
                break;
            }
        }
    }

    for (unsigned int i = 0; i < A.size(); i++) {
        match.add(i, i);
    }
    match.calculateBestFit(A, B);
    return match.rmsd();
}


double SuperBB::calcRmsdCentroidsWithMap(const SuperBB &other, std::map<int, int> &bbIdMap) const {
    // calculate on centroids first
    Match cmMatch;
    Molecule<Vector3> cmA, cmB;

    // create transformed cmA and cmB mols
    for (unsigned int i = 0; i < size_; i++) {
        Vector3 a = trans_[i] * bbs_[i]->cm_;
        cmA.add(a);
        int otherBBId = (bbIdMap.count(bbs_[i]->getID()) > 0) ? bbIdMap[bbs_[i]->getID()] : bbs_[i]->getID();

        // look for the same BB in other
        for (unsigned int j = 0; j < size_; j++) {
            if (otherBBId == other.bbs_[j]->getID()) { // same BB
                Vector3 b = other.trans_[j] * other.bbs_[j]->cm_;
                cmB.add(b);
                break;
            }
        }
    }

    for (unsigned int i = 0; i < size_; i++) {
        cmMatch.add(i, i);
    }
    cmMatch.calculateBestFit(cmA, cmB);
    return cmMatch.rmsd();
}

double SuperBB::calcRmsd(const SuperBB &other, std::vector<std::vector<unsigned int>> &identGroups, bool verbose) const {
    // create a list of present ident groups
    std::vector<std::vector<unsigned int>> presentIdentGroups;
    std::set<unsigned int> foundBitIds;
    for (std::vector<unsigned int> identGroup : identGroups) {
        std::vector<unsigned int> possiblyRelevantIdentGroup;
        for (unsigned int i : identGroup) {
            if ((bitIDS_ & (1 << i)) != 0){
                possiblyRelevantIdentGroup.push_back(i);
                foundBitIds.insert(i);
            }
            
        }
        if (possiblyRelevantIdentGroup.size() >= 2)
            presentIdentGroups.push_back(possiblyRelevantIdentGroup);
    }

    if (presentIdentGroups.size() == 0) {
        std::map<int, int> emptyMap;
        return calcRmsdWithMap(other, emptyMap);
    }

    // add the remaining BBs to presentIdentGroups as a single element groups 
    for (unsigned int i = 0; i < bbs_.size(); i++) {
        if (foundBitIds.count(bbs_[i]->getID()) == 0) {
            std::vector<unsigned int> possiblyRelevantIdentGroup;
            possiblyRelevantIdentGroup.push_back(bbs_[i]->getID());
            presentIdentGroups.push_back(possiblyRelevantIdentGroup);
        }
    }

    // get the smallest ident group
    unsigned int smallestIdentGroupIndex = 0;
    // int otherBaselineBBId = -1;
    for (unsigned int i = 1; i < presentIdentGroups.size(); i++) {
        if (presentIdentGroups[i].size() < presentIdentGroups[smallestIdentGroupIndex].size()) {
            smallestIdentGroupIndex = i;
            // otherBaselineBBId = presentIdentGroups[i][0];
        }
    }

    // get all indexes of the BBs in the smallest ident group
    std::vector<std::vector<unsigned int>> identGroupToThisBBIndexes;
    std::vector<std::vector<unsigned int>> identGroupToOtherBBIndexes;
    for (unsigned int i = 0; i < presentIdentGroups.size(); i++) {
        std::vector<unsigned int> thisBBIndexes;
        std::vector<unsigned int> otherBBIndexes;
        for (unsigned int j = 0; j < presentIdentGroups[i].size(); j++) {
            unsigned int bitId = presentIdentGroups[i][j];
            for (unsigned int k = 0; k < size_; k++) {
                if (bitId == bbs_[k]->getID()) {
                    thisBBIndexes.push_back(k);
                    break;
                }
            }
            for (unsigned int k = 0; k < other.size_; k++) {
                if (bitId == other.bbs_[k]->getID()) {
                    otherBBIndexes.push_back(k);
                    break;
                }
            }
        }
        identGroupToThisBBIndexes.push_back(thisBBIndexes);
        identGroupToOtherBBIndexes.push_back(otherBBIndexes);
    }


    unsigned int otherBaselineBBIndex = identGroupToOtherBBIndexes[smallestIdentGroupIndex][0];
    double bestRmsd = 9999;
    std::map<int, int> bestBBImposeIndexMap;
    for (unsigned int i = 0; i < identGroupToThisBBIndexes[smallestIdentGroupIndex].size(); i++) {
        unsigned int thisBaselineBBIndex = identGroupToThisBBIndexes[smallestIdentGroupIndex][i];

        RigidTrans3 matchTrans = !other.trans_[otherBaselineBBIndex] * trans_[thisBaselineBBIndex];

        // apply matchTrans on each centorid in other BB
        // for each bb in this, find the closest centroid in other, remove centroids already matched
        std::map<int, int> BBImposeIndexMap;
        for (unsigned int j = 0; j < presentIdentGroups.size(); j++) {
            std::set<unsigned int> usedOtherBitIndexes;
            for (unsigned int thisIndex : identGroupToThisBBIndexes[j]) {
                double minDist = 9999;
                unsigned int minIndex = -1;
                Vector3 a = trans_[thisIndex] * bbs_[thisIndex]->cm_;

                for (unsigned int otherIndex : identGroupToOtherBBIndexes[j]) {
                    if (usedOtherBitIndexes.count(otherIndex) != 0)
                        continue;
                    Vector3 b;
                    b = matchTrans * other.trans_[otherIndex] * other.bbs_[otherIndex]->cm_;
                    double dist = (a - b).norm();
                    if (dist < minDist) {
                        minDist = dist;
                        minIndex = otherIndex;
                    }
                }
                BBImposeIndexMap[bbs_[thisIndex]->getID()] = other.bbs_[minIndex]->getID();
                usedOtherBitIndexes.insert(minIndex);
            }
        }

        // get rmsd using mapping and calcRmsdWithMap
        double foundRMSD = calcRmsdCentroidsWithMap(other, BBImposeIndexMap);
        if (verbose) {
            std::cout << "current best: " << bestRmsd << ", foundRMSD " << foundRMSD << std::endl;
            std::cout << "BBImposeIndexMap: " << std::endl;
            for (auto it = BBImposeIndexMap.begin(); it != BBImposeIndexMap.end(); it++) {
                std::cout << it->first << "->" << it->second << ", ";
            }
            std::cout << std::endl;
        }
        if (foundRMSD < bestRmsd) {
            bestRmsd = foundRMSD;
            bestBBImposeIndexMap = BBImposeIndexMap;
        }
        if (foundRMSD < 0.1)
            break;
    }
    return calcRmsdWithMap(other, bestBBImposeIndexMap);
}

double SuperBB::calcRmsd(const SuperBB &other) const {
    std::vector<std::vector<unsigned int>> emptyIdentGroups;
    return calcRmsd(other, emptyIdentGroups);
}

void SuperBB::fullReport(std::ostream &s) {
    HierarchicalFold::countResults_++;
    RigidTrans3 tr;
    s << "size_ " << size_ << " transScore_ " << transScore_ << " multPen_ " << multPen_ << " singlePen_ " 
      << singlePen_ << " diffPen " << singlePen_ - multPen_ << " backBonePen_ " << backBonePen_ << " maxPen_ " 
      << maxPen_ << " restraintsRatio_ " << restraintsRatio_ << " weightedTransScore " << weightedTransScore_;
    s << " [";
    int i = 0;
    for (auto it = bbs_.begin(); it != bbs_.end(); it++, i++) {
        s << (*it)->id_ << "(" << trans_[i] << ")";
        if (it + 1 != bbs_.end())
            s << ",";
    }
    s << "]  " << tr;
    FoldStep::outputFoldSteps(s, foldSteps_);
    s << std::endl;
}

std::ostream &operator<<(std::ostream &s, const SuperBB &sbb) {
    // scores
    s << "size_ " << sbb.size_ << " backBonePen_ " << sbb.backBonePen_ << " restraintsRatio_ " << sbb.restraintsRatio_;

    // transformations
    s << " [";
    int i = 0;
    for (auto it = sbb.bbs_.begin(); it != sbb.bbs_.end(); it++, i++) {
        s << (*it)->getID() << "(" << sbb.trans_[i] << ")";
        if (it + 1 != sbb.bbs_.end())
            s << ",";
    }
    s << "]  ";

    FoldStep::outputFoldSteps(s, sbb.foldSteps_);
    s << std::endl;
    return s;
}

TransIterator2::TransIterator2(const SuperBB &bb1, const SuperBB &bb2, int pbb1, int pbb2)
    : bb1_(bb1), bb2_(bb2), index1_(0), index2_(0) {

    std::shared_ptr<const BB> b1 = bb1_.bbs_[index1_];
    while (b1->getID() != pbb1) {
        b1 = bb1_.bbs_[++index1_];
    }
    std::shared_ptr<const BB> b2 = bb2_.bbs_[index2_];
    while (b2->getID() != pbb2) {
        b2 = bb2_.bbs_[++index2_];
    }

    mediatorTrans_ = !bb2_.trans_[index2_];
    trans_ = &(b1->getTransformations(b2->getID())); // BB::trans(*b1, *b2);
    it_ = trans_->begin();
    if (!isAtEnd()) {
        generateTransformation();
    }
}
