#include "ComplexDistanceConstraint.h"

#include "BB.h"

void ComplexDistanceConstraint::addConstraint(int suInd1, int suInd2, Vector3 receptorAtom, Vector3 ligandAtom,
                                              float maxDistance, float minDistance) {
    int ind = suInd1 * noOfSUs_ + suInd2;
    constraints_[ind].push_back(DistanceRestraint(receptorAtom, ligandAtom, maxDistance, minDistance));
}

bool ComplexDistanceConstraint::areConstraintsSatisfied(int suInd1, int suInd2, RigidTrans3 &trans2) const {
    int ind = suInd1 * noOfSUs_ + suInd2;
    for (unsigned int index = 0; index < constraints_[ind].size(); index++) {
        bool sat = constraints_[ind][index].isSatisfied(trans2);
        if (!sat) {
            return false;
        }
    }
    return true;
}

int ComplexDistanceConstraint::addChainConnectivityConstraints() {
    for (int suInd1 = 0; suInd1 < noOfSUs_; suInd1++) {
        for (int suInd2 = suInd1 + 1; suInd2 < noOfSUs_; suInd2++) {
            std::vector<std::pair<char, std::pair<int, int>>> constraints;
            int counter = 0;
            float epsilon = 20.0; // for shorter connections
            bbs_[suInd1]->getChainConnectivityConstraints(*bbs_[suInd2], constraints);
            for (int i = 0; i < (int)constraints.size(); i++) {
                int su1resIndex = constraints[i].second.first;
                int su2resIndex = constraints[i].second.second;

                int receptorAtomIndex = bbs_[suInd1]->allAtoms_.getFirstAtomEntryForResIndex(
                    constraints[i].first, std::to_string(su1resIndex));
                int ligandAtomIndex = bbs_[suInd2]->allAtoms_.getFirstAtomEntryForResIndex(constraints[i].first,
                                                                                           std::to_string(su2resIndex));

                Vector3 receptorAtom = bbs_[suInd1]->getChemAtomByIndex(receptorAtomIndex).position();
                Vector3 ligandAtom = bbs_[suInd2]->getChemAtomByIndex(ligandAtomIndex).position();

                int sequenceDist = std::fabs(su1resIndex - su2resIndex);
                if (sequenceDist < 100) {
                    float maxDist = sequenceDist * 2.0;
                    if (sequenceDist <= 5)
                        maxDist += epsilon;
                    addConstraint(suInd1, suInd2, receptorAtom, ligandAtom, maxDist);
                    addConstraint(suInd2, suInd1, ligandAtom, receptorAtom, maxDist);
                    std::cerr << "DISTANCE CONSTRAINT ADDED: " << constraints[i].first << " su1 endpoint "
                              << su1resIndex << " su2 endpoint " << su2resIndex << " " << maxDist << std::endl;
                    counter++;
                }
            }
            numberOfConstraints_ += counter;
        }
    }
    return numberOfConstraints_;
}

std::set<int> ComplexDistanceConstraint::getSUs(const std::string residueSequenceID, const std::string chains,
                                                int maxOffset) const {

    int rOffset = 0;
    std::set<int> ret;
    // iterate chains
    for (unsigned int chainIndex = 0; chainIndex < chains.size(); chainIndex++) {
        // for each chain find corresponding SU
        for (int suIndex = 0; suIndex < noOfSUs_; suIndex++) {
            int atomIndex = bbs_[suIndex]->allAtoms_.getClosestAtomEntryForResIndex(
                chains[chainIndex], residueSequenceID, maxOffset, rOffset);
            if (atomIndex != -1) { // found
                ret.insert(suIndex);
                continue;
            }
        }
    }
    return ret;
}

int ComplexDistanceConstraint::readRestraintsFile(const std::string fileName) {
    // read cross links
    std::vector<CrossLink> crosslinks;
    int xnum = readCrossLinkFile(fileName, crosslinks);
    Logger::infoMessage() << "# of xlinks " << xnum << " read from file" << fileName << std::endl;
    int MAX_OFFSET = 20;

    crosslinkIndToWeight_ = std::vector<float>(crosslinks.size());

    // map them to BB pairs
    for (unsigned int i = 0; i < crosslinks.size(); i++) {
        crosslinkIndToWeight_[i] = crosslinks[i].getWeight();

        // SUs that have xlinks endpoints (sus1, sus2)
        std::string residueSequenceID1 = std::to_string(crosslinks[i].getResidue1());
        std::set<int> sus1 = getSUs(residueSequenceID1, crosslinks[i].getChain1());
        if (sus1.size() == 0)
            sus1 = getSUs(residueSequenceID1, crosslinks[i].getChain1(), MAX_OFFSET);
        std::string residueSequenceID2 = std::to_string(crosslinks[i].getResidue2());
        std::set<int> sus2 = getSUs(residueSequenceID2, crosslinks[i].getChain2());
        if (sus2.size() == 0)
            sus2 = getSUs(residueSequenceID2, crosslinks[i].getChain2(), MAX_OFFSET);

        // iterate over all pairs of SUs
        for (auto suIndexIter1 = sus1.begin(); suIndexIter1 != sus1.end(); suIndexIter1++) {
            for (unsigned int chainIndex1 = 0; chainIndex1 < crosslinks[i].getChain1().size(); chainIndex1++) {
                char chainId1 = crosslinks[i].getChain1()[chainIndex1];
                int rOffset = 0;
                int receptorAtomIndex = bbs_[*suIndexIter1]->allAtoms_.getClosestAtomEntryForResIndex(
                    chainId1, residueSequenceID1, MAX_OFFSET, rOffset);
                if (rOffset < 0)
                    rOffset = -1 * rOffset;
                if (receptorAtomIndex == -1)
                    continue;
                Vector3 rcoord = bbs_[*suIndexIter1]->getChemAtomByIndex(receptorAtomIndex).position();

                for (auto suIndexIter2 = sus2.begin(); suIndexIter2 != sus2.end(); suIndexIter2++) {
                    for (unsigned int chainIndex2 = 0; chainIndex2 < crosslinks[i].getChain2().size(); chainIndex2++) {
                        char chainId2 = crosslinks[i].getChain2()[chainIndex2];
                        int lOffset = 0;
                        int ligandAtomIndex = bbs_[*suIndexIter2]->allAtoms_.getClosestAtomEntryForResIndex(
                            chainId2, residueSequenceID2, MAX_OFFSET, lOffset);
                        if (ligandAtomIndex == -1)
                            continue;
                        if (lOffset < 0)
                            lOffset = -1 * lOffset;

                        Vector3 lcoord = bbs_[*suIndexIter2]->getChemAtomByIndex(ligandAtomIndex).position();

                        float maxDistance = crosslinks[i].getMaxDistance() + (rOffset + lOffset) * 3.0;
                        float minDistance = crosslinks[i].getMinDistance();

                        int ind1 = *suIndexIter1 * noOfSUs_ + *suIndexIter2;
                        DistanceRestraint d1(rcoord, lcoord, maxDistance, minDistance, crosslinks[i].getWeight());
                        restraintIndsToCrosslinkInds_[ind1].push_back(i);
                        restraints_[ind1].push_back(d1);

                        int ind2 = *suIndexIter2 * noOfSUs_ + *suIndexIter1;
                        DistanceRestraint d2(lcoord, rcoord, maxDistance, minDistance, crosslinks[i].getWeight());
                        restraintIndsToCrosslinkInds_[ind2].push_back(i);
                        restraints_[ind2].push_back(d2);

                        std::cout << " adding restraint to " << *suIndexIter1 << "x" << *suIndexIter2 << " :"
                                  << residueSequenceID1 << chainId1 << " : " << residueSequenceID2 << chainId2
                                  << " dist " << maxDistance << " " << minDistance << " indexes " << ind1 << " : "
                                  << ind2 << std::endl;
                    }
                }
            }
        }
    }
    return crosslinks.size();
}

float ComplexDistanceConstraint::getRestraintsRatio(const std::vector<std::shared_ptr<const BB>> &bbs,
                                                    const std::vector<RigidTrans3> &trans) const {
    std::set<unsigned int> totalSeenCrosslinks;
    std::set<unsigned int> satisfiedCrosslinks;
    for (int suInd1 = 0; suInd1 < (int)bbs.size(); suInd1++) {
        for (int suInd2 = 0; suInd2 < (int)bbs.size(); suInd2++) {
            int index1 = bbs[suInd1]->getID();
            int index2 = bbs[suInd2]->getID();
            if (index1 == index2)
                continue;
            int suPairIndex = index1 * noOfSUs_ + index2;
            for (unsigned int restraintIndex = 0; restraintIndex < restraints_[suPairIndex].size(); restraintIndex++) {
                unsigned int crosslinkIndex = restraintIndsToCrosslinkInds_[suPairIndex][restraintIndex];
                if (satisfiedCrosslinks.find(crosslinkIndex) != satisfiedCrosslinks.end()) {
                    continue;
                }
                totalSeenCrosslinks.insert(crosslinkIndex);
                if (restraints_[suPairIndex][restraintIndex].isSatisfied(trans[suInd1], trans[suInd2]))
                    satisfiedCrosslinks.insert(crosslinkIndex);
            }
        }
    }
    if (totalSeenCrosslinks.size() == 0)
        return 1.0;

    float satisfiedWeight = 0.0;
    float totalWeight = 0.0;
    for (auto it = satisfiedCrosslinks.begin(); it != satisfiedCrosslinks.end(); it++) {
        satisfiedWeight += crosslinkIndToWeight_[*it];
    }
    for (auto it = totalSeenCrosslinks.begin(); it != totalSeenCrosslinks.end(); it++) {
        totalWeight += crosslinkIndToWeight_[*it];
    }

    return satisfiedWeight / totalWeight;
    // return (float)satisfiedCrosslinks.size() / (float)totalSeenCrosslinks.size();
}