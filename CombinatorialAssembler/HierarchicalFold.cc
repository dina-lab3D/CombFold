#include "HierarchicalFold.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

Timer HierarchicalFold::timer_;
Timer HierarchicalFold::timerAll_;
unsigned int HierarchicalFold::countResults_(0);

std::vector<std::vector<unsigned int>> createIdentGroups(unsigned int N_, BestKContainer &bestKContainer_) {
    std::vector<bool> addedToGroup(N_, false);
    std::vector<std::vector<unsigned int>> identGroups;
    for (unsigned int i = 0; i < N_; i++) {
        std::shared_ptr<SuperBB> sbbI = *(bestKContainer_[BitId(i)].begin());
        const std::shared_ptr<const BB> bbI = sbbI->bbs_[0];

        if (addedToGroup[i])
            continue;
        addedToGroup[i] = true;
        std::vector<unsigned int> identical;

        for (unsigned int j = i + 1; j < N_; j++) {
            std::cout << "checking ident " << i << " & " << j << std::endl;
            std::shared_ptr<SuperBB> sbbJ = *(bestKContainer_[BitId(j)].begin());

            const std::shared_ptr<const BB> bbJ = sbbJ->bbs_[0];

            if (bbI->isIdent(*bbJ)) {
                addedToGroup[j] = true;
                std::cout << "found ident " << i << " " << j << std::endl;
                identical.push_back(j);
            }
        }
        if (identical.size() > 0) {
            identical.insert(identical.begin(), i);
            identGroups.push_back(identical);
        }
    }
    std::cout << "---- ident groups " << std::endl;
    for (std::vector<unsigned int> identGroup : identGroups) {
        std::cout << "ident group: ";
        for (unsigned int i : identGroup) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }

    return identGroups;
}

std::map<unsigned int, std::vector<unsigned int>> createAssemblyGroupsMap(unsigned int N_,
                                                                          BestKContainer &bestKContainer_) {
    std::map<unsigned int, std::vector<unsigned int>> assemblyGroupsMap;
    for (unsigned int i = 0; i < N_; i++) {
        std::shared_ptr<SuperBB> sbbI = *(bestKContainer_[BitId(i)].begin());
        const std::shared_ptr<const BB> bbI = sbbI->bbs_[0];
        if (assemblyGroupsMap.count(bbI->groupId()) == 0) {
            std::vector<unsigned int> newGroup;
            assemblyGroupsMap[bbI->groupId()] = newGroup;
        }
        assemblyGroupsMap[bbI->groupId()].push_back(i);
    }
    std::cout << "---- assembly groups " << std::endl;
    // std::vector<std::vector<unsigned int>> assemblyGroups;
    for (const auto &[groupId, groupBBIds] : assemblyGroupsMap) {
        std::cout << "assembly group " << groupId << ":";
        for (unsigned int i : groupBBIds) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
        // assemblyGroups.push_back(groupBBIds);
    }
    return assemblyGroupsMap;
}

void printBestK(unsigned int N_, BestK *bestK) {
    // just for print
    std::map<unsigned int, unsigned int> count_new_kept_by_bb;
    std::unordered_map<BitId, unsigned int> count_new_kept_by_resSet;
    for (unsigned int bit_index = 0; bit_index < N_; bit_index++)
        count_new_kept_by_bb[bit_index] = 0;
    for (auto it1 = bestK->begin(); it1 != bestK->end(); it1++) {
        if (count_new_kept_by_resSet.count((**it1).bitIds()) == 0)
            count_new_kept_by_resSet[(**it1).bitIds()] = 0;
        count_new_kept_by_resSet[(**it1).bitIds()] += 1;
        for (unsigned int bit_index = 0; bit_index < N_; bit_index++)
            if ((**it1).bitIds().test(bit_index)) {
                count_new_kept_by_bb[bit_index] += 1;
            }
    }

    std::cout << "scores of saved " << bestK->minScore() << ":" << bestK->maxScore() << std::endl;
    std::cout << "new kept results by chain ";
    for (const auto &elem : count_new_kept_by_bb)
        std::cout << elem.first << ":" << elem.second << ", ";
    std::cout << std::endl << "new kept results by resSet ";
    for (const auto &elem : count_new_kept_by_resSet)
        std::cout << elem.first << ":" << elem.second << ", ";
    std::cout << std::endl;
}

bool isValidBasedOnAssembly(std::map<unsigned int, std::vector<unsigned int>> assemblyGroupsMap,
                            BitId currResSet) {
    if (assemblyGroupsMap.size() <= 1)
        return true;

    bool seenPartialGroup = false;
    int numberOfGroups = 0;

    bool seenGroupZero = false;

    for (const auto &[groupId, assemblyGroup] : assemblyGroupsMap) {
        bool seenZero = false;
        bool seenOne = false;

        for (unsigned int i = 0; i < assemblyGroup.size(); i++) {
            if (currResSet.test(assemblyGroup[i]))
                seenOne = true;
            else
                seenZero = true;
        }
        if (groupId == 0) {
            // assembly group 0 is special, this group should not be assembled before joining others
            // so we never treat it as partial, but it can't be with another partial group
            if (seenOne)
                seenGroupZero = true;
        } else {
            if (seenOne)
                numberOfGroups++;

            if (seenZero && seenOne) {
                seenPartialGroup = true;
            }
        }
    }
    if (numberOfGroups >= 1 && seenPartialGroup && seenGroupZero)
        return false;


    if (numberOfGroups > 1 && seenPartialGroup)
        return false;

    return true;
}

void HierarchicalFold::fold(const std::string &outFileNamePrefix) {
    std::vector<std::vector<unsigned int>> identGroups = createIdentGroups(N_, bestKContainer_);
    std::map<unsigned int, std::vector<unsigned int>> assemblyGroupsMap = createAssemblyGroupsMap(N_, bestKContainer_);

    std::map<unsigned int, BestK *> precomputedResults;
    for (unsigned int i = 2; i <= N_; i++)
        precomputedResults[i] = new BestK(K_);

    // populate with homomers subunits
    for (std::vector<unsigned int> identGroup : identGroups) {
        for (unsigned int groupDivider = 1; (identGroup.size() / groupDivider) >= 5; groupDivider++) {
        if ((identGroup.size() % groupDivider) != 0)
            continue;
        unsigned int groupSize = identGroup.size() / groupDivider;
        std::vector<std::shared_ptr<SuperBB>> groupSBBs;
        for (unsigned int j = 0; j < groupSize; j++) {
            groupSBBs.push_back(*(bestKContainer_[BitId(identGroup[j])].begin()));
        }
        std::cout << "searching for size " << groupSBBs.size() << " has " << precomputedResults.count(groupSBBs.size())
                  << std::endl;
        createSymmetry(groupSBBs, *precomputedResults[groupSBBs.size()]);
        }
    }

    // Hierarchical Assembly
    for (unsigned int length = 2; length <= N_; length++) { // # subunits iteration
        std::cout << "*** running iteration " << length
                  << " prev kept results: " << keptResultsByLength[length - 1]->size() << std::endl;
        std::unordered_map<BitId, BestK *> best_k_by_id;

        // populate with precomputedResults
        for (auto it1 = precomputedResults[length]->begin(); it1 != precomputedResults[length]->end(); it1++) {
            BitId currResSet = (**it1).bitIds();
            if (best_k_by_id.count(currResSet) == 0)
                best_k_by_id[currResSet] = new BestK(K_);
            best_k_by_id[currResSet]->push(*it1);
        }

        // try to assemble from each pair of kept results that together have (length) subunits
        for (unsigned int firstResultSize = 1; firstResultSize <= length / 2; firstResultSize++) {
            unsigned int secondResultSize = length - firstResultSize;
            std::cout << "** running sub-iteration " << firstResultSize << " " << secondResultSize << std::endl;
            std::cout << "counters " << countFilterTrasSkipped_ << "/" << countFilterTras_ << std::endl;

            for (auto it1 = keptResultsByLength[firstResultSize]->begin();
                 it1 != keptResultsByLength[firstResultSize]->end(); it1++) {
                SuperBB sbb1 = **it1;
                BitId setA = sbb1.bitIds();

                for (auto it2 = keptResultsByLength[secondResultSize]->begin();
                     it2 != keptResultsByLength[secondResultSize]->end(); it2++) {

                    if (firstResultSize == secondResultSize &&
                        std::distance(keptResultsByLength[firstResultSize]->begin(), it1) >
                            std::distance(keptResultsByLength[secondResultSize]->begin(), it2))
                        continue; // Since in this case there are 2 identical loops, don't do things twice

                    // If there are identical subunits in both results, rewrite the second result to not have the same
                    std::shared_ptr<SuperBB> sbb2Pointer = getMatchingSBB(sbb1, **it2, identGroups);
                    if (sbb2Pointer == NULL)
                        continue;
                    SuperBB sbb2 = *sbb2Pointer;

                    // make sure that the two results can be connected
                    BitId setB = sbb2.bitIds();
                    if ((setA & setB) != 0)
                        continue;
                    BitId currResSet = setA | setB;
                    if(!isValidBasedOnAssembly(assemblyGroupsMap, currResSet)){
                        std::cout << "invalid assembly " << currResSet << std::endl;
                        continue;
                    }

                    // connect the two results and add all new combined results to best_k_by_id[currResSet]
                    if (best_k_by_id.count(currResSet) == 0)
                        best_k_by_id[currResSet] = new BestK(K_);

                    unsigned int resCountBefore = best_k_by_id[currResSet]->size();
                    float minScoreBefore = best_k_by_id[currResSet]->minScore();

                    std::promise<int> promise1;
                    this->tryToConnect(1, sbb1, sbb2, *best_k_by_id[currResSet], (length < N_), promise1, identGroups);

                    if (resCountBefore < best_k_by_id[currResSet]->size() ||
                        minScoreBefore != best_k_by_id[currResSet]->minScore())
                        std::cout << "found more for " << currResSet << " based on " << setA << " and " << setB
                                  << " before: " << resCountBefore << " after: " << best_k_by_id[currResSet]->size()
                                  << " scores " << best_k_by_id[currResSet]->minScore() << ":"
                                  << best_k_by_id[currResSet]->maxScore() << std::endl;
                }
            }
        }

        // cluster results and save them
        std::map<unsigned int, BestK *> bestForSubunitId;
        keptResultsByLength[length] = new BestK(K_);

        for (const auto &[currResSet, currBestK] : best_k_by_id) {
            if (currBestK->size() > 0) {
                BestK *clusteredBestK = bestKContainer_.newBestK(currResSet);
                currBestK->cluster(*clusteredBestK, 1.0, identGroups);

                std::cerr << "clustering resSet " << currResSet << " before: " << currBestK->size() << " after "
                          << bestKContainer_[currResSet].size() << " scores " << bestKContainer_[currResSet].minScore()
                          << ":" << bestKContainer_[currResSet].maxScore() << std::endl;

                for (unsigned int i = 0; i < N_; i++) {
                    if (currResSet.test(i)) {
                        if (bestForSubunitId.count(i) == 0) {
                            bestForSubunitId[i] = new BestK(1);
                        }
                        bestForSubunitId[i]->push(*clusteredBestK->rbegin());
                    }
                }

                unsigned int count = 0;
                for (auto it1 = bestKContainer_[currResSet].rbegin(); it1 != bestKContainer_[currResSet].rend(); it1++) {
                    keptResultsByLength[length]->push(*it1);
                    count += 1;
                    if (count >= maxResultPerResSet)
                        break;
                }
            }
            delete currBestK;
        }

        // save best from each subunit
        keptResultsByLength[length]->setK(K_ + bestForSubunitId.size());
        for (const auto &[subunitId, currBestK] : bestForSubunitId) {
            keptResultsByLength[length]->push(*currBestK->rbegin());
            delete currBestK;
        }

        printBestK(N_, keptResultsByLength[length]);
    }

    // output fully assembled results or largest subsets
    if (keptResultsByLength[N_]->size() != 0) {
        std::string outFileName = outFileNamePrefix + ".res";
        std::ofstream outFile(outFileName);
        std::ofstream outFileClustered(outFileNamePrefix + "_clustered.res");
        BestK clusteredBestK(finalSizeLimit_); // TODO: this should also change on the best_k_by_id level

        for (auto it = keptResultsByLength[N_]->rbegin(); it != keptResultsByLength[N_]->rend(); it++)
            (*it)->fullReport(outFile);
        // output after clustering
        keptResultsByLength[N_]->cluster(clusteredBestK, 5.0, identGroups);
        for (auto it = clusteredBestK.rbegin(); it != clusteredBestK.rend(); it++)
            (*it)->fullReport(outFileClustered);
        outFile.close();
        outFileClustered.close();
    } else {
        // output largest subsets
        for (unsigned int i = N_; i > 1; i--) {
            if (keptResultsByLength[i]->size() == 0)
                continue;
            std::string outFileName = "cb_" + std::to_string(i) + "_" + outFileNamePrefix + ".res";
            std::ofstream outFile(outFileName);
            for (auto it = keptResultsByLength[i]->rbegin(); it != keptResultsByLength[i]->rend(); it++)
                (*it)->fullReport(outFile);
            outFile.close();
            break;
        }
    }

    // cleanup
    for (const auto &[length, currBestK] : precomputedResults) {
        delete currBestK;
    }
    for (const auto &[length, currBestK] : keptResultsByLength) {
        delete currBestK;
    }
}

void HierarchicalFold::tryToConnect(int id, const SuperBB &sbb1, const SuperBB &sbb2, BestK &results, bool toAdd,
                                    std::promise<int> &output, std::vector<std::vector<unsigned int>> &identGroups) {
    // iterate over pairs of BBs os SuperBB1 and SuperBB2
    for (int i = 0; i < (int)sbb1.bbs_.size(); i++) {
        int firstBB = sbb1.bbs_[i]->getID();
        for (int j = 0; j < (int)sbb2.bbs_.size(); j++) {
            int secondBB = sbb2.bbs_[j]->getID();

            // loop over possible transformations between BBs
            for (TransIterator2 it(sbb1, sbb2, firstBB, secondBB); !it.isAtEnd(); it++) {
                // optimization - check that the score is not lower than the minimum in the current bestK
                if((it.getScore() + sbb1.transScore_ + sbb2.transScore_) < results.minScore()){
                    continue;
                }

                // discard any invalid transformations
                bool filtered = filterTrans(sbb1, sbb2, it.transformation());
                if (filtered)
                    continue;

                FoldStep step(firstBB, secondBB, it.getScore());
                std::shared_ptr<SuperBB> theNew = createJoined(sbb1, sbb2, it.transformation(), 0, step, it.getScore());

                if (theNew->getRestraintsRatio() < restraintsRatioThreshold_) {
                    //            std::cout << "not enough restraints " << theNew->getRestraintsRatio() << " : " <<
                    //            complexConst_.getDistanceRestraintsRatioThreshold();
                    continue;
                }

                // results.push(theNew);
                results.push_cluster(theNew, 1, identGroups);
            }
        }
    }
    output.set_value(1);
}

bool HierarchicalFold::filterTrans(const SuperBB &sbb1, const SuperBB &sbb2, const RigidTrans3 &trans) const {

    // check distance constraints & restraints
    for (unsigned int i = 0; i < sbb1.size_; i++) {
        const BB &bb1 = *sbb1.bbs_[i];
        RigidTrans3 t = (!sbb1.trans_[i]) * trans;
        for (unsigned int j = 0; j < sbb2.size_; j++) {
            const BB &bb2 = *sbb2.bbs_[j];
            RigidTrans3 t2 = t * sbb2.trans_[j];
            // check constraints first
            if (!complexConst_.areConstraintsSatisfied(bb1.getID(), bb2.getID(), t2))
                return true;
        }
    }

    // backbone penetrations for each pair of BBs
    for (unsigned int i = 0; i < sbb1.size_; i++) {
        const BB &bb1 = *sbb1.bbs_[i];
        RigidTrans3 t = (!sbb1.trans_[i]) * trans;

        for (unsigned int j = 0; j < sbb2.size_; j++) {
            const BB &bb2 = *sbb2.bbs_[j];
            RigidTrans3 t2 = t * sbb2.trans_[j];

            const BB *pBB1 = &bb1, *pBB2 = &bb2;
            if (bb2.getSurfaceSize() > bb1.getSurfaceSize()) {
                pBB1 = &bb2;
                pBB2 = &bb1;
                t2 = !t2;
            }

            countFilterTras_ = countFilterTras_ + 1;

            // optimization - check if radiuses are too far apart and if so, skip check
            if ((pBB1->getRadius() + pBB2->getRadius()) < (pBB1->getCM() - t2*pBB2->getCM()).norm()) {
                countFilterTrasSkipped_ = countFilterTrasSkipped_ + 1;
                continue;
            }

            unsigned int bbPenetrations = 0;
            unsigned int totalUsedAtoms = 0;

            // TODO: maybe should save Weighted bbPen using pBB1->grid_->getDist(v) as weight
            for (Molecule<Atom>::const_iterator it = pBB2->caAtoms_.begin(); it != pBB2->caAtoms_.end(); it++) {
                if (it->getTempFactor() < minTemperatureToConsiderCollision) {
                    continue;
                }
                totalUsedAtoms++;

                Vector3 v = t2 * it->position();
                if (pBB1->getDistFromSurface(v) < 0) {
                    // getResidueEntry(v) when used in BBGrid.h will return -1*res_index if res_index is backbone
                    if (pBB1->grid_->getResidueEntry(v) < 0 && pBB1->grid_->getDist(v) < penetrationThreshold_) {
                        int resEntry = pBB1->grid_->getResidueEntry(v) * -1;
                        if (pBB1->getAtomByResId(resEntry).getTempFactor() < minTemperatureToConsiderCollision)
                            continue;

                        bbPenetrations++;
                    }
                }
            }
            float bbPenChangePercent = (float)(bbPenetrations) / (float)totalUsedAtoms;
            if (bbPenChangePercent > maxBackboneCollisionPercentPerChain) {
                return true;
            }
        }
    }

    return false;
}

void HierarchicalFold::createSymmetry(std::vector<std::shared_ptr<SuperBB>> identBBs, BestK &results) {
    unsigned int transCount = 0;
    std::cout << "started trans check, bb_size:" << identBBs.size() << std::endl;
    for (TransIterator2 it(*identBBs[0], *identBBs[1], identBBs[0]->bbs_[0]->getID(), identBBs[1]->bbs_[0]->getID());
         !it.isAtEnd(); it++)
        transCount++;

    std::cout << "number of transformations:" << transCount << std::endl;
    unsigned int addedSymCount = 0;
    for (unsigned int transNum = 0; transNum < transCount; transNum++) {
        std::cout << "checking trans indexed" << transNum << std::endl;
        // create symSBB for a trans, this is cumbersome because I don't really know how to handle transformations
        std::shared_ptr<SuperBB> symSBB = identBBs[0];
        for (unsigned int i = 1; i < identBBs.size(); i++) {
            unsigned int count = -1;

            // There is a memory issue here, I create TransIterator2 with symSBB as bb1, but then I override it
            // and freeing the memory(?) of bb1_, so it is important to break after chanigng symSBB
            for (TransIterator2 it(*symSBB, *identBBs[i], symSBB->bbs_[i - 1]->getID(), identBBs[i]->bbs_[0]->getID());
                 !it.isAtEnd(); it++) {
                count++;
                if (count != transNum)
                    continue;
                float transScore = it.getScore() + it.getScore() * ((100 - it.getScore()) / 100);

                FoldStep step(symSBB->bbs_[i - 1]->getID(), identBBs[i]->bbs_[0]->getID(), transScore);
                std::cout << "adding trans " << it.transformation() << " **** " << it.getScore() << std::endl;
                symSBB = createJoined(*symSBB, *identBBs[i], it.transformation(), 0, step, transScore);
                break;
            }
        }
        std::cout << "created possibly symSBB" << std::endl;

        // check bb penetration between each 2 chains
        double maxPenetration = 0;
        bool shouldContinuePen = false;
        for (unsigned int i = 0; i < symSBB->size_; i++) {
            const BB &bb1 = *symSBB->bbs_[i];
            RigidTrans3 t1 = (!symSBB->trans_[i]);
            for (unsigned int j = i + 1; j < symSBB->size_; j++) {
                const BB &bb2 = *symSBB->bbs_[j];
                RigidTrans3 t2 = t1 * symSBB->trans_[j];

                const BB *pBB1 = &bb1, *pBB2 = &bb2;
                unsigned int totalUsedAtoms = 0;
                unsigned int bbPenetrations = 0;

                for (Molecule<Atom>::const_iterator it = pBB2->caAtoms_.begin(); it != pBB2->caAtoms_.end(); it++) {
                    if (it->getTempFactor() < minTemperatureToConsiderCollision) {
                        continue;
                    }
                    totalUsedAtoms++;

                    Vector3 v = t2 * it->position();
                    if (pBB1->getDistFromSurface(v) < 0) {
                        // getResidueEntry(v) when used in BBGrid.h will return -1*res_index if res_index is backbone
                        if (pBB1->grid_->getResidueEntry(v) < 0 && pBB1->grid_->getDist(v) < -1.0) {
                            int resEntry = pBB1->grid_->getResidueEntry(v) * -1;
                            if (pBB1->getAtomByResId(resEntry).getTempFactor() < minTemperatureToConsiderCollision)
                                continue;
                            bbPenetrations++;
                        }
                    }
                }

                if ((bbPenetrations / (1.0 * totalUsedAtoms)) > 0.2) {
                    std::cout << "dropping " << identBBs.size() << " because penetration "
                              << bbPenetrations / (1.0 * totalUsedAtoms) << std::endl;
                    shouldContinuePen = true;
                    break;
                }

                maxPenetration = std::max(maxPenetration, bbPenetrations / (1.0 * totalUsedAtoms));
            }
            if (shouldContinuePen)
                break;
        }
        if (shouldContinuePen)
            continue;
        std::cout << "checked penetrations ratio max: " << maxPenetration << std::endl;
        // if above some TH (for everything, not per chain) (20%) - drop

        // if last and first centers are the farthest - drop
        std::vector<Vector3> centroids;
        for (unsigned int i = 0; i < symSBB->size_; i++) {
            centroids.push_back(symSBB->trans_[i] * symSBB->bbs_[i]->getCM());
        }
        std::cout << "centroids distance " << (centroids[0] - centroids[1]).norm2() << " : "
                  << (centroids[0] - centroids.back()).norm2() << std::endl;

        float allowedDistFactor = 1.5 + (symSBB->size_ - 3) * 0.25;
        if ((centroids[0] - centroids[1]).norm2() * allowedDistFactor < (centroids[0] - centroids.back()).norm2()) {
            std::cout << "dropping " << identBBs.size() << " because centroids distance "
                      << (centroids[0] - centroids[1]).norm2() << " : " << (centroids[0] - centroids.back()).norm2()
                      << std::endl;
            continue;
        }

        // verify that centroids are first all increasing distance from first centroid and then all decreasing distance
        // from first centroid
        bool increasing = true;
        bool shouldContinue = false;
        for (unsigned int i = 1; i < centroids.size(); i++) {
            if (increasing) {
                if ((centroids[i] - centroids[0]).norm2() < (centroids[i - 1] - centroids[0]).norm2()) {
                    increasing = false;
                }
            } else {
                if ((centroids[i] - centroids[0]).norm2() > (centroids[i - 1] - centroids[0]).norm2()) {
                    std::cout << "dropping " << identBBs.size() << " because centroids not increasing and decreasing"
                              << std::endl;
                    shouldContinue = true;
                    break;
                }
            }
        }
        if (shouldContinue)
            continue;
        if (increasing) {
            std::cout << "dropping " << identBBs.size() << " because centroids only increasing " << std::endl;
            continue;
        }

        std::cout << "added with score " << symSBB->weightedTransScore_ << std::endl;
        results.push(symSBB);
        addedSymCount++;
    }

    BitId groupIdentifier;
    for (unsigned int i = 0; i < identBBs.size(); i++) {
        groupIdentifier |= identBBs[i]->bbs_[0]->bitId();
    }
    std::cout << "Created " << addedSymCount << " Symmetrical for " << groupIdentifier << std::endl;
}

// utils
std::shared_ptr<SuperBB> HierarchicalFold::createJoined(const SuperBB &sbb1, const SuperBB &sbb2, RigidTrans3 &trans,
                                                        int bbPen, FoldStep &step, float transScore) const {
    std::shared_ptr<SuperBB> theNew = std::make_shared<SuperBB>(sbb1);
    theNew->join(trans, sbb2, bbPen, step, transScore);
    theNew->setRestraintsRatio(complexConst_.getRestraintsRatio(theNew->bbs_, theNew->trans_));
    return theNew;
}

std::shared_ptr<SuperBB> HierarchicalFold::getMatchingSBB(SuperBB sbb1, SuperBB sbb2,
                                                          std::vector<std::vector<unsigned int>> &identGroups) {
    /*
    This function recieves two SuperBBs and checks if they have common BBs that are a part of the same ident group.
    If so, it checks wether the total amount of BBs from the same ident group is smaller than the size of the ident 
    group. If so, It return a new SuperBB based on sbb2, but with the BBs from the ident group that are not in sbb1.
    
    It also validates that in both sbb1&sbb2, the BBs from the ident group are the smallest BBs in the ident group (This
    prevents duplications of results). If they are not valid - returns NULL.
    */
    std::map<unsigned int, unsigned int> bbIdToNewId;
    for (std::vector<unsigned int> identGroup : identGroups) {
        // verify sbb1 is valid (mostly needed to ignore initial structures of BBs that are not first in group) 
        bool flag = false;
        int maxIdInSbb1 = -1;
        for (unsigned int i = 0; i < identGroup.size(); i++) {
            if (!sbb1.bitIds().test(identGroup[i])) // ident_group[i] not in currResSet
                flag = true;
            else if (flag) {
                if (sbb1.bbs_.size() != 1)
                    std::cout << "sbb1 not valid " << sbb1.bitIds() << ":" << sbb2.bitIds() << std::endl;
                return NULL;
            } else {
                maxIdInSbb1 = i;
            }
        }

        // compute mapping from sbb2 bb ids to new bb ids
        int maxIdInSbb2 = -1;
        flag = false;
        for (unsigned int i = 0; i < identGroup.size(); i++) {
            if (!sbb2.bitIds().test(identGroup[i])) // ident_group[i] not in currResSet
                flag = true;
            else if (flag) {
                if (sbb2.bbs_.size() != 1)
                    std::cout << "sbb2 not valid " << sbb1.bitIds() << ":" << sbb2.bitIds() << std::endl;
                return NULL;
            } else {
                maxIdInSbb2 = i;
            }
        }

        if (maxIdInSbb1 == -1 || maxIdInSbb2 == -1)
            continue;

        // check if there are more copies in sbb1 and sbb2 than the size of the ident group
        if ((maxIdInSbb1 + 1) + (maxIdInSbb2 + 1) > identGroup.size()) { 
            return NULL;
        }

        for (unsigned int i = 0; i < maxIdInSbb2 + 1; i++) {
            bbIdToNewId[identGroup[i]] = identGroup[i + maxIdInSbb1 + 1];
        }
    }

    // create new SuperBB with new ids
    if (bbIdToNewId.size() == 0)
        return std::make_shared<SuperBB>(sbb2);

    // std::cout << "converting " << sbb1.bitIds() << ":" << sbb2.bitIds() << " with " << bbIdToNewId.size() <<
    // std::endl;
    std::shared_ptr<SuperBB> newSbb = std::make_shared<SuperBB>(sbb2);

    for (auto iter = bbIdToNewId.rbegin(); iter != bbIdToNewId.rend(); ++iter) {
        unsigned int oldId = iter->first;
        unsigned int newId = iter->second;
        newSbb->replaceIdentBB(BitId(oldId), (*(bestKContainer_[BitId(newId)].begin()))->bbs_[0]);
    }

    return newSbb;
}

void HierarchicalFold::checkConnectivity() const {
    typedef boost::adjacency_list<boost::vecS,        // edge list
                                  boost::vecS,        // vertex list
                                  boost::undirectedS, // directedness
                                  float>              // property associated with vertices
        Graph;

    Graph g(N_);
    for (unsigned int suIndex = 0; suIndex < N_; suIndex++) {
        BitId index = BitId(suIndex);
        std::cerr << "isEmpty " << suIndex << " " << bestKContainer_.isEmpty(index) << std::endl;
        std::cerr << "SBB " << suIndex << " size " << bestKContainer_[index].size() << std::endl;
        if (bestKContainer_[index].size() >= 1) {
            std::shared_ptr<SuperBB> sbb = *(bestKContainer_[index].rbegin());
            std::shared_ptr<const BB> bb = sbb->bbs_[0];

            for (unsigned int suIndex2 = 0; suIndex2 < N_; suIndex2++) {
                unsigned int transSize = bb->getTransformations(suIndex2).size();
                std::cerr << suIndex << " " << suIndex2 << " trans size " << transSize << std::endl;
                if (transSize > 0)
                    boost::add_edge(suIndex, suIndex2, g);
                std::cerr << "done add_edge" << std::endl;
            }
        }
    }

    std::vector<int> component(boost::num_vertices(g));
    size_t num_components = boost::connected_components(g, &component[0]);
    std::cerr << "Num of components " << num_components << std::endl;
    if (num_components != 1) {
        for (size_t i = 0; i < boost::num_vertices(g); ++i)
            std::cerr << "SU " << i << " component " << component[i] << std::endl;
        std::cerr << "Not enough transformations between subunits, there should be one connected component!"
                  << std::endl;
        exit(1);
    }
}

void HierarchicalFold::outputConnectivityGraph(std::string outFileName) const {
    std::ofstream outFile(outFileName);
    outFile << "SU1 Prot1 size1 SU2 Prot2 size2 Restraints" << std::endl;
    for (unsigned int suIndex = 0; suIndex < N_; suIndex++) {
        BitId set = BitId(suIndex);
        std::cerr << "SBB " << suIndex << " set " << set << " isEmpty " << bestKContainer_.isEmpty(set) << std::endl;
        std::cerr << "SBB " << suIndex << " size " << bestKContainer_[set].size() << std::endl;
        if (bestKContainer_[set].size() >= 1) {
            std::shared_ptr<SuperBB> sbb = *(bestKContainer_[set].rbegin());
            std::shared_ptr<const BB> bb = sbb->bbs_[0];

            for (unsigned int suIndex2 = suIndex + 1; suIndex2 < N_; suIndex2++) {
                unsigned int transSize = bb->getTransformations(suIndex2).size();
                BitId set2 = BitId(suIndex2);
                std::cerr << suIndex << " " << suIndex2 << " trans size " << transSize << " " << set2 << std::endl;
                std::shared_ptr<SuperBB> sbb2 = *(bestKContainer_[set2].rbegin());
                std::shared_ptr<const BB> bb2 = sbb2->bbs_[0];

                if (transSize > 0) {
                    std::vector<std::string> results1, results2;
                    std::string PDBFileName1 = bb->getPDBFileName();
                    boost::split(results1, PDBFileName1, [](char c) { return c == '_'; });
                    std::string PDBFileName2 = bb2->getPDBFileName();
                    boost::split(results2, PDBFileName2, [](char c) { return c == '_'; });

                    unsigned int xlinkNumber = complexConst_.numberOfRestraints(suIndex, suIndex2);
                    unsigned int connectivityConstraintsNumber = complexConst_.numberOfConstraints(suIndex, suIndex2);

                    if (xlinkNumber > 0) {
                        outFile << bb->getPDBFileName() << " " << results1[0] << " " << bb->backBone_.size() / 4 << " ";
                        outFile << bb2->getPDBFileName() << " " << results2[0] << " " << bb2->backBone_.size() / 4
                                << " ";
                        outFile << complexConst_.numberOfRestraints(suIndex, suIndex2) << std::endl;
                    }
                    if (connectivityConstraintsNumber > 0) {
                        outFile << bb->getPDBFileName() << " " << results1[0] << " " << bb->backBone_.size() / 4 << " ";
                        outFile << bb2->getPDBFileName() << " " << results2[0] << " " << bb2->backBone_.size() / 4
                                << " ";
                        outFile << complexConst_.numberOfConstraints(suIndex, suIndex2) << std::endl;
                    }
                    /*outFile << bb->getPDBFileName() << " "
                            << complexConst_.numberOfRestraints(suIndex, suIndex2) << " "
                            << bb2->getPDBFileName() << std::endl;*/
                }
            }
        }
    }
    outFile.close();
}
