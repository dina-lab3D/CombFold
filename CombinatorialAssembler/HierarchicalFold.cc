/*
  This class creates superBBs in the fold function
*/

#include "HierarchicalFold.h"
#include "SubSet.h"

#include "ctpl_stl.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

Timer HierarchicalFold::timer_;
Timer HierarchicalFold::timerAll_;
unsigned int HierarchicalFold::countResults_(0);

void HierarchicalFold::checkConnectivity() const {
    typedef boost::adjacency_list<boost::vecS,        // edge list
                                  boost::vecS,        // vertex list
                                  boost::undirectedS, // directedness
                                  float>              // property associated with vertices
        Graph;

    Graph g(N_);
    for (unsigned int suIndex = 0; suIndex < N_; suIndex++) {
        unsigned long one = 1;
        unsigned long index = one << suIndex;
        std::cerr << "isEmpty " << suIndex << " " << clusteredSBBS_.isEmpty(index) << std::endl;
        std::cerr << "SBB " << suIndex << " size " << clusteredSBBS_[index].size() << std::endl;
        if (clusteredSBBS_[index].size() >= 1) {
            std::shared_ptr<SuperBB> sbb = *(clusteredSBBS_[index].rbegin());
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
        unsigned long one = 1;
        unsigned long set = one << suIndex;
        std::cerr << "SBB " << suIndex << " set " << set << " isEmpty " << clusteredSBBS_.isEmpty(set) << std::endl;
        std::cerr << "SBB " << suIndex << " size " << clusteredSBBS_[set].size() << std::endl;
        if (clusteredSBBS_[set].size() >= 1) {
            std::shared_ptr<SuperBB> sbb = *(clusteredSBBS_[set].rbegin());
            std::shared_ptr<const BB> bb = sbb->bbs_[0];

            for (unsigned int suIndex2 = suIndex + 1; suIndex2 < N_; suIndex2++) {
                unsigned int transSize = bb->getTransformations(suIndex2).size();
                unsigned long set2 = one << suIndex2;
                std::cerr << suIndex << " " << suIndex2 << " trans size " << transSize << " " << set2 << std::endl;
                std::shared_ptr<SuperBB> sbb2 = *(clusteredSBBS_[set2].rbegin());
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

int get_available_concurrency() {
    if (std::getenv("SLURM_CPUS_PER_TASK") != NULL)
        return std::atoi(std::getenv("SLURM_CPUS_PER_TASK"));

    return std::thread::hardware_concurrency();
}

unsigned nChoosek(unsigned int n, unsigned int k) {
    if (k > n)
        return 0;
    if (k * 2 > n)
        k = n - k;
    if (k == 0)
        return 1;

    int result = n;
    for (unsigned int i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

bool isValidBasedOnIdent(std::vector<std::vector<unsigned int>> &identGroups, unsigned long currResSet) {
    // Validating jobs will prevent running the same jobs for different ident (for example if A=B=C,D=E) we wouldn't
    // need to check both AD and AE since they are the same so only AD will be valid.
    // we will enable BE (but not AE/CE) so we can still have fold steps of AD+BE=ADBE
    // we still have some repetition, for example we will create both AB and BC to enable later AD+BC

    /*
     * def is_valid(ident_groups, resSet):
    ...:     baseline = 999
    ...:     for ident_group in ident_groups:
    ...:         for i in range(len(ident_group)):
    ...:             if ident_group[i] in resSet:
    ...:                 baseline = min(baseline, i)
    ...:                 break
    ...:     if baseline == 999:
    ...:         return True
    ...:     for ident_group in ident_groups:
    ...:         flag = False
    ...:         for i in range(baseline, len(ident_group)):
    ...:             if ident_group[i] not in resSet:
    ...:                 flag = True
    ...:             elif flag:
    ...:                 return False
    ...:     return True
    ...:
     */
    unsigned int baseline = 999;
    for (std::vector<unsigned int> identGroup : identGroups) {
        for (unsigned int i = 0; i < identGroup.size(); i++) {
            if (((1 << identGroup[i]) & currResSet) != 0)
                if (baseline > i)
                    baseline = i;
        }
    }
    if (baseline == 999)
        return true;

    for (std::vector<unsigned int> identGroup : identGroups) {
        bool flag = false;
        for (unsigned int i = baseline; i < identGroup.size(); i++) {
            if (((1 << identGroup[i]) & currResSet) == 0) // ident_group[i] not in currResSet
                flag = true;
            else if (flag) {
                // std::cout << "not valid " << currResSet << " based on " << i << std::endl;
                return false;
            }
        }
    }
    return true;
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
                // float transScore = 95 + it.getScore() / 20;
                float transScore = it.getScore() + it.getScore() / 5;

                FoldStep step(symSBB->bbs_[i - 1]->getID(), identBBs[i]->bbs_[0]->getID(), transScore);
                std::cout << "adding trans " << it.transformation() << " **** " << it.getScore() << std::endl;
                // symSBB = createJoined(*symSBB, *identBBs[i], it.transformation(), 0, 0, 0, step, it.getScore());
                symSBB = createJoined(*symSBB, *identBBs[i], it.transformation(), 0, 0, 0, step, transScore);
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


    unsigned long groupIdentifier = 0;
    for (unsigned int i = 0; i < identBBs.size(); i++) {
        groupIdentifier += identBBs[i]->bbs_[0]->bitId();
    }
    std::cout << "Created " << addedSymCount << " Symmetrical for " << groupIdentifier << std::endl;
}

void HierarchicalFold::fold(const std::string &outFileNamePrefix) {
    const size_t nthreads = get_available_concurrency();
    std::cout << "parallel (" << nthreads << " threads):" << std::endl;
    ctpl::thread_pool threadPool(nthreads);
    bool output = false; // did we find something for N or N-1 subunits
    std::string outFileName = outFileNamePrefix + ".res";
    std::ofstream outFile(outFileName);
    std::ofstream outFileClustered(outFileNamePrefix + "_clustered.res");

    std::map<unsigned int, BestK *> keptResultsByLength;
    keptResultsByLength[1] = new BestK(N_);

    //  BestK* keptResults = new BestK(N_);
    for (unsigned int i = 0; i < N_; i++) {
        BestK &singleChainSBBS = clusteredSBBS_[1 << i];
        if (singleChainSBBS.size() != 1) {
            std::cerr << "started with " << singleChainSBBS.size() << " SBBs for " << (1 << i) << std::endl;
        }

        for (auto it1 = singleChainSBBS.begin(); it1 != singleChainSBBS.end(); it1++) {
            keptResultsByLength[1]->push(*it1);
            //          keptResults->push(*it1);
        }
        std::cout << "added to kept results chain " << i << std::endl;
    }

    // create ident chains
    std::vector<bool> addedToGroup(N_, false);
    std::vector<std::vector<unsigned int>> identGroups;
    for (unsigned int i = 0; i < N_; i++) {
        if (addedToGroup[i])
            continue;
        addedToGroup[i] = true;
        std::vector<unsigned int> identical;
        std::shared_ptr<SuperBB> sbbI = *(clusteredSBBS_[1 << i].begin());
        for (unsigned int j = i + 1; j < N_; j++) {
            std::cout << "checking ident " << i << " & " << j << std::endl;
            std::shared_ptr<SuperBB> sbbJ = *(clusteredSBBS_[1 << j].begin());

            const std::shared_ptr<const BB> bbI = sbbI->bbs_[0];
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

    // populate with homomers subunits
    for (unsigned int i = 2; i <= N_; i++)
        keptResultsByLength[i] = new BestK(K_);

    for (std::vector<unsigned int> identGroup : identGroups) {
        for (unsigned int groupDivider = 1; identGroup.size() / groupDivider >= 5; groupDivider++) {
            if ((identGroup.size() % groupDivider) != 0)
                continue;
            unsigned int groupSize = identGroup.size() / groupDivider;
            for (unsigned int i = 0; i < groupDivider; i++) {
                std::vector<std::shared_ptr<SuperBB>> groupSBBs;
                for (unsigned int j = 0; j < groupSize; j++) {
                    groupSBBs.push_back(*(clusteredSBBS_[1 << identGroup[i * groupSize + j]].begin()));
                }
                std::cout << "searching for size " << groupSBBs.size() << " has "
                          << keptResultsByLength.count(groupSBBs.size()) << std::endl;
                createSymmetry(groupSBBs, *keptResultsByLength[groupSBBs.size()]);
            }
        }
    }

    for (unsigned int length = 2; length <= N_; length++) { // # subunits iteration
        std::cout << "*** running iteration " << length
                  << " prev kept results: " << keptResultsByLength[length - 1]->size() << std::endl;
        BestK *newKept = keptResultsByLength[length];
        std::map<unsigned long, BestK *> best_k_by_id;

        for (auto it1 = newKept->begin(); it1 != newKept->end(); it1++) {
            unsigned long currResSet = (**it1).bitIds();
            if (best_k_by_id.count(currResSet) == 0)
                best_k_by_id[currResSet] = new BestK(K_);
            best_k_by_id[currResSet]->push(*it1);
        }
        for (unsigned int firstResultSize = 1; firstResultSize <= length / 2; firstResultSize++) {
            unsigned int secondResultSize = length - firstResultSize;
            std::cout << "** running sub-iteration " << firstResultSize << " " << secondResultSize << std::endl;

            for (auto it1 = keptResultsByLength[firstResultSize]->begin();
                 it1 != keptResultsByLength[firstResultSize]->end(); it1++) {
                unsigned long setA = (**it1).bitIds();
                for (auto it2 = keptResultsByLength[secondResultSize]->begin();
                     it2 != keptResultsByLength[secondResultSize]->end(); it2++) {
                    unsigned long setB = (**it2).bitIds();
                    if ((setA & setB) != 0)
                        continue;
                    if (firstResultSize == secondResultSize && setA < setB)
                        continue; // Since those are 2 identical loops, don't do things twice

                    unsigned long currResSet = setA | setB;

                    if (!isValidBasedOnIdent(identGroups, currResSet))
                        continue;

                    if (best_k_by_id.count(currResSet) == 0)
                        best_k_by_id[currResSet] = new BestK(K_);

                    unsigned int resCountBefore = best_k_by_id[currResSet]->size();

                    std::promise<int> promise1;
                    this->tryToConnect(1, **it1, **it2, *best_k_by_id[currResSet], (length < N_), promise1,
                                       identGroups);

                    if (resCountBefore < best_k_by_id[currResSet]->size())
                        std::cout << "found more for " << currResSet << " based on " << setA << " and " << setB
                                  << " before: " << resCountBefore << " after: " << best_k_by_id[currResSet]->size()
                                  << " min_score: " << best_k_by_id[currResSet]->minScore() << std::endl;
                }
            }
        }
        for (const auto &[currResSet, currBestK] : best_k_by_id) {
            if (currBestK->size() > 0) {
                BestK *clusteredBestK = clusteredSBBS_.newBestK(currResSet);

                if (length == N_) { // output before clustering
                    output = true;
                    clusteredBestK->setK(finalSizeLimit_);
                    for (auto it = currBestK->rbegin(); it != currBestK->rend(); it++) {
                        (*it)->calcFinalScore();
                        (*it)->setRestraintsRatio(complexConst_.getRatio((*it)->bbs_, (*it)->trans_));
                        (*it)->fullReport(outFile);
                    }
                    // output after clustering
                    currBestK->cluster(*clusteredBestK, 5.0, identGroups);
                    for (auto it = clusteredBestK->rbegin(); it != clusteredBestK->rend(); it++) {
                        (*it)->calcFinalScore();
                        (*it)->setRestraintsRatio(complexConst_.getRatio((*it)->bbs_, (*it)->trans_));
                        (*it)->fullReport(outFileClustered);
                    }
                } else {
                    currBestK->cluster(*clusteredBestK, 1.0, identGroups);
                    // for(auto it=currBestK->rbegin(); it != currBestK->rend(); it++)
                    //     clusteredBestK->insert(*it);
                }
                std::cerr << "clustering resSet " << currResSet << " before: " << currBestK->size() << " after "
                          << clusteredSBBS_[currResSet].size() << std::endl;

                unsigned int count = 0;
                for (auto it1 = clusteredSBBS_[currResSet].begin(); it1 != clusteredSBBS_[currResSet].end(); it1++) {
                    newKept->push(*it1);
                    count += 1;
                    if (count >= maxResultPerResSet)
                        break;
                }
            }
            // TODO: delete currBestK and dictionary
            //        delete currBestK;
        }

        // just for print
        std::map<unsigned int, unsigned int> count_new_kept_by_bb;
        std::map<unsigned int, unsigned int> count_new_kept_by_resSet;
        for (unsigned int bit_index = 0; bit_index < N_; bit_index++)
            count_new_kept_by_bb[bit_index] = 0;
        for (auto it1 = newKept->begin(); it1 != newKept->end(); it1++) {
            if (count_new_kept_by_resSet.count((**it1).bitIds()) == 0)
                count_new_kept_by_resSet[(**it1).bitIds()] = 0;
            count_new_kept_by_resSet[(**it1).bitIds()] += 1;
            for (unsigned int bit_index = 0; bit_index < N_; bit_index++)
                if (((**it1).bitIds() & (1 << bit_index)) > 0) {
                    count_new_kept_by_bb[bit_index] += 1;
                }
        }

        std::cout << "new kept results by chain ";
        for (const auto &elem : count_new_kept_by_bb)
            std::cout << elem.first << ":" << elem.second << ", ";
        std::cout << std::endl << "new kept results by resSet ";
        for (const auto &elem : count_new_kept_by_resSet)
            std::cout << elem.first << ":" << elem.second << ", ";
        std::cout << std::endl;
    }

    outFile.close();

    if (!output)
        outputSubsets();
}

void HierarchicalFold::outputSubsets() const {
    bool found = false;
    for (unsigned int length = N_ - 1; length > 2; length--) { // # subunits iteration

        SubSet subsetGen(N_, length); // generates subsets of size length

        do { // iterate subsets of length out of N
            // unsigned long max = subsetGen.maxK() / 2;
            unsigned long one = 1;
            unsigned long currResSet = subsetGen.getSubset(one) | subsetGen.getComplementry(one);

            if (!clusteredSBBS_.isEmpty(currResSet)) {
                found = true;
                std::string outFileName = "cb_" + std::to_string(length) + "_" + std::to_string(currResSet) + ".res";
                std::cout << outFileName << " missing " << findMissingSubunits(currResSet) << std::endl;
                std::ofstream oFile(outFileName);
                for (auto it = clusteredSBBS_[currResSet].rbegin(); it != clusteredSBBS_[currResSet].rend(); it++) {
                    (*it)->calcFinalScore();
                    (*it)->setRestraintsRatio(complexConst_.getRatio((*it)->bbs_, (*it)->trans_));
                    (*it)->fullReport(oFile);
                }
                oFile.close();
            }

        } while (subsetGen.increase()); // end subset iteration

        if (found)
            break;
    }
}

std::string HierarchicalFold::findMissingSubunits(unsigned long currResSet) const {
    std::string ret;
    for (unsigned int suIndex = 0; suIndex < N_; suIndex++) {
        unsigned long one = 1;
        unsigned long set = one << suIndex;
        if ((currResSet | set) != currResSet) {
            std::shared_ptr<SuperBB> sbb = *(clusteredSBBS_[set].rbegin());
            std::shared_ptr<const BB> bb = sbb->bbs_[0];
            ret += bb->getPDBFileName();
            ret += " ";
        }
    }
    // std::cerr << "Can't find missing subunit " << currResSet << std::endl;
    return ret;
}

void HierarchicalFold::tryToConnect(int id, const SuperBB &sbb1, const SuperBB &sbb2, BestK &results, bool toAdd,
                                    std::promise<int> &output, std::vector<std::vector<unsigned int>> &identGroups) {

    float penetrationThreshold_ = BB::penetrationThreshold_;

    unsigned int totalCa = 0;
    for (int i = 0; i < (int)sbb1.bbs_.size(); i++)
        totalCa += sbb1.bbs_[i]->caAtoms_.size();
    for (int i = 0; i < (int)sbb2.bbs_.size(); i++)
        totalCa += sbb2.bbs_[i]->caAtoms_.size();

    //  int backbonePenetrationThreshold = totalCa / 20;
    //    int backbonePenetrationThreshold = (int) (0.01 * BB::backbonePenetrationMultiplierThreshold_ * totalCa);
    int backbonePenetrationThreshold = (int)totalCa;

    // iterate BBs os SuperBB1
    for (int i = 0; i < (int)sbb1.bbs_.size(); i++) {
        int firstBB = sbb1.bbs_[i]->getID();

        // iterate BBs os SuperBB1
        for (int j = 0; j < (int)sbb2.bbs_.size(); j++) {
            int secondBB = sbb2.bbs_[j]->getID();
            // if (!isValid(sbb1, sbb2, firstBB, secondBB)) continue;

            // check backbone penetration
            //      int bbPen = sbb1.backBonePen_ + sbb2.backBonePen_;
            //      if (bbPen > backbonePenetrationThreshold){
            //          std::cout << "filtered because joined bbPen" << std::endl;
            //          continue;
            //      }

            // loop over possible transformations between BBs
            for (TransIterator2 it(sbb1, sbb2, firstBB, secondBB); !it.isAtEnd(); it++) {
                //        std::cout << "trans check started... ";

                // check that the score is not lower than the minimum in the current bestK
                //        int posibleScore = (int) it.getScore() + sbb1.transScore_ + sbb2.transScore_;
                //        if (posibleScore < results.minScore()){
                //            // std::cout << "filtered because score too low" << std::endl;
                //            continue;
                //        }

                int bbPen = sbb1.backBonePen_ + sbb2.backBonePen_;
                float newScore = 0.0;
                // discard any invalid transformations
                bool filtered = filterTrans(sbb1, sbb2, it.transformation(), penetrationThreshold_, bbPen,
                                            backbonePenetrationThreshold, newScore, i, j);
                if (filtered)
                    continue;

                int numBurried = 0; // numBurried2;// + numBurried1;
                FoldStep step(firstBB, secondBB, it.getScore());
                std::shared_ptr<SuperBB> theNew =
                    createJoined(sbb1, sbb2, it.transformation(), numBurried, bbPen, newScore, step, it.getScore());

                if (theNew->getRestraintsRatio() < complexConst_.getDistanceRestraintsRatio()) {
                    //            std::cout << "not enough restraints " << theNew->getRestraintsRatio() << " : " <<
                    //            complexConst_.getDistanceRestraintsRatio();
                    continue;
                }

                results.push(theNew);
                // results.push_cluster(theNew, 1, identGroups);
            }
        }
    }
    output.set_value(1);
}

std::shared_ptr<SuperBB> HierarchicalFold::createJoined(const SuperBB &sbb1, const SuperBB &sbb2, RigidTrans3 &trans,
                                                        int numBurried, int bbPen, float newScore, FoldStep &step,
                                                        int transScore) const { // first esume i < j
    std::shared_ptr<SuperBB> theNew = std::make_shared<SuperBB>(sbb1);
    theNew->join(trans, sbb2, numBurried, bbPen, newScore, step, transScore);
    theNew->setRestraintsRatio(complexConst_.getRatio(theNew->bbs_, theNew->trans_));
    return theNew;
}

bool HierarchicalFold::filterTrans(const SuperBB &sbb1, const SuperBB &sbb2, const RigidTrans3 &trans,
                                   float penThreshold, int &bbPenetrations, int maxBBPenetrations, float &newScore,
                                   unsigned int sbb1Index, unsigned int sbb2Index) const {

    // check distance constraints & restraints
    for (unsigned int i = 0; i < sbb1.size_; i++) {
        const BB &bb1 = *sbb1.bbs_[i];
        RigidTrans3 t = (!sbb1.trans_[i]) * trans;
        for (unsigned int j = 0; j < sbb2.size_; j++) {
            const BB &bb2 = *sbb2.bbs_[j];
            RigidTrans3 t2 = t * sbb2.trans_[j];
            // check constraints first
            if (!complexConst_.isSatisfied(bb1.getID(), bb2.getID(), t2))
                return true;
            // check restraints
            //      if (!complexConst_.isRatioSatisfied(bb1.getID(), bb2.getID(), t2)) return true;
        }
    }

    // backbone penetrations
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

            unsigned int bbPenetrationsBefore = bbPenetrations;

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
                    if (pBB1->grid_->getResidueEntry(v) < 0 && pBB1->grid_->getDist(v) < -1.0) {
                        int resEntry = pBB1->grid_->getResidueEntry(v) * -1;
                        if (pBB1->getAtomByResId(resEntry).getTempFactor() < minTemperatureToConsiderCollision)
                            continue;
                        bbPenetrations++;
                        //                      if (bbPenetrations > maxBBPenetrations) {
                        //                          // std::cout << " bbPenetrations " << bbPenetrations << " " <<
                        //                          maxBBPenetrations; return true;
                        //                      }
                    }
                }
            }

            // TODO: should this threshold be configurable
            // float bbPenChangePercent = (float)(bbPenetrations - bbPenetrationsBefore) / (float)pBB2->caAtoms_.size();
            float bbPenChangePercent = (float)(bbPenetrations - bbPenetrationsBefore) / (float)totalUsedAtoms;
            if (bbPenChangePercent > maxBackboneCollisionPercentPerChain) {
                return true;
            }
        }
    }

    return false;
}
