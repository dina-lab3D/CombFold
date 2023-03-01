#ifndef BESTK_H
#define BESTK_H

#include "SuperBB.h"
#include <mutex>

/**
 *This struct is used for comparison of two SuperBBs
 */
struct comp {
    bool operator()(const std::shared_ptr<SuperBB> &lhs, const std::shared_ptr<SuperBB> &rhs) const {

        // return lhs->weightedTransScore_ < rhs->weightedTransScore_;
        return lhs->transScore_ < rhs->transScore_;
    }
};

/**
   This class stores the best k permutations of a specific size
   implemented as inheriting from a multiset with shapred_ptr and
   the comp struct
*/
class BestK : public std::multiset<std::shared_ptr<SuperBB>, comp> {

  public:
    BestK(unsigned int k, bool toDel = true) : k_(k), toDel_(toDel), curMinScore(-1) {}

    BestK() : k_(500), toDel_(true) {}

    BestK(BestK &&other) : std::multiset<std::shared_ptr<SuperBB>, comp>(other) {
        std::lock_guard<std::mutex> locker(_mu);
        toDel_ = other.toDel_;
        other.toDel_ = false;
        k_ = other.k_;
        curMinScore = other.curMinScore;
    }

    BestK &operator=(BestK &&other) {
        std::lock_guard<std::mutex> locker(_mu);
        *(std::multiset<std::shared_ptr<SuperBB>, comp> *)this =
            *(std::multiset<std::shared_ptr<SuperBB>, comp> *)(&other);
        toDel_ = other.toDel_;
        other.toDel_ = false;
        k_ = other.k_;
        curMinScore = other.curMinScore;
        return *this;
    }
    
    // Note: if you change this also needs to change in HierarchicalFold::tryToConnect
    float score(const std::shared_ptr<SuperBB> &sbb) const { return sbb->transScore_; }
    // float score(const std::shared_ptr<SuperBB> &sbb) const { return sbb->weightedTransScore_; }

    float minScore() const { return curMinScore; }
    float maxScore() const { 
        if (size() == 0)
            return 0;
        return score(*rbegin()); 
    }

    void setK(int k) { k_ = k; }

    bool push(std::shared_ptr<SuperBB> in) {
        std::lock_guard<std::mutex> locker(_mu);
        if (size() < k_) {
            insert(in);
            return true;
        } else {

            if (internalMinScore() < score(in)) {
                erase(begin());
                insert(in);
                curMinScore = score(*begin());
                return true;
            }
        }

        return false;
    }

    bool push_cluster(std::shared_ptr<SuperBB> in, double rmsd, std::vector<std::vector<unsigned int>> &identGroups) {
        std::lock_guard<std::mutex> locker(_mu);

        if (internalMinScore() > score(in))
            return false;

        for (auto it = begin(); it != end(); it++) {
            if (score(in) <= score(*it) && (*it)->calcRmsd(*in, identGroups) < rmsd)
                return false;
        }

        for (auto it = begin(); it != end();) {
            if (score(in) > score(*it) && (*it)->calcRmsd(*in, identGroups) < rmsd) {
                erase(it++);
            } else {
                ++it;
            }
        }

        if (size() >= k_) {
            erase(begin());
        }

        insert(in);
        curMinScore = score(*begin());
        return true;
    }

    void cluster(BestK &clusteredBest, double rmsd, std::vector<std::vector<unsigned int>> &identGroups) const {

        if (size() == 0)
            return;
        const std::shared_ptr<SuperBB> firstSBB = *rbegin();
        if (firstSBB->size() == 2) { // don't cluster
            for (auto it = rbegin(); it != rend(); it++) {
                clusteredBest.insert(*it);
            }
        } else {
            std::vector<bool> clustered(size(), false);
            int i = 0;

            for (auto it = rbegin(); it != rend(); it++, i++) {
                const std::shared_ptr<SuperBB> refSBB = *it;
                // std::cout << "clustering " << i << " clsutered:" << clustered[i] << std::endl;
                if (!clustered[i]) {
                    clustered[i] = true;
                    clusteredBest.insert(refSBB);
                }

                // TODO: maybe shouldn't cluster more if already clustered (can lead to drift)

                // cluster to other SBBs
                int j = 0;
                // for(auto it2 = begin(); it2!= end(); it2++, j++) {
                for (auto it2 = rbegin(); it2 != rend(); it2++, j++) {
                    if (!clustered[j] && (*it2)->calcRmsd(*refSBB, identGroups) < rmsd) {
                        clustered[j] = true;
                    }
                }
            }
        }
    }

    virtual ~BestK() {}

  private:
    float internalMinScore() {
        if (size() < k_)
            return -1;
        return score(*begin());
    }

    BestK(const BestK &other);
    BestK &operator=(const BestK &other);

  private:
    std::mutex _mu;
    unsigned int k_;
    bool toDel_;
    float curMinScore;
};

#endif /* BESTK_H */
