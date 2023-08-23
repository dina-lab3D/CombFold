#ifndef BESTK_H
#define BESTK_H

#include "SuperBB.h"
#include <mutex>

static float scoreSuperBB(const std::shared_ptr<SuperBB> &sbb) {
    // return sbb->getRestraintsRatio();
    // return sbb->weightedTransScore_;

    // Note: if you change this, you also needs to change in HierarchicalFold::tryToConnect which optimizes by 
    // summing and comparing trans scores before trying to connect
    // return sbb->transScore_;
    return sbb->transScore_ * sbb->getRestraintsRatio();
}
struct comp {
    bool operator()(const std::shared_ptr<SuperBB> &lhs, const std::shared_ptr<SuperBB> &rhs) const {
        return scoreSuperBB(lhs) < scoreSuperBB(rhs);
    }
};

/**
   This class stores the best k permutations of a specific size
   implemented as inheriting from a multiset with shapred_ptr and
   the comp struct
*/
class BestK : public std::multiset<std::shared_ptr<SuperBB>, comp> {

  public:
    BestK(unsigned int k, bool toDel = true) : k_(k), curMinScore(-1) {}

    float score(const std::shared_ptr<SuperBB> &sbb) const { return scoreSuperBB(sbb); }

    float minScore() const { return curMinScore; }
    float maxScore() const { 
        if (size() == 0)
            return 0;
        return score(*rbegin()); 
    }

    void setK(int k) { k_ = k; }

    bool push(std::shared_ptr<SuperBB> in);

    bool push_cluster(std::shared_ptr<SuperBB> in, double rmsd, std::vector<std::vector<unsigned int>> &identGroups);

    void cluster(BestK &clusteredBest, double rmsd, std::vector<std::vector<unsigned int>> &identGroups) const;
    virtual ~BestK() {}

  private:
    float internalMinScore() {
        if (size() < k_)
            return -1;
        return score(*begin());
    }

  private:
    std::mutex _mu;
    unsigned int k_;
    float curMinScore;
};

#endif /* BESTK_H */
