/**
 * \file BestKContainer.h
 * \brief
 *
 * \authors Dina Schneidman
 *
 */
#ifndef BESTKCONTAINER_H
#define BESTKCONTAINER_H

#include <unordered_map>

class BestKContainer : private std::vector<BestK *> {
  public:
    // N - subunit number
    BestKContainer(int K) : K_(K) {}

    ~BestKContainer() {
        for (auto it = begin(); it != end(); it++)
            delete *it;
    }

    BestK *newBestK(unsigned long set) {
        BestK *best = new BestK(K_);
        push_back(best);
        set2index_[set] = size() - 1;
        return best;
    }

    bool isEmpty(const unsigned long set) const { return (set2index_.find(set) == set2index_.end()); }

    // assumes BestK for set exists, can be checked with isEmpty
    const BestK &operator[](const unsigned long set) const {
        unsigned int index = set2index_.find(set)->second;
        return *((std::vector<BestK *>)(*this))[index];
    }

    // assumes BestK for set exists
    BestK &operator[](const unsigned long set) { return *((std::vector<BestK *>)(*this))[set2index_[set]]; }

  private:
    int K_;
    // The long is a unique representation key for each SBB
    // And it is the long type because we want to support more than 32 subunits
    // int corresponds to an index in the vector that holds the BestK
    std::unordered_map<unsigned long, unsigned int> set2index_;
};

#endif /* BESTKCONTAINER_H */
