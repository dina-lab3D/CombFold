#ifndef SUBSET_H
#define SUBSET_H

#include <vector>

class SubSet {
  public:
    SubSet(unsigned int N, unsigned int K);
    bool increase();
    unsigned long getSubset(unsigned long n);
    unsigned long getComplementry(unsigned long n);
    unsigned long max() { return MAX_; }
    unsigned long maxK() { return MAX_K; }

  protected:
    unsigned int N_, K_;
    unsigned long MAX_, MAX_K;
    std::vector<unsigned long> index_;
    std::vector<unsigned long> max_;
};

#endif
