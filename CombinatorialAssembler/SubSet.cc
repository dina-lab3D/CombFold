#include "SubSet.h"

SubSet::SubSet(unsigned int N, unsigned int K) : N_(N), K_(K) {
  index_.push_back(1);
  for (unsigned int i = 1; i < K_; i++) {
    index_.push_back(index_[i-1]*2);
  }
  unsigned long max = 1;

  for (unsigned int i = 0; i < N_; i++) {
    if (i >= N_- K_) {
      max_.push_back(max);
      //cout << "max " << max << endl;
    }
    max *= 2;
  }
  index_.push_back(max);
  MAX_ = max;
  MAX_K = index_[K-1] * 2;
}

bool SubSet::increase() {
  int i = K_-1;
  while(i >= 0) {
    index_[i] *= 2;
    if (index_[i] <= max_[i])
      break;
    if (i == 0)
      return false;
    i--;
  }
  //cout << "i " << i << " index_[i] " << index_[i] << " max_[i] " << max_[i] <<endl;
  for (i++;i < (int)K_; i++) {
    index_[i] = index_[i - 1]*2;
  }
  return true;
}

unsigned long SubSet::getSubset(unsigned long n) {
  long ret = 0;
  for (unsigned long i = 1,j=0; i <= n; i*=2) {
    if ((n&i) > 0) ret |= index_[j];
    j++;
  }
  return ret;
}

unsigned long SubSet::getComplementry(unsigned long n) {
  long ret = 0;
  for (unsigned long i = 1,j=0; i < MAX_K; i*=2) {
    if ((n&i) == 0) ret |= index_[j];
    j++;
  }
  return ret;
}
