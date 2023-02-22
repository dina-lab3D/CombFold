#ifndef COMBDOCK_H
#define COMBDOCK_H

#include "BB.h"
#include "SuperBB.h"

#include <unordered_map>

#define MAXIMIZE(a,b) if (a < b) a = b;
#define MINIMIZE(a,b) if (a > b) a = b;

class ClusterSBBs {

  class Cluster {
  public:
    Cluster() : representive_(NULL), count_(0){}
    Cluster(SuperBB* sbb) : representive_(sbb), count_(1) {}

    bool add(SuperBB* sbb);

    static void max(SuperBB& sbb1, SuperBB& sbb2) {
      if (sbb1.singlePen_ - sbb1.multPen_ > sbb2.singlePen_ - sbb2.multPen_) {
        sbb1.singlePen_ = sbb2.singlePen_;
        sbb1.multPen_ = sbb2.multPen_;
      }
      if (((float)sbb1.numOfHydBuried_) / (sbb1.numOfBuried_ - sbb1.numOfHydBuried_) >
          ((float)sbb1.numOfHydBuried_) / (sbb1.numOfBuried_ - sbb1.numOfHydBuried_)) {
        sbb1.numOfBuried_ = sbb2.numOfBuried_;
        sbb1.numOfHydBuried_ = sbb2.numOfHydBuried_;
      }
      MAXIMIZE(sbb1.transScore_, sbb2.transScore_);
      //MINIMIZE(sbb1.rmsd_, sbb2.rmsd_);
      MAXIMIZE(sbb1.maxPen_, sbb2.maxPen_);
    }
    SuperBB* representive_;
    int count_;
  };

public:
  ClusterSBBs(float distance = 4.0) : size_(0), distance_(distance), numOfRmsdCMs_(0),  numOfRmsdAll_(0){}
  bool addSBB(SuperBB *sbb);
  void reportAll(std::ostream& s);
  bool withinDist(const SuperBB &sbb1, const SuperBB &sbb2, float distance) const;
  bool addSBB(std::vector<Cluster> &vec, SuperBB *sbb);
  void merge();

private:
  //data members
  std::unordered_map<std::string, std::vector<Cluster>> contactHash_; // key = contact map, value = clusters
  std::vector<Cluster*> finalClusters_;
  unsigned int size_;
  float distance_;
  mutable unsigned int numOfRmsdCMs_;
  mutable unsigned int numOfRmsdAll_;
};
#endif
