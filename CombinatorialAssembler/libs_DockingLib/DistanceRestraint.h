#ifndef DISTANCE_RESTRAINT_H
#define DISTANCE_RESTRAINT_H

#include <Vector3.h>
#include <RigidTrans3.h>

#include <vector>
#include <limits>

class DistanceRestraint {
public:
  DistanceRestraint(const std::vector<Vector3>& p1,
                    const std::vector<Vector3>& p2,
                    float maxdist, float mindist = 0.0, float weight=1.0)
    : p1_(p1), p2_(p2), mindist_(mindist), maxdist_(maxdist),
      mindist2_(mindist*mindist), maxdist2_(maxdist*maxdist),
      weight_(weight) {}

  DistanceRestraint(const std::vector<Vector3>& p1,
                    const std::vector<Vector3>& p2,
                    const std::vector<std::pair<std::string, std::string>>& p1Info,
                    const std::vector<std::pair<std::string, std::string>>& p2Info,
                    float maxdist, float mindist = 0.0, float weight=1.0)
    : DistanceRestraint(p1, p2, maxdist, mindist, weight) {
    p1Info_ = p1Info;
    p2Info_ = p2Info;
  }

  DistanceRestraint(const Vector3& p1, const Vector3& p2,
                    float maxdist, float mindist = 0.0, float weight=1.0)
    :  mindist_(mindist), maxdist_(maxdist),
      mindist2_(mindist*mindist), maxdist2_(maxdist*maxdist),
      weight_(weight)
  {
    p1_.push_back(p1);
    p2_.push_back(p2);
  }

  DistanceRestraint() = default;

  float getWeight() const { return weight_; }
  float getMinDistance() const { return mindist_; }
  float getMaxDistance() const { return maxdist_; }
  const std::vector<Vector3>& getPoints1() const { return p1_; }
  const std::vector<Vector3>& getPoints2() const { return p2_; }

  bool isViolated() const {
    // calculate the distance and compare to mindist_ and maxdist_
    for(unsigned int i = 0; i < p1_.size(); i++) {
      for(unsigned int j = 0; j < p2_.size(); j++) {
        float dist2 = p1_[i].dist2(p2_[j]);
        if(dist2 <= maxdist2_) return false; // found one in range
      }
    }
    return true;
  }

  bool isViolated(const std::vector<RigidTrans3>& trans) const {
    // calculate the distance and compare to mindist_ and maxdist_
    for(unsigned int i = 0; i < p1_.size(); i++) {
      for(unsigned int j = 0; j < p2_.size(); j++) {
        float dist2 = p1_[i].dist2(trans[j]*p2_[j]);
        if(dist2 <= maxdist2_) return false; // found one in range
      }
    }
    return true;
  }

  bool isViolated(const std::vector<RigidTrans3>& trans1,
                  const std::vector<RigidTrans3>& trans2) const {
    // calculate the distance and compare to mindist_ and maxdist_
    for(unsigned int i = 0; i < p1_.size(); i++) {
      Vector3 transP1 = trans1[i]*p1_[i];
      for(unsigned int j = 0; j < p2_.size(); j++) {
        float dist2 = transP1.dist2(trans2[j]*p2_[j]);
        if(dist2 <= maxdist2_) return false; // found one in range
      }
    }
    return true;
  }

  // backward compatability
  bool isSatisfied() const { return !isViolated(); }
  bool isSatisfied(const RigidTrans3& trans) const {
    std::vector<RigidTrans3> t(p2_.size(), trans);
    return !isViolated(t);
  }
  bool isSatisfied(const RigidTrans3& trans1,
                   const RigidTrans3& trans2) const {
    std::vector<RigidTrans3> t1(p1_.size(), trans1);
    std::vector<RigidTrans3> t2(p2_.size(), trans2);
    return !isViolated(t1, t2);
  }

  //calculate violationDistance
  float violationDistance() const {
    float bestdist = distance();
    if(bestdist > maxdist_) return bestdist - maxdist_;
    return 0.0;
  }

  float violationDistance(const std::vector<RigidTrans3>& trans) const {
    float bestdist = distance(trans);
    if(bestdist > maxdist_) return bestdist - maxdist_;
    return 0.0;
  }

  float violationDistance(const std::vector<RigidTrans3>& trans1,
                          const std::vector<RigidTrans3>& trans2) const {
    float bestdist = distance(trans1, trans2);
    if(bestdist > maxdist_) return bestdist - maxdist_;
    return 0.0;
  }

  float distance2() const {
    float bestdist2 = std::numeric_limits<float>::max();
    for(unsigned int i = 0; i < p1_.size(); i++) {
      for(unsigned int j = 0; j < p2_.size(); j++) {
        float dist2 = p1_[i].dist2(p2_[j]);
        if(dist2 < bestdist2) bestdist2 = dist2;
      }
    }
    return bestdist2;
  }

  float distance() const { return sqrt(distance2()); }

  float distance2(const std::vector<RigidTrans3>& trans) const {
    float bestdist2 = std::numeric_limits<float>::max();
    for(unsigned int i = 0; i < p1_.size(); i++) {
      for(unsigned int j = 0; j < p2_.size(); j++) {
        float dist2 = p1_[i].dist2(trans[j]*p2_[j]);
        if(dist2 < bestdist2) bestdist2 = dist2;
      }
    }
    return bestdist2;
  }

  float distance(const std::vector<RigidTrans3>& trans) const {
    return sqrt(distance2(trans));
  }

  float distance(const RigidTrans3& trans) const {
    std::vector<RigidTrans3> t(1, trans);
    return distance(t);
  }

  float distance2(const std::vector<RigidTrans3>& trans1,
                  const std::vector<RigidTrans3>& trans2) const {
    float bestdist2 = std::numeric_limits<float>::max();
    for(unsigned int i = 0; i < p1_.size(); i++) {
      Vector3 transP1 = trans1[i]*p1_[i];
      for(unsigned int j = 0; j < p2_.size(); j++) {
        float dist2 = transP1.dist2(trans2[j]*p2_[j]);
        if(dist2 < bestdist2) bestdist2 = dist2;
      }
    }
    return bestdist2;
  }

  float distance(const std::vector<RigidTrans3>& trans1,
                 const std::vector<RigidTrans3>& trans2) const {
    return sqrt(distance2(trans1, trans2));
  }

  // TODO: add trans versions
  float minDistanceIndices(unsigned int& index1, unsigned int& index2) const;

  std::string getChimeraXPseudoBond() const;

  bool operator < (const DistanceRestraint& d) const { return (maxdist_ < d.maxdist_); }

  friend std::ostream& operator<<(std::ostream& s, const DistanceRestraint& d) {
    // TODO
    return s << d.p1_[0] << ' ' << d.p2_[0] << ' ' << d.mindist_ <<' ' << d.maxdist_;
  }

private:
  std::vector<Vector3> p1_; // end point1, vector for ambiguity
  std::vector<Vector3> p2_; // end point2, vector for ambiguity
  std::vector<std::pair<std::string, std::string>> p1Info_, p2Info_; // residueSequenceID and chain for p1 and p2 endpoints
  float mindist_; // minimal distance threshold
  float maxdist_; // maximal distance threshold
  float mindist2_; // minimal distance threshold (squared)
  float maxdist2_; // maximal distance threshold (squared)
  float weight_;
};

#endif /* DISTANCE_RESTRAINT_H */
