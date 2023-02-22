#include "DistanceRestraint.h"
#include <boost/algorithm/string/trim.hpp>

float DistanceRestraint::minDistanceIndices(unsigned int& index1, unsigned int& index2) const {
  float bestdist = std::numeric_limits<float>::max();
  for(unsigned int i = 0; i < p1_.size(); i++) {
    for(unsigned int j = 0; j < p2_.size(); j++) {
      float dist = p1_[i].dist(p2_[j]);
      if(dist < bestdist) {
        bestdist = dist;
        index1 = i;
        index2 = j;
      }
    }
  }
  return bestdist;
}

std::string DistanceRestraint::getChimeraXPseudoBond() const {
  unsigned int best_i = 0, best_j = 0;
  float bestdist = minDistanceIndices(best_i, best_j);

  if(p1Info_.size() == p1_.size() && p2Info_.size() == p2_.size()) {
    // /E:157@ca /B:58@ca  red
    std::string res1 = p1Info_[best_i].first;
    boost::algorithm::trim(res1);
    std::string res2 = p2Info_[best_j].first;
    boost::algorithm::trim(res2);

    // skip self cross link
    if(res1 == res2 && p1Info_[best_i].second == p2Info_[best_j].second)
      return std::string();

    std::string ret = "/" + p1Info_[best_i].second + ":" + res1 + "@ca /" +
      p2Info_[best_j].second + ":" + res2 + "@ca";
    if(bestdist > maxdist_) ret += " red";
    return ret;
  }
  return std::string();
}
