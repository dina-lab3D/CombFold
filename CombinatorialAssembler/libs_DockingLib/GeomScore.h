#ifndef GEOM_SCORE_H
#define GEOM_SCORE_H

#include <stack>
#include <utility>
#include <vector>

#include <Surface.h>
#include <MoleculeGrid.h>
#include "MultiResolution.h"

#include <Match.h>
#include <Particle.h>
//#include "MMonteCarlo.h"

class MyMatch : public Match {
 public:
  void setScore(unsigned int pairIndex, float score) { pairs[pairIndex].score = score; }
};

class RangeParams {
 public:
  RangeParams(const std::vector<int>& wts);
  unsigned int findRange(float dist) const {
    for(unsigned int i=0; i< ranges.size(); i++)
      if(dist <= ranges[i]) return i;
    return ranges.size();
  }
 public:
  std::vector<float> ranges;
  std::vector<int> weights;
};

class ScoreData {
 public:
  ScoreData() { score=0; asScore=0; interfaceArea=0;}
  int score;
  int asScore;
  float interfaceArea;
  float maxPenetration;
 public:
  //  int getSurfaceRange() {return pointsInRanges[3];}
  friend std::ostream& operator<<(std::ostream& s, const ScoreData& data) {
    s.width(5); s << data.score << "|";
    s.width(5); s << data.asScore << "|";
    s.width(5); s << data.maxPenetration << "|";
    s.width(5); s << data.interfaceArea << "|";
    return s;
  }
};

class ScorePair : public std::pair<ScoreData, ScoreData> {
 public:
  ScorePair(ScoreData f, ScoreData s, RigidTrans3 trans) : std::pair<ScoreData, ScoreData> (f,s) {
    trans_ = trans;
    score_ = f.score + s.score - abs(f.score - s.score);
  }
  int score() const { return score_; }
  const RigidTrans3& rigidTrans() const { return trans_; }
  friend std::ostream& operator<<(std::ostream& s, const ScorePair& p) {
    s << p.first << "|";
    s << std::endl << "           ";
    s << p.second << "|";
    return s;
  }
 protected:
  RigidTrans3 trans_;
  int score_;
};

class GeomScore {
 public:
  typedef MultiResolution<SurfacePoint>::Node Node;
  typedef MultiResolution<SurfacePoint>::Level Level;

  GeomScore(const Surface& lowLevel, MoleculeGrid* grid,
	    const std::vector<int>& wts, float penetration_thr, float ns_thr, float density=10.0, float gridMargins=7.0);

  GeomScore(const Surface& lowLevel,
	    const std::vector<int>& wts, float penetration_thr, float ns_thr, float density=10.0, float gridMargins=7.0);

  GeomScore(MoleculeGrid* grid,
	    const std::vector<int>& wts, float penetration_thr, float ns_thr=0.5, float density=10.0, float gridMargins=7.0);

  void setGrid(const MoleculeGrid* grid) { grid_ = grid; }

  void setSurface(const Surface& lowLevel) {
    tree.setPointSet(lowLevel);
    buildTree();
  }

  void buildTree();
  void buildTree(const std::vector<bool>& as);

  //// Penetration checks
  float maxPenetration(const RigidTrans3& trans);
  bool isPenetrating(const RigidTrans3& trans);

  //// fast check, but can return false for penetrating trans
  bool fastIsPenetrating(const RigidTrans3& trans);

  void getInterface(const RigidTrans3& trans, float low_thr, float high_thr, std::vector<const SurfacePoint*>& interface);

  //// Scoring functions
  int score(const RigidTrans3& trans);
  int score(const RigidTrans3& trans, float& asRatio, int& asScore);
  int score(const RigidTrans3& trans, const Surface& surface, float& penetration);
  ScoreData fullScore(const RigidTrans3& trans, bool as=false);

  int refineTrans(const RigidTrans3& trans, RigidTrans3& newTrans);

  void printTree(std::ofstream& outFile);
  float computePropensity(const RigidTrans3& trans);
 private:
  void scoreRanges(const RigidTrans3& trans, std::vector<unsigned int>& pointsInRanges);
  void scoreRanges(const RigidTrans3& trans, std::vector<unsigned int>& pointsInRanges,
		   std::vector<unsigned int>& asInRanges);
  void node2interface(const Node* node, std::vector<const SurfacePoint*>& interface);
  void countActiveSitePoints(const std::vector<bool>& as);
 private:
  const MoleculeGrid* grid_;
  const RangeParams rangeParams;
  float penetrationThr_;
  float divNsThr; // 1/nsThr
  float density;
  MultiResolution<SurfacePoint> tree;
};

#endif
