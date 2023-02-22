#include "GeomScore.h"
#include "Logger.h"

RangeParams::RangeParams(const std::vector<int>& wts) {
  ranges.push_back(-3.6);
  ranges.push_back(-2.2);
  ranges.push_back(-1.0);
  ranges.push_back(1.0);
  weights.insert(weights.begin(), wts.begin(), wts.end());
}

//-------------------------------- GeomScore class -----------------------------

GeomScore::GeomScore(const Surface& lowLevel, MoleculeGrid* grid,
		     const std::vector<int>& wts, float penetration_thr, float ns_thr, float d, float gridMargins) :
  grid_(grid), rangeParams(wts),
  penetrationThr_(penetration_thr), divNsThr(1.0/ns_thr),
  density(d), tree(lowLevel, density, gridMargins)
{}

GeomScore::GeomScore(const Surface& lowLevel,
		     const std::vector<int>& wts, float penetration_thr, float ns_thr, float d, float gridMargins) :
  grid_(NULL), rangeParams(wts),
  penetrationThr_(penetration_thr), divNsThr(1.0/ns_thr),
  density(d), tree(lowLevel, density, gridMargins)
{}

GeomScore::GeomScore(MoleculeGrid* grid,
		     const std::vector<int>& wts, float penetration_thr, float ns_thr, float d, float gridMargins) :
  grid_(grid), rangeParams(wts),
  penetrationThr_(penetration_thr), divNsThr(1.0/ns_thr),
  density(d)
{}

void GeomScore::buildTree() {
  tree.buildTree();
}

void GeomScore::buildTree(const std::vector<bool>& as) {
  tree.buildTree(as);
  //  countActiveSitePoints(as);
}

float GeomScore::maxPenetration(const RigidTrans3& trans) {
  float maxPenetration=MAX_FLOAT;
  Level& highLevel = tree.getLevel(tree.size()-1);
  std::stack<const Node *> st;
  for (unsigned int i=0; i< highLevel.size(); i++)
    st.push(&highLevel[i]);
  while (!st.empty()) {
    const Node *node = st.top();
    st.pop();
    float dist = grid_->getDist(trans*(node->position()));
    if (dist != MAX_FLOAT) {
      if (node->getLevel() == 0) {
	if (dist < maxPenetration)
	  maxPenetration = dist;
	continue;
      }
      if (dist - node->getRadius() <= maxPenetration) {
	const std::vector<const Node *>& lowLevelPointers = node->getChildren();
	for (unsigned int i=0; i<lowLevelPointers.size(); i++)
	  st.push(lowLevelPointers[i]);
      }
    }
  }
  return maxPenetration;
}

bool GeomScore::isPenetrating(const RigidTrans3& trans) {
  Level& highLevel = tree.getLevel(tree.size()-1);
  std::stack<const Node *> st;
  for(unsigned int i=0; i< highLevel.size(); i++) {
    float dist = grid_->getDist(trans*(highLevel[i].position()));
    if(dist < penetrationThr_) {
      return true;
    } else {
      st.push(&highLevel[i]);
    }
  }
  while(!st.empty()) {
    const Node *node = st.top();
    st.pop();
    float dist = grid_->getDist(trans*(node->position()));
    if(dist < penetrationThr_)
      return true;
    if(node->getLevel() != 0) {
      // check if dist - pointRadius is more then thr
      if(dist - node->getRadius() <= penetrationThr_) {
	const std::vector<const Node *>& lowLevelPointers = node->getChildren();
	for(unsigned int i=0; i< lowLevelPointers.size(); i++)
	  st.push(lowLevelPointers[i]);
      }
    }
  }
  return false;
}

bool GeomScore::fastIsPenetrating(const RigidTrans3& trans) {
  Level& highLevel = tree.getLevel(tree.size()-1);
  std::stack<const Node *> st;
  for(unsigned int i=0; i< highLevel.size(); i++) {
    float dist = grid_->getDist(trans*(highLevel[i].position()));
    if(dist < penetrationThr_) {
      return true;
    }
  }
  return false;
}


void GeomScore::getInterface(const RigidTrans3& trans, float low_thr, float high_thr,
                             std::vector<const SurfacePoint*>& interface) {
  Level& highLevel = tree.getLevel(tree.size()-1);
  std::stack<const Node *> st;
  for (unsigned int i=0; i< highLevel.size(); i++)
    st.push(&highLevel[i]);
  while (!st.empty()) {
    const Node *node = st.top();
    st.pop();
    float dist = grid_->getDist(trans*(node->position()));
    if (dist == MAX_FLOAT) continue;
    if (node->getLevel() == 0) {
      if (dist > low_thr && dist < high_thr)
	interface.push_back(node->getPoint());
      continue;
    }
    if (dist - node->getRadius() > low_thr && dist + node->getRadius() < high_thr) {
      node2interface(node, interface);
      continue;
    }
    const std::vector<const Node *>& lowLevelPointers = node->getChildren();
    for (unsigned int i=0; i<lowLevelPointers.size(); i++)
      st.push(lowLevelPointers[i]);
  }
}

void GeomScore::node2interface(const Node* node, std::vector<const SurfacePoint*>& interface) {
  std::stack<const Node *> st;
  st.push(node);
  while (!st.empty()) {
    const Node *node = st.top();
    st.pop();
    if (node->getLevel() == 0) {
      interface.push_back(node->getPoint());
      continue;
    }
    const std::vector<const Node *>& lowLevelPointers = node->getChildren();
    for (unsigned int i=0; i<lowLevelPointers.size(); i++)
      st.push(lowLevelPointers[i]);
  }
}

int GeomScore::score(const RigidTrans3& trans) {
  std::vector<unsigned int> pointsInRanges(rangeParams.weights.size(),0);
  scoreRanges(trans, pointsInRanges);
  //calculate total, based on weights
  int score=0;
  for(unsigned int i=0; i<rangeParams.weights.size(); i++)
    score+= (pointsInRanges[i]*rangeParams.weights[i]);

  if(score <= 0)
    return -1;

  if(((float)pointsInRanges[pointsInRanges.size()-2])/score > divNsThr) {
    return -1;
  }
  return score;
}

int GeomScore::score(const RigidTrans3& trans, const Surface& surface, float& penetration) {
  std::vector<unsigned int> pointsInRanges(rangeParams.weights.size(),0);
  penetration = MAX_FLOAT;
  // count number of points
  for(unsigned int i=0; i<surface.size(); i++) {
    float dist = grid_->getDist(trans * (surface[i].position()));
    if(dist < penetrationThr_)
      return -1;
    if(dist < penetration)
      penetration = dist;
    unsigned int range = rangeParams.findRange(dist);
    pointsInRanges[range]++;
  }

  //calculate total, based on weights
  int score=0;
  for(unsigned int i=0; i<rangeParams.weights.size(); i++)
    score+= (pointsInRanges[i]*rangeParams.weights[i]);

  if(score <= 0)
    return -1;

  if(((float)pointsInRanges[pointsInRanges.size()-2])/score > divNsThr) {
    return -1;
  }
  return score;
}

int GeomScore::score(const RigidTrans3& trans, float& asRatio, int& asScore) {
  std::vector<unsigned int> pointsInRanges(rangeParams.weights.size(),0);
  std::vector<unsigned int> asInRanges(rangeParams.weights.size(),0);
  scoreRanges(trans, pointsInRanges, asInRanges);
  //calculate total, based on weights
  int score=0;
  int interfaceCount = 0;
  int asInterfaceCount = 0;
  asScore=0;
  for(unsigned int i=0; i<rangeParams.weights.size(); i++) {
    score+= (pointsInRanges[i]*rangeParams.weights[i]);
    asScore+= (asInRanges[i]*rangeParams.weights[i]);
    if(i < rangeParams.weights.size() -1 ) {
      interfaceCount+= pointsInRanges[i];
      asInterfaceCount+= asInRanges[i];
    }
  }

  if(interfaceCount == 0) {
    asRatio = 0;
  } else {
    asRatio = ((float)asInterfaceCount)/interfaceCount;
  }

  if(score <= 0)
    return -1;

  if(((float)pointsInRanges[pointsInRanges.size()-2])/score > divNsThr) {
    return -1;
  }

  return score;
}

void GeomScore::scoreRanges(const RigidTrans3& trans, std::vector<unsigned int>& pointsInRanges) {
  Level& highLevel = tree.getLevel(tree.size()-1);
  std::stack<const Node *> st;
  for (unsigned int i=0; i< highLevel.size(); i++)
    st.push(&highLevel[i]);
  while (!st.empty()) {
    const Node *node = st.top();
    st.pop();
    float dist = grid_->getDist(trans*(node->position()));
    if (dist == MAX_FLOAT) { //add to last range
      pointsInRanges[pointsInRanges.size()-1]+=node->getSubtreeSize();
      continue;
    }
    unsigned int range = rangeParams.findRange(dist);
    if (node->getLevel() == 0) {
      pointsInRanges[range]++;
      continue;
    }
    unsigned int range1 = rangeParams.findRange(dist + node->getRadius());
    unsigned int range2 = rangeParams.findRange(dist - node->getRadius());
    if (range == range1 && range == range2) {
      pointsInRanges[range] += node->getSubtreeSize();
      continue;
    }
    const std::vector<const Node *>& lowLevelPointers = node->getChildren();
    for (unsigned int i=0; i<lowLevelPointers.size(); i++)
      st.push(lowLevelPointers[i]);
  }
}

void GeomScore::scoreRanges(const RigidTrans3& trans, std::vector<unsigned int>& pointsInRanges,
                            std::vector<unsigned int>& asInRanges) {
  Level& highLevel = tree.getLevel(tree.size()-1);
  std::stack<const Node *> st;
  for (unsigned int i=0; i< highLevel.size(); i++)
    st.push(&highLevel[i]);
  while (!st.empty()) {
    const Node *node = st.top();
    st.pop();
    Vector3 point = trans*(node->position());
    float dist = grid_->getDist(point);
    if (dist == MAX_FLOAT) { //add to last range
      pointsInRanges[pointsInRanges.size()-1]+=node->getSubtreeSize();
      asInRanges[asInRanges.size()-1]+=node->getAsNum();
      continue;
    }
    unsigned int range = rangeParams.findRange(dist);
    if (node->getLevel() == 0) {
      pointsInRanges[range]++;
      asInRanges[range]+=node->getAsNum();
      continue;
    }
    unsigned int range1 = rangeParams.findRange(dist + node->getRadius());
    unsigned int range2 = rangeParams.findRange(dist - node->getRadius());
    if (range == range1 && range == range2) {
      pointsInRanges[range]+=node->getSubtreeSize();
      asInRanges[range]+=node->getAsNum();
      continue;
    }
    const std::vector<const Node *>& lowLevelPointers = node->getChildren();
    for (unsigned int i=0; i< lowLevelPointers.size(); i++)
      st.push(lowLevelPointers[i]);
  }
}

ScoreData GeomScore::fullScore(const RigidTrans3& trans, bool as) {
  ScoreData scoreData;
  std::vector<unsigned int> pointsInRanges(rangeParams.weights.size(),0);
  if (as) {
    std::vector<unsigned int> asInRanges(rangeParams.weights.size(),0);
    scoreRanges(trans, pointsInRanges, asInRanges);
    for (unsigned int i=0; i<rangeParams.weights.size()-1; i++) {
      scoreData.score+= (pointsInRanges[i] * rangeParams.weights[i]);
      scoreData.asScore+= (asInRanges[i] * rangeParams.weights[i]);
      scoreData.interfaceArea+= pointsInRanges[i];
    }
    asInRanges.clear();
  } else {
    scoreRanges(trans, pointsInRanges);
    for (unsigned int i=0; i<rangeParams.weights.size()-1; i++) {
      scoreData.score+= (pointsInRanges[i] * rangeParams.weights[i]);
      scoreData.interfaceArea+= pointsInRanges[i];
    }
  }
  scoreData.interfaceArea/=density;
  scoreData.maxPenetration = maxPenetration(trans);
  pointsInRanges.clear();
  return scoreData;
}

float GeomScore::computePropensity(const RigidTrans3& trans) {
  std::vector<unsigned int> pointsInRanges(rangeParams.weights.size(),0);
  std::vector<unsigned int> asInRanges(rangeParams.weights.size(),0);
  scoreRanges(trans, pointsInRanges, asInRanges);
  //calculate total, based on weights
  int score=0;
  int asScore=0;
  for (unsigned int i=0; i<rangeParams.weights.size()-1; i++) {
    score+= pointsInRanges[i];
    asScore+= asInRanges[i];
  }
  float prop = ((float)asScore/score)/((float)tree.getAsPointsNumber()/tree.getLevel(0).size());
  return prop;
}

int GeomScore::refineTrans(const RigidTrans3& trans, RigidTrans3& newTrans) {
/*
    float Tmag = 0.1;
    float Rmag =  0.05;


    // compute ligand interface points
    std::vector<const SurfacePoint*> interface;
    getInterface(trans, -2.0, 1.5, interface);  //add thr - which??

    Logger::debugMessage() << "interface size: " << interface.size() << " trans " << trans << endl;

    std::vector<Vector3> interfacePoints(interface.size());
    std::vector<Vector3> virtualPoints(interface.size());
    MyMatch matchList;

    Vector3 centerOfRotation;
    for (unsigned int i=0; i<interface.size(); i++) {
        interfacePoints[i] = interface[i]->position();
        centerOfRotation += interfacePoints[i];
        matchList.add(i,i);
    }
    centerOfRotation /= interface.size();

    RigidTrans3 bestTrans;
    RigidTrans3 currTrans = trans;
    int currScore=1000000;
    int oldScore = score(trans);
    MMonteCarlo::init(oldScore);
    newTrans = trans;
    int nums = 0;
    while (nums < 1000) {

        currTrans = newTrans;
        nums++;
        //Logger::debugMessage() << "oldScore: " << oldScore;
        for (unsigned int i=0; i<interface.size(); i++) {
            Vector3 point = newTrans * interfacePoints[i];
            float dist = grid_->getDist(point);
            virtualPoints[i] = point + dist*(interface[i]->normal());
            if (dist == 0.00)
                matchList.setScore(i,10.0);
            else
                matchList.setScore(i,1/fabs(dist));
        }
        matchList.calculateBestFit(virtualPoints, interfacePoints);
        newTrans = matchList.rigidTrans();
        //oldScore = currScore;
        currScore = score(newTrans);
        MMonteCarlo::MMC_acceptance mmc_accept = MMonteCarlo::step(currScore);

        if (mmc_accept == MMonteCarlo::ACCEPT_LOWEST) {
            bestTrans = newTrans;
        }

        if (mmc_accept == MMonteCarlo::REJECT) {


            float delta_TR[6] = {0.0};
            delta_TR[0] = MMonteCarlo::getGaussianRandomNumber() * Tmag;
            delta_TR[1] = MMonteCarlo::getGaussianRandomNumber() * Tmag;
            delta_TR[2] = MMonteCarlo::getGaussianRandomNumber() * Tmag;
            delta_TR[3] = MMonteCarlo::getGaussianRandomNumber() * Rmag;
            delta_TR[4] = MMonteCarlo::getGaussianRandomNumber() * Rmag;
            delta_TR[5] = MMonteCarlo::getGaussianRandomNumber() * Rmag;

            Vector3 rotVec(delta_TR[3],delta_TR[4],delta_TR[5]);
            Vector3 transVec(delta_TR[0],delta_TR[1],delta_TR[2]);

            RigidTrans3 delta_rt3(rotVec,transVec);

            const Matrix3 &delta_mat3 = delta_rt3.rotation();
            RigidTrans3 delta_rt3_COR(delta_mat3,centerOfRotation-delta_mat3*centerOfRotation+delta_rt3.translation());
            newTrans = delta_rt3_COR*currTrans;

        }

        //Logger::debugMessage() << " RMSD: " << matchList.rmsd() << " newTrans " << newTrans << " newScore " << currScore << endl;
    }

    newTrans = bestTrans;

    //Logger::debugMessage() << "oldScore: " << oldScore << "  newScore " << currScore << endl;
    return currScore;*/
return 0;
}
