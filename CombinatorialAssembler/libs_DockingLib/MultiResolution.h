#ifndef MULT_RES_H
#define MULT_RES_H

#include <Vector3.h>
#include <Atom.h>
#include <GeomHash.h>
#include <numerics.h>

template<class PointT>
class MultiResolution {
 public:

  //// TreeNode
  class Node {

   public:
    // GROUP: Constructors
    Node(const PointT* point, short level=0, float lowLevelRadius=0.0,
	 unsigned int subtreeSize=1, unsigned int asNum=0):
      point_(point), level_(level), lowLevelRadius_(lowLevelRadius), subtreeSize_(subtreeSize), asNum_(asNum) {}

    // GROUP: Modifiers
    //// update asNum
    void incrementAsNum(unsigned int num) { asNum_ += num; }

    //// add children to node
    void addChildren(const std::vector<const Node *>& children) {
      lowLevelPointers_.insert(lowLevelPointers_.begin(), children.begin(), children.end());
    }

    // GROUP: Queries
    //// return point pointer
    const PointT* getPoint() const { return point_; }

    //// return point position
    Vector3 position() const { return point_->position(); }

    //// get node level
    short getLevel() const { return level_; }

    //// get radius
    float getRadius() const { return lowLevelRadius_; }

    //// get subtree size
    unsigned int getSubtreeSize() const { return subtreeSize_; }

    //// get number of active nodes in the subtree
    unsigned int getAsNum() const { return asNum_; }

    //// get pointers to node children
    const std::vector<const Node *>& getChildren() const { return lowLevelPointers_; }

   protected:
    // point
    const PointT* point_;
    // level
    short level_;
    // radius of lower level points
    float lowLevelRadius_;
    // number of subtree nodes
    unsigned int subtreeSize_;
    // number of active site nodes
    unsigned int asNum_;
    // pointers to lower level
    std::vector<const Node *> lowLevelPointers_;

  };

  typedef std::vector<Node> Level;

  // GROUP: Contructors
  //// The initial bin size is computed based on density and number of tree
  // levels is computed based on the size of the grid margins
  MultiResolution(const std::vector<PointT>& pointSet, float density=10.0, float gridMargins=6.0);
  MultiResolution(float density=10.0, float gridMargins=6.0);

  // GROUP: Modifiers
  void setPointSet(const std::vector<PointT>& pointSet) { pointSet_ = pointSet; }
  void buildTree();
  void buildTree(const std::vector<bool>& as);
  void countActiveSitePoints(const std::vector<bool>& as);

  // GROUP: Queries
  Level& getLevel(unsigned int num) { return tree_[num]; }
  unsigned int size() const { return tree_.size(); }
  int getAsPointsNumber() const { return asPointsNumber_; }
  void printTree(std::ofstream& outFile);

 private:
  void createHighLevel(const int lowLevelIndex, float radius);
  void findMaxRadiusForPoint(const Vector3& point, const Node* currNode, float& currMaxRadius2);

 private:
  std::vector<PointT> pointSet_;
  std::vector<Level> tree_;
  int asPointsNumber_;
  float binSize_;
};

template<class PointT>
MultiResolution<PointT>::MultiResolution(const std::vector<PointT>& pointSet, float density, float gridMargins) :
  pointSet_(pointSet) {
  unsigned int POINTS_PER_CUBE = 5;
  binSize_ = POINTS_PER_CUBE/density;
  Level l;
  for(float cs=binSize_; cs <= gridMargins*2; cs*=2) { //-1
    tree_.push_back(l);
  }
  //  cout << "Tree size: " << tree_.size() << endl;
}

template<class PointT>
MultiResolution<PointT>::MultiResolution(float density, float gridMargins) {
  unsigned int POINTS_PER_CUBE = 5;
  binSize_ = POINTS_PER_CUBE/density;
  Level l;
  for(float cs=binSize_; cs <= gridMargins*2; cs*=2) { //-1
    tree_.push_back(l);
  }
  // cout << "Tree size: " << tree_.size() << endl;
}


template<class PointT>
void MultiResolution<PointT>::createHighLevel(const int lowLevelIndex, float binSize) {
  const Level& lowLevel = tree_[lowLevelIndex];
  // insert to ghash
  GeomHash<Vector3, const Node* > gHash(3, binSize);
  for(unsigned int i=0; i< lowLevel.size(); i++)
    gHash.insert(lowLevel[i].position(), &lowLevel[i]);

  // init Level
  short levelIndex = lowLevelIndex+1;
  Level& highLevel = tree_[levelIndex];

  // compute high resolution - iterate over each Bucket
  typename GeomHash<Vector3, const Node*>::BucketsPointerList *bucketList = gHash.getBuckets();
  typename GeomHash<Vector3, const Node*>::BucketsPointerList::const_iterator bIter, bEndIter = bucketList->end();
  for(bIter = bucketList->begin(); bIter != bEndIter; bIter++) {
    const std::vector<const Node*>& currBucket = **bIter;
    // compute average bucket point
    // compute subtree size
    Vector3 average;
    unsigned int subTreeSize=0;
    typename std::vector<const Node*>::const_iterator currBucketIter, bucketEndIter = currBucket.end();
    for(currBucketIter = currBucket.begin(); currBucketIter != bucketEndIter; currBucketIter++) {
      average+= (*currBucketIter)->position();
      subTreeSize+=(*currBucketIter)->getSubtreeSize();
    }
    average/=currBucket.size();

    // find point closest to average
    float minDist2 = MAX_FLOAT;
    const Node* selectedNode=NULL;
    for(currBucketIter = currBucket.begin(); currBucketIter != bucketEndIter; currBucketIter++) {
      float currDist2 = average.dist2((*currBucketIter)->position());
      if(currDist2 < minDist2) {
	selectedNode = *currBucketIter;
	minDist2 = currDist2;
      }
    }

    // find radius of pts == maxDist from average
    float maxDist2 = MIN_FLOAT;
    for(currBucketIter = currBucket.begin(); currBucketIter != bucketEndIter; currBucketIter++)
      findMaxRadiusForPoint(selectedNode->position(), *currBucketIter, maxDist2);

    // create new node
    Node newNode(selectedNode->getPoint(), levelIndex, sqrt(maxDist2), subTreeSize);
    newNode.addChildren(currBucket);
    highLevel.push_back(newNode);
  }
  delete bucketList;
  // cout << ": " << highLevel.size() << "  ";
}

template<class PointT>
void MultiResolution<PointT>::findMaxRadiusForPoint(const Vector3& point, const Node* currNode, float& currMaxRadius2) {
  int levelIndex = currNode->getLevel();
  float dist2 = currNode->position().dist2(point);
  if(levelIndex==0) {
    if(dist2 > currMaxRadius2)
      currMaxRadius2=dist2;
    return;
  }
  // check if distance to (point+radius)^2 < currMaxRadius2
  float radius = currNode->getRadius();
  float dist = sqrt(dist2);
  if(sqr(radius+dist) < currMaxRadius2) return;
  const std::vector<const Node*>& lowLevelPointersForPoint = currNode->getChildren();
  typename std::vector<const Node*>::const_iterator iter, endIter=lowLevelPointersForPoint.end();
  if(levelIndex==1) {
    for(iter=lowLevelPointersForPoint.begin(); iter!=endIter; iter++) {
      float dist2=point.dist2((*iter)->position());
      if(dist2 > currMaxRadius2)
	currMaxRadius2=dist2;
    }
    return;
  }
  if(levelIndex > 1) {
    for(iter=lowLevelPointersForPoint.begin(); iter!=endIter; iter++) {
      findMaxRadiusForPoint(point, *iter, currMaxRadius2);
    }
    return;
  }
  return;
}

template<class PointT>
void MultiResolution<PointT>::buildTree() {
  // initiate lowest level
  Level& lowLevel = tree_[0];
  for(unsigned int i=0; i<pointSet_.size(); i++) {
    Node node(&pointSet_[i], 0, 0, 1);
    lowLevel.push_back(node);
  }

  // create the rest of the tree
  // cout << "start building levels ";
  float radius = binSize_;
  for(unsigned int i=0; i< tree_.size()-1; i++) {
    // cout << i+1;
    createHighLevel(i, radius);
    radius*=2;
  }
  // cout << "done " << endl;
}

template<class PointT>
void MultiResolution<PointT>::buildTree(const std::vector<bool>& as) {
  buildTree();
  countActiveSitePoints(as);
}

template<class PointT>
void MultiResolution<PointT>::countActiveSitePoints(const std::vector<bool>& as) {
  Level& lowLevel = tree_[0];
  if(as.size() != lowLevel.size()) {
    std::cerr << "Error in active site array" << std::endl;
    return;
  }
  // mark low level leaves
  asPointsNumber_=0;
  for(unsigned int i=0; i<lowLevel.size(); i++) {
    lowLevel[i].incrementAsNum(as[i]);
    if(as[i]==1) asPointsNumber_++;
  }
  // go up the tree and mark nodes
  for(unsigned int levelIndex=1; levelIndex<tree_.size(); levelIndex++) {
    Level& currLevel = tree_[levelIndex];
    for(unsigned int nodeIndex=0; nodeIndex<currLevel.size(); nodeIndex++) {
      Node& currNode = currLevel[nodeIndex];
      const std::vector<const Node *>& lowLevelPointers = currNode.getChildren();
      for(unsigned int lowLevelIndex=0; lowLevelIndex<lowLevelPointers.size(); lowLevelIndex++)
	currNode.incrementAsNum(lowLevelPointers[lowLevelIndex]->getAsNum());
    }
  }
}

template<class PointT>
void MultiResolution<PointT>::printTree(std::ofstream& outFile)
{
  int count=0;
  for(unsigned int levelIndex=0; levelIndex<tree_.size(); levelIndex++) {
    Level& currLevel = tree_[levelIndex];
    for(unsigned int nodeIndex=0; nodeIndex<currLevel.size(); nodeIndex++) {
      Node& currNode = currLevel[nodeIndex];
      if(currNode.asNum_ > 0) {
	Atom atom(currNode.position(), 'A', count++, levelIndex, "INT", 'X');
	outFile << atom << std::endl;
      } else {
	Atom atom(currNode.position(), 'A', count++, levelIndex, "PRB", 'X');
	outFile << atom << std::endl;
      }
    }
  }
}

#endif
