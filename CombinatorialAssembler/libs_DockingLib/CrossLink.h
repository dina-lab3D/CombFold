#ifndef CROSS_LINK_H
#define CROSS_LINK_H

#include <fstream>
#include <iostream>
#include <vector>

class CrossLink {
 public:
  CrossLink()
      : residueNumber1_(0), chainId1_(" "),
        residueNumber2_(0), chainId2_(" "),
        minDistance_(0.0), maxDistance_(0.0), weight_(1.0) {}

  CrossLink(int residueNumber1, std::string chainId1,
            int residueNumber2, std::string chainId2,
            float minDistance, float maxDistance, float weight=1.0)
    : residueNumber1_(residueNumber1), residueSequenceID1_(std::to_string(residueNumber1)), chainId1_(chainId1),
      residueNumber2_(residueNumber2), residueSequenceID2_(std::to_string(residueNumber2)), chainId2_(chainId2),
        minDistance_(minDistance), maxDistance_(maxDistance), weight_(weight) {}


  int getResidue1() const { return residueNumber1_; }
  int getResidue2() const { return residueNumber2_; }

  std::string getResidueSequenceID1() const { return residueSequenceID1_; }
  std::string getResidueSequenceID2() const { return residueSequenceID2_; }

  std::string getChain1() const { return chainId1_; }
  std::string getChain2() const { return chainId2_; }

  float getMinDistance() const { return minDistance_; }
  float getMaxDistance() const { return maxDistance_; }
  float getWeight() const { return weight_; }

  bool operator==(const CrossLink& cl) const;

  bool isInList(const std::vector<CrossLink>& crossLinks) const;

  CrossLink getReversed() const;

  friend std::ostream& operator<<(std::ostream& q, const CrossLink& cl);
  friend std::istream& operator>>(std::istream& s, CrossLink& cl);

 protected:
  int residueNumber1_;
  std::string residueSequenceID1_; // to support numbering of type 99A
  std::string chainId1_;
  int residueNumber2_;
  std::string residueSequenceID2_;
  std::string chainId2_;
  float minDistance_;
  float maxDistance_;
  float weight_;
};


int readCrossLinkFile(const std::string& fileName,
                      std::vector<CrossLink>& crossLinks,
                      bool addReverse = false);


void writeCrossLinkFile(const std::string& fileName,
                        const std::vector<CrossLink>& crossLinks);

void writeXLAnalizerFile(const std::string& fileName,
                         const std::vector<CrossLink>& crossLinks);

#endif
