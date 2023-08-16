//
// Created by dina on 8/1/18.
//

#ifndef BB_CONTAINER_H
#define BB_CONTAINER_H

#include "BB.h"
#include <memory>

class BBContainer {
  public:
    // Constructor
    BBContainer(std::string SUFileName, std::string chemLibFileName, float minTempFactor);

    // Group: access
    std::shared_ptr<const BB> getBB(unsigned int bbIndex) const { return bbs_[bbIndex]; }
    const std::vector<std::shared_ptr<const BB>> &getBBs() const { return bbs_; }
    unsigned int getBBsNumber() const { return numOfBBs_; }

    void readTransformationFiles(std::string transFilePrefix, unsigned int transNumToRead);

  private:
    int readSUFile(const std::string SUFileName);
    bool readTrans(TransformationAndScore &trans, std::ifstream &transFile);

  private:
    // PDB filenames
    std::vector<std::string> pdbs_;

    // BBs
    std::vector<std::shared_ptr<const BB>> bbs_;

    // total bbs num
    unsigned int numOfBBs_;

    // BBs can be grouped according to a number provided by the user in SUlist
    // BBs in the same group will be assembled first
    // if no number is given, all the BBs are considered a group
    std::vector<int> groupIDs_;
};

#endif // BBCONTAINER_H
