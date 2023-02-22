// Created by dina on 8/1/18.

#include "BBContainer.h"

#include <boost/algorithm/string.hpp>

namespace {
    std::string trim_extension(const std::string file_name) {
        if (file_name[file_name.size() - 4] == '.')
            return file_name.substr(0, file_name.size() - 4);
        return file_name;
    }
}

BBContainer::BBContainer(const std::string SUFileName, float penetrationThr, std::string chemLibFileName) {
  readSUFile(SUFileName);

  // init parameters
  BB::scoreParams_.weights[0] = 0;
  BB::scoreParams_.weights[1] = 0;
  BB::scoreParams_.weights[2] = 0;
  BB::scoreParams_.weights[3] = 1;
  BB::scoreParams_.weights[4] = 0;
  BB::scoreParams_.max_penetration = BB::penetrationThreshold_ = penetrationThr;
  BB::gridResolution_ = 0.5; //p.getFloat("GRID_RESOLUTION");
  BB::gridMargins_ = 5.0; //p.getFloat("GRID_MARGINS");

  // prepare ChemLib
  ChemLib chemLib(chemLibFileName);

  // read the building blocks
  bbs_.reserve(numOfBBs_);
  for (unsigned int i = 0; i < numOfBBs_; i++) {
    bbs_.push_back(std::make_shared<BB>(i, pdbs_[i], groupIDs_[i], groupID2size_[groupIDs_[i]], chemLib));
  }
}

void BBContainer::readTransformationFiles(std::string transFilePrefix, unsigned int transNumToRead) {
  // init transformations vector
  for (unsigned int i = 0; i < numOfBBs_; i++) {
    bbs_[i]->initTrans(numOfBBs_, {});
  }

  // read the files
  for (size_t i = 0; i < numOfBBs_; i++) {
    for (size_t j = 0; j < numOfBBs_; j++) {
      if (i == j) continue;
      std::string inFileName =
        transFilePrefix + trim_extension(pdbs_[i]) + "_plus_" + trim_extension(pdbs_[j]);
      std::ifstream inS(inFileName);
      std::cerr << "Opening transformation file " << inFileName << std::endl;
      if (!inS) {
        std::cerr << "Problem opening transformation file " << inFileName << std::endl;
        continue;
      }

      unsigned int counter = 0;
      while (!inS.eof() && counter < transNumToRead) {
        TransformationAndScore *t1 = new TransformationAndScore();
        if (readTrans(*t1, inS)) {
          std::shared_ptr<TransformationAndScore> t2 = std::make_shared<TransformationAndScore>(*t1);
          t2->refFrame_ = !t2->refFrame_;
          bbs_[i]->putTransWith(j, std::make_shared<TransformationAndScore>(*t1), {});
          bbs_[j]->putTransWith(i, t2, {});
          counter++;
        }
      }
      std::cerr << counter << " transforms were read from file " << inFileName << std::endl;
      inS.close();
    }
  }

  computeRanges(1);
}

int BBContainer::readSUFile(const std::string SUFileName) {
  numOfBBs_ = 0;
  std::ifstream SUFile(SUFileName);
  if (!SUFile) {
    std::cerr << "Can't open SU file" << SUFileName << std::endl;
    exit(1);
  }
  while (!SUFile.eof()) {
    std::string line;
    getline(SUFile, line);
    boost::trim(line);
    if (line.length() > 0) {
      std::vector<std::string> split_results;
      boost::split(split_results, line, boost::is_any_of(" "),
                   boost::token_compress_on);
      std::string pdbName(split_results[0]);
      std::cout << "PDBname " << pdbName << ":" << line <<  std::endl;
      numOfBBs_++;
      pdbs_.push_back(pdbName);
      // group assignment
      int groupID = 1;
      if (split_results.size() == 2) {
        groupID = stoi(split_results[1]);
      }
      groupIDs_.push_back(groupID);
      if(groupID2size_.find(groupID) == groupID2size_.end()) {
        groupID2size_[groupID] = 1;
      } else {
        groupID2size_[groupID]++;
      }
    }
  }
  return numOfBBs_;
}

bool BBContainer::readTrans(TransformationAndScore &trans, std::ifstream &transFile) {

  std::string line;
  if (!transFile.eof()) {
    getline(transFile, line);
    boost::trim(line);  // remove all spaces
    // skip comments
    if (line[0] == '#' || line[0] == '\0') return false;
    std::vector<std::string> split_results;
    boost::split(split_results, line, boost::is_any_of(" "),
                 boost::token_compress_on);
    //std::cerr << split_results.size() << std::endl;
    if (split_results.size() < 7) return false;

    // extract trans
    unsigned int trans_index = split_results.size() - 6;
    RigidTrans3 tr(Vector3(std::stof(split_results[trans_index].c_str()),
                           std::stof(split_results[trans_index+1].c_str()),
                           std::stof(split_results[trans_index+2].c_str())),
                   Vector3(std::stof(split_results[trans_index+3].c_str()),
                           std::stof(split_results[trans_index+4].c_str()),
                           std::stof(split_results[trans_index+5].c_str())));
    trans.refFrame_ = tr;
    trans.score_.totalScore_ = 1.0;
    trans.dist_ = 1.0;

    try {
      // try to split and get rmsd and score
      if (split_results.size() == 9) { // input of num, score, rmsd, trans
        trans.score_.totalScore_ = std::stof(split_results[1].c_str());
        std::cout << "loaded score " << trans.score_.totalScore_ << std::endl;
      } else {
        std::vector<std::string> split_results2;
        boost::split(split_results2, line, boost::is_any_of(":|\t"), boost::token_compress_on);
        if(split_results2.size() > 3) {
          trans.score_.totalScore_ = std::stof(split_results2[1].c_str());
        }
      }
    } catch(std::invalid_argument& e) {
      std::cerr << "Couldn't get score from trans file, continue without\n";
    }

    try {
      // try to split and get rmsd and score
      if (split_results.size() == 9) { // input of num, score, rmsd, trans
        trans.dist_ = std::stof(split_results[2].c_str());
      } else {
        std::vector<std::string> split_results2;
        boost::split(split_results2, line, boost::is_any_of(":|\t"), boost::token_compress_on);
        if(split_results2.size() > 3) {
          trans.dist_ = std::stof(split_results2[2].c_str());
        }
      }
    } catch(std::invalid_argument& e) {
      std::cerr << "Couldn't get RMSD from trans file, continue without\n";
    }
  }
  return !transFile.eof();
}

void BBContainer::computeRanges(int margin) {
  for (size_t i = 0; i < numOfBBs_; i++) {
    bbs_[i]->computeRanges(bbs_, margin);
  }
}
