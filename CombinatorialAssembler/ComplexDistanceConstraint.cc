#include "ComplexDistanceConstraint.h"

#include "BB.h"

void ComplexDistanceConstraint::addConstraint(int suInd1, int suInd2, Vector3 receptorAtom, Vector3 ligandAtom,
                                              float maxDistance, float minDistance) {
    int ind = suInd1 * noOfSUs_ + suInd2;
    constraints_[ind].push_back(DistanceRestraint(receptorAtom, ligandAtom, maxDistance, minDistance));
}

int ComplexDistanceConstraint::addChainConnectivityConstraints() {
    for (int suInd1 = 0; suInd1 < noOfSUs_; suInd1++) {
        for (int suInd2 = suInd1 + 1; suInd2 < noOfSUs_; suInd2++) {
            std::vector<std::pair<char, std::pair<int, int>>> constraints;
            int counter = 0;
            float epsilon = 20.0; // for shorter connections
            bbs_[suInd1]->getChainConnectivityConstraints(*bbs_[suInd2], constraints);
            for (int i = 0; i < (int)constraints.size(); i++) {
                int su1resIndex = constraints[i].second.first;
                int su2resIndex = constraints[i].second.second;

                int receptorAtomIndex = bbs_[suInd1]->allAtoms_.getFirstAtomEntryForResIndex(
                    constraints[i].first, std::to_string(su1resIndex));
                int ligandAtomIndex = bbs_[suInd2]->allAtoms_.getFirstAtomEntryForResIndex(constraints[i].first,
                                                                                           std::to_string(su2resIndex));

                Vector3 receptorAtom = bbs_[suInd1]->getChemAtomByIndex(receptorAtomIndex).position();
                Vector3 ligandAtom = bbs_[suInd2]->getChemAtomByIndex(ligandAtomIndex).position();

                int sequenceDist = std::fabs(su1resIndex - su2resIndex);
                if (sequenceDist < 100) {
                    float maxDist = sequenceDist * 2.0;
                    if (sequenceDist <= 5)
                        maxDist += epsilon;
                    addConstraint(suInd1, suInd2, receptorAtom, ligandAtom, maxDist);
                    addConstraint(suInd2, suInd1, ligandAtom, receptorAtom, maxDist);
                    std::cerr << "DISTANCE CONSTRAINT ADDED: " << constraints[i].first << " su1 endpoint "
                              << su1resIndex << " su2 endpoint " << su2resIndex << " " << maxDist << std::endl;
                    counter++;
                }
            }
            numberOfConstraints_ += counter;
        }
    }
    return numberOfConstraints_;
}

std::set<int> ComplexDistanceConstraint::getSUs(const std::string residueSequenceID, const std::string chains,
                                                int maxOffset) const {

    int rOffset = 0;
    std::set<int> ret;
    // iterate chains
    for (unsigned int chainIndex = 0; chainIndex < chains.size(); chainIndex++) {
        // for each chain find corresponding SU
        for (int suIndex = 0; suIndex < noOfSUs_; suIndex++) {
            int atomIndex = bbs_[suIndex]->allAtoms_.getClosestAtomEntryForResIndex(
                chains[chainIndex], residueSequenceID, maxOffset, rOffset);
            if (atomIndex != -1) { // found
                ret.insert(suIndex);
                continue;
            }
        }
    }
    return ret;
}

int ComplexDistanceConstraint::readRestraintsFile(const std::string fileName) {
    // read cross links
    std::vector<CrossLink> crosslinks;
    int xnum = readCrossLinkFile(fileName, crosslinks);
    Logger::infoMessage() << "# of xlinks " << xnum << " read from file" << fileName << std::endl;
    int MAX_OFFSET = 20;

    // map them to BB pairs
    for (unsigned int i = 0; i < crosslinks.size(); i++) {

        // SUs that have xlinks endpoints (su1, su2)
        std::string residueSequenceID1 = std::to_string(crosslinks[i].getResidue1());
        std::set<int> su1 = getSUs(residueSequenceID1, crosslinks[i].getChain1());
        if (su1.size() == 0)
            su1 = getSUs(residueSequenceID1, crosslinks[i].getChain1(), MAX_OFFSET);
        std::string residueSequenceID2 = std::to_string(crosslinks[i].getResidue2());
        std::set<int> su2 = getSUs(residueSequenceID2, crosslinks[i].getChain2());
        if (su2.size() == 0)
            su2 = getSUs(residueSequenceID2, crosslinks[i].getChain2(), MAX_OFFSET);
        int maxOffset = 0;

        for (auto siter1 = su1.begin(); siter1 != su1.end(); siter1++) {
            std::vector<Vector3> p1;
            int offset = 0;
            // find all chains for su1
            for (unsigned int chainIndex1 = 0; chainIndex1 < crosslinks[i].getChain1().size(); chainIndex1++) {
                char chainId1 = crosslinks[i].getChain1()[chainIndex1];
                int rOffset = 0;
                // int maxOffset = 0;

                int receptorAtomIndex = bbs_[*siter1]->allAtoms_.getClosestAtomEntryForResIndex(
                    chainId1, residueSequenceID1, MAX_OFFSET, rOffset);
                if (receptorAtomIndex == -1)
                    continue;
                offset += std::abs(rOffset);
                Vector3 rcoord = bbs_[*siter1]->getChemAtomByIndex(receptorAtomIndex).position();
                p1.push_back(rcoord);
            }

            for (auto siter2 = su2.begin(); siter2 != su2.end(); siter2++) {
                if (*siter1 == *siter2)
                    continue; // same SU xlink
                std::vector<Vector3> p2;
                for (unsigned int chainIndex2 = 0; chainIndex2 < crosslinks[i].getChain2().size(); chainIndex2++) {
                    char chainId2 = crosslinks[i].getChain2()[chainIndex2];
                    int lOffset = 0;
                    int ligandAtomIndex = bbs_[*siter2]->allAtoms_.getClosestAtomEntryForResIndex(
                        chainId2, residueSequenceID2, MAX_OFFSET, lOffset);
                    if (ligandAtomIndex == -1)
                        continue;
                    offset += std::abs(lOffset);
                    Vector3 lcoord = bbs_[*siter2]->getChemAtomByIndex(ligandAtomIndex).position();
                    p2.push_back(lcoord);
                    if (offset > maxOffset)
                        maxOffset = offset;
                }

                if (p1.size() > 0 && p2.size() > 0) {
                    // use 3.0 instead of 3.8A for CA-CA distance as it is not likely to be a straight line
                    DistanceRestraint d1(p1, p2,
                                         crosslinks[i].getMaxDistance(), // + maxOffset*3.0,
                                         crosslinks[i].getMinDistance(), crosslinks[i].getWeight());
                    int ind1 = *siter1 * noOfSUs_ + *siter2;
                    restraints_[ind1].push_back(d1);

                    DistanceRestraint d2(p2, p1, crosslinks[i].getMaxDistance() + maxOffset * 3.0,
                                         crosslinks[i].getMinDistance(), crosslinks[i].getWeight());
                    int ind2 = *siter2 * noOfSUs_ + *siter1;
                    restraints_[ind2].push_back(d2);
                    // counter++;
                    std::cerr << " adding restraint to " << *siter1 << "x" << *siter2 << " :" << p1.size() << " W "
                              << p2.size() << " R1 " << residueSequenceID1 << crosslinks[i].getChain1()[0] << " R2 "
                              << residueSequenceID2 << crosslinks[i].getChain2()[0] << " dist "
                              << d1.getMaxDistance() + maxOffset * 3.0 << " " << d1.getMinDistance() << std::endl;
                    crossLinks_.push_back(crosslinks[i]);
                }
            }
        }
    }
    numberOfRestraints_ = crossLinks_.size();
    std::string readFileName = fileName + ".read";
    writeCrossLinkFile(readFileName, crossLinks_);
    return crossLinks_.size();
}
/*
int ComplexDistanceConstraint::readRestraintsFile(const std::string fileName) {
  std::ifstream file(fileName);
  if(!file) {
    std::cerr << "No constraints file " << fileName << std::endl;
    exit(0);
  }
  std::ofstream xlAnalizerFile("xlAnalizer.csv");
  xlAnalizerFile << "id,Protein1,Protein2,AbsPos1,AbsPos2,score" << std::endl;
  int id = 1;
  for (int suInd1 = 0; suInd1 < noOfSUs_; suInd1++) {
    for (int suInd2 = suInd1+1; suInd2 < noOfSUs_; suInd2++) {
      int receptorAtomIndex, ligandAtomIndex;
      float minDist(0), maxDist(0);
      Vector3 receptorAtom,ligandAtom;
      int counter = 0;
      int maxOffset = 20;

      std::string line;
      while(!file.eof()) {
        getline(file, line);
        boost::trim(line); // remove all spaces
        if(line.length()==0) continue;
        // skip comments
        if (line[0] == '#' || line[0] == '\0' || !isdigit(line[0])) continue;

        std::vector<std::string> split_results;
        boost::split(split_results, line, boost::is_any_of("\t "),
                     boost::token_compress_on);

        // Option 1: read atom indices
        if (split_results.size() == 3 || split_results.size() == 4) {
          receptorAtomIndex = atoi(split_results[0].c_str());
          ligandAtomIndex = atoi(split_results[1].c_str());

          receptorAtom = bbs_[suInd1]->getChemAtom(receptorAtomIndex).position();
          ligandAtom = bbs_[suInd2]->getChemAtom(ligandAtomIndex).position();

          float dist = atof(split_results[2].c_str());
          if (split_results.size() == 4) {
            maxDist = atof(split_results[3].c_str());
            minDist = dist;
          } else {
            maxDist = dist;
          }
        }

        // Option 2: read residue index and chain ID
        if (split_results.size() == 5 || split_results.size() == 6) {
          std::string residueSequenceID1 = split_results[0].c_str();
          char chainId1 = split_results[1][0];
          std::string residueSequenceID2 = split_results[2].c_str();
          char chainId2 = split_results[3][0];

          int offset = 0;
          int rOffset = 0;
          receptorAtomIndex =
            bbs_[suInd1]->allAtoms_.getClosestAtomEntryForResIndex(chainId1,
                                                                   residueSequenceID1,
                                                                   maxOffset, rOffset);
          if(receptorAtomIndex == -1) continue;
          offset += std::abs(rOffset);

          int lOffset = 0;
          ligandAtomIndex =
            bbs_[suInd2]->allAtoms_.getClosestAtomEntryForResIndex(chainId2,
                                                                   residueSequenceID2,
                                                                   maxOffset, lOffset);
          if(ligandAtomIndex == -1) continue;
          offset += std::abs(lOffset);

          receptorAtom = bbs_[suInd1]->getChemAtomByIndex(receptorAtomIndex).position();
          ligandAtom = bbs_[suInd2]->getChemAtomByIndex(ligandAtomIndex).position();

          float dist = atof(split_results[4].c_str());
          if (split_results.size() == 6) {
            maxDist = atof(split_results[5].c_str());
            minDist = dist;
          } else {
            maxDist = dist;
          }

          if(offset > 0) {
            // use 3.0 instead of 3.8A for CA-CA distance as it is not likely to be a straight line
            maxDist += offset * 3.0;
          }

          if(offset <= 2) {
            xlAnalizerFile << id << "," << chainId1 << "," << chainId2 <<  ","
                           << residueSequenceID1 << "," << residueSequenceID2 << ",100" << std::endl;
            id++;
          }

          std::cerr << "DISTANCE RESTRAINT ADDED: " << chainId1 << residueSequenceID1
                <<  bbs_[suInd1]->getChemAtomByIndex(receptorAtomIndex).residueName()
                << " " << chainId2 << residueSequenceID2
                <<  bbs_[suInd2]->getChemAtomByIndex(ligandAtomIndex).residueName()
                << " " << maxDist << std::endl;
        }

        addRestraint(suInd1, suInd2, receptorAtom, ligandAtom, maxDist, minDist);
        addRestraint(suInd2, suInd1, ligandAtom, receptorAtom, maxDist, minDist);
        counter++;
        std::cerr << " adding restraint to " << suInd1 << "x" << suInd2 << " :" << receptorAtom << " W " << ligandAtom
<< " dist " << maxDist << " " << minDist << std::endl;
      }
      file.clear();
      file.seekg(0, std::ios::beg);
      numberOfRestraints_ += counter;
      std::cout << "INFO " << counter << " distances were found for "
                << suInd1 << " " << suInd2 << std::endl;
    }
  }
  file.close();
  xlAnalizerFile.close();
  std::cout << "Total number of restraints: " << numberOfRestraints_ << std::endl;
  return numberOfRestraints_;
}

int ComplexDistanceConstraint::readCrossLinksFile(const std::string fileName) {
  std::cout << "start readling xlinks file " << fileName << std::endl;
  std::vector<CrossLink> crossLinks;
  int xnum = readCrossLinkFile(fileName, crossLinks);
  std::cout << "# of xlinks " << xnum << std::endl;
  writeXLAnalyzerFile("xlAnalizer.csv");

  // find SUs for each cross link
  for (unsigned int i = 0; i < crossLinks.size(); i++) {

  int id = 1;
  for (int suInd1 = 0; suInd1 < noOfSUs_; suInd1++) {
    for (int suInd2 = suInd1+1; suInd2 < noOfSUs_; suInd2++) {
      int receptorAtomIndex, ligandAtomIndex;
      float minDist(0), maxDist(0);
      Vector3 receptorAtom,ligandAtom;
      int counter = 0;
      int maxOffset = 20;



          int offset = 0;
          int rOffset = 0;
          receptorAtomIndex =
            bbs_[suInd1]->allAtoms_.getClosestAtomEntryForResIndex(chainId1,
                                                                   residueSequenceID1,
                                                                   maxOffset, rOffset);
          if(receptorAtomIndex == -1) continue;
          offset += std::abs(rOffset);

          int lOffset = 0;
          ligandAtomIndex =
            bbs_[suInd2]->allAtoms_.getClosestAtomEntryForResIndex(chainId2,
                                                                   residueSequenceID2,
                                                                   maxOffset, lOffset);
          if(ligandAtomIndex == -1) continue;
          offset += std::abs(lOffset);

          receptorAtom = bbs_[suInd1]->getChemAtomByIndex(receptorAtomIndex).position();
          ligandAtom = bbs_[suInd2]->getChemAtomByIndex(ligandAtomIndex).position();

          float dist = atof(split_results[4].c_str());
          if (split_results.size() == 6) {
            maxDist = atof(split_results[5].c_str());
            minDist = dist;
          } else {
            maxDist = dist;
          }

          if(offset > 0) {
            // use 3.0 instead of 3.8A for CA-CA distance as it is not likely to be a straight line
            maxDist += offset * 3.0;
          }



          std::cerr << "DISTANCE RESTRAINT ADDED: " << chainId1 << residueSequenceID1
                <<  bbs_[suInd1]->getChemAtomByIndex(receptorAtomIndex).residueName()
                << " " << chainId2 << residueSequenceID2
                <<  bbs_[suInd2]->getChemAtomByIndex(ligandAtomIndex).residueName()
                << " " << maxDist << std::endl;
        }

        addRestraint(suInd1, suInd2, receptorAtom, ligandAtom, maxDist, minDist);
        addRestraint(suInd2, suInd1, ligandAtom, receptorAtom, maxDist, minDist);
        counter++;
        std::cerr << " adding restraint to " << suInd1 << "x" << suInd2 << " :" << receptorAtom << " W " << ligandAtom
<< " dist " << maxDist << " " << minDist << std::endl;
      }
      file.clear();
      file.seekg(0, std::ios::beg);
      numberOfRestraints_ += counter;
      std::cout << "INFO " << counter << " distances were found for "
                << suInd1 << " " << suInd2 << std::endl;
    }
  }
  file.close();

  std::cout << "Total number of restraints: " << numberOfRestraints_ << std::endl;
  return numberOfRestraints_;
}
*/
