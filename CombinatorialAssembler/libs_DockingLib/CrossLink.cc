#include "CrossLink.h"

#include <boost/algorithm/string.hpp>

#include <iostream>

int readCrossLinkFile(const std::string& fileName,
                      std::vector<CrossLink>& crossLinks,
                      bool addReverse) {
  std::ifstream s(fileName);
  if (!s) {
    std::cerr << "Can't find cross links file " << fileName << std::endl;
    exit(0);
  }
  CrossLink cl;
  while (s >> cl) {
    if(! cl.isInList(crossLinks)) { // check for duplicates

      if(!addReverse) {
        // check reversed for duplicates
        CrossLink clReversed = cl.getReversed();
        if(! clReversed.isInList(crossLinks)) {
          crossLinks.push_back(cl);
        }
      } else {
        crossLinks.push_back(cl);
      }

    } else {
      // TODO: check if shorter distance restraint and replace if needed
      std::cerr << "Duplicate cross link " << cl << std::endl;
    }

    // add reverse
    if(addReverse) {
      CrossLink clReversed = cl.getReversed();
      if(! clReversed.isInList(crossLinks)) {
        crossLinks.push_back(clReversed);
      }
    }
  }
  return crossLinks.size();
}

void writeCrossLinkFile(const std::string& fileName,
                        const std::vector<CrossLink>& crossLinks) {
  std::ofstream ofile(fileName);
  for (unsigned int i = 0; i < crossLinks.size(); i++) {
    ofile << crossLinks[i];
    ofile << std::endl;
  }
  ofile.close();
}

void writeXLAnalizerFile(const std::string& fileName,
                         const std::vector<CrossLink>& crossLinks) {
  std::ofstream ofile(fileName);
  ofile << "id,Protein1,Protein2,AbsPos1,AbsPos2,score" << std::endl;
  for (unsigned int i = 0; i < crossLinks.size(); i++) {
    ofile << i+1 << "," << crossLinks[i].getChain1() << "," << crossLinks[i].getChain2() <<  ","
          << crossLinks[i].getResidue1() << "," << crossLinks[i].getResidue2()
          << ",100" << std::endl;
    ofile << std::endl;
  }
  ofile.close();
}

std::ostream& operator<<(std::ostream& s, const CrossLink& cl) {
  s << cl.residueNumber1_ << " ";
  if (cl.chainId1_ == " ")
    s << "-";
  else
    s << cl.chainId1_;

  s << " " << cl.residueNumber2_ << " ";
  if (cl.chainId2_ == " ")
    s << "-";
  else
    s << cl.chainId2_;
  s << " " << cl.minDistance_ << " " << cl.maxDistance_;
  return s;
}

std::istream& operator>>(std::istream& s, CrossLink& cl) {

  std::string line;
  std::getline(s, line);
  boost::trim(line); // remove spaces at the beginning/end of the line
  if(line.length()==0) return s;
  // skip comments
  if (line[0] == '#' || line[0] == '\0' || !isdigit(line[0])) return s;

  std::vector<std::string> split_results;
  boost::split(split_results, line, boost::is_any_of("\t "),
               boost::token_compress_on);

  if (split_results.size() >= 5 && split_results.size() <= 7) {
    cl.residueSequenceID1_ = split_results[0];
    cl.chainId1_ = split_results[1];
    cl.residueSequenceID2_ = split_results[2];
    cl.chainId2_ = split_results[3];
    cl.residueNumber1_ = std::stoi(cl.residueSequenceID1_);
    cl.residueNumber2_ = std::stoi(cl.residueSequenceID2_);
    cl.maxDistance_ = std::stof(split_results[4]);

    if(split_results.size() >= 6) {
      cl.minDistance_ = std::stof(split_results[4]);
      cl.maxDistance_ = std::stof(split_results[5]);
    }

    if(split_results.size() == 7) {
      cl.weight_ = std::stof(split_results[6]);
    }
  }

  if (cl.chainId1_ == "-") cl.chainId1_ = " ";
  if (cl.chainId2_ == "-") cl.chainId2_ = " ";
  //  std::cerr << cl << std::endl;
  return s;
}

bool CrossLink::operator==(const CrossLink& cl) const {
  return (residueNumber1_ == cl.residueNumber1_ && residueNumber2_ == cl.residueNumber2_ &&
          residueSequenceID1_ == cl.residueSequenceID1_ && residueSequenceID2_ == cl.residueSequenceID2_ &&
          chainId1_ == cl.chainId1_ && chainId2_ == cl.chainId2_);
}

bool CrossLink::isInList(const std::vector<CrossLink>& crossLinks) const {
  for(unsigned int i = 0; i<crossLinks.size(); i++)
    if(*this == crossLinks[i])
      return true;
  return false;
}

CrossLink CrossLink::getReversed() const {
  CrossLink reversed(residueNumber2_, chainId2_,
                     residueNumber1_, chainId1_,
                     minDistance_, maxDistance_, weight_);
  return reversed;
}
