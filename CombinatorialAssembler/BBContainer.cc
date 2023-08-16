#include "BBContainer.h"
#include <boost/algorithm/string.hpp>

namespace {
std::string trim_extension(const std::string file_name) {
    if (file_name[file_name.size() - 4] == '.')
        return file_name.substr(0, file_name.size() - 4);
    return file_name;
}
} // namespace

BBContainer::BBContainer(const std::string SUFileName, std::string chemLibFileName, float minTempFactor) {
    readSUFile(SUFileName);

    // prepare ChemLib
    ChemLib chemLib(chemLibFileName);

    // read the building blocks
    bbs_.reserve(numOfBBs_);
    for (unsigned int i = 0; i < numOfBBs_; i++) {
        bbs_.push_back(std::make_shared<BB>(i, pdbs_[i], groupIDs_[i], chemLib, 0.5, 5.0, minTempFactor));
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
            if (i == j)
                continue;
            std::string inFileName = transFilePrefix + trim_extension(pdbs_[i]) + "_plus_" + trim_extension(pdbs_[j]);
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
            boost::split(split_results, line, boost::is_any_of(" "), boost::token_compress_on);
            std::string pdbName(split_results[0]);
            std::cout << "PDBname " << pdbName << ":" << line << std::endl;
            numOfBBs_++;
            pdbs_.push_back(pdbName);
            // group assignment
            int groupID = 0;
            if (split_results.size() == 2) {
                groupID = stoi(split_results[1]);
                if (groupID <= 0) {
                    std::cerr << "Group ID must be positive" << std::endl;
                }
            }
            groupIDs_.push_back(groupID);
        }
    }
    return numOfBBs_;
}

bool BBContainer::readTrans(TransformationAndScore &trans, std::ifstream &transFile) {
    std::string line;
    if (!transFile.eof()) {
        getline(transFile, line);
        boost::trim(line); // remove all spaces
        // skip comments
        if (line[0] == '#' || line[0] == '\0')
            return false;

        std::vector<std::string> split_results;
        boost::split(split_results, line, boost::is_any_of(":|\t"), boost::token_compress_on);
        if (split_results.size() != 4) {
            std::cerr << "Wrong number of fields in transformation file, should be: "
                      << "index(int) | score(float) | comment | transformation(space seperated 6 floats)" << std::endl;
            return false;
        }

        boost::trim(split_results[3]);
        std::vector<std::string> splitted_transformation;
        boost::split(splitted_transformation, split_results[3], boost::is_any_of(" "), boost::token_compress_on);

        if (splitted_transformation.size() != 6)
            return false;

        // extract trans
        RigidTrans3 tr(Vector3(std::stof(splitted_transformation[0].c_str()),
                               std::stof(splitted_transformation[1].c_str()),
                               std::stof(splitted_transformation[2].c_str())),
                       Vector3(std::stof(splitted_transformation[3].c_str()),
                               std::stof(splitted_transformation[4].c_str()),
                               std::stof(splitted_transformation[5].c_str())));
        trans.refFrame_ = tr;
        trans.score_.totalScore_ = 1.0;
        trans.dist_ = 1.0;

        trans.score_.totalScore_ = std::stof(split_results[1].c_str());
        std::cout << "loaded score: " << trans.score_.totalScore_ << std::endl;
    }
    return !transFile.eof();
}
