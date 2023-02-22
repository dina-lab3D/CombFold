/**
 * \file MainCombDock.cc
 *
 * Author: Dina Schneidman, Yuval Inbar
 *
 */

#include "ClusterSBBs.h"
#include "HierarchicalFold.h"
#include "BBContainer.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <vector>
#include <stdio.h>
#include <string.h>
//#include <memory>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

//#include <gperftools/profiler.h>

// forward declarations
void clusterResults(const std::string resFileName,
                    const std::vector<std::shared_ptr<const BB>>& bbs,
                    double clusterRMSD);


int main(int argc, char *argv[]) {
  // output arguments
  for (int i = 0; i < argc; i++) std::cout << argv[i] << " ";
  std::cout << std::endl;

  // input parsing
  double penetrationThr = -5.0;
  float maxBackboneCollisionPerChain = 0.1;
  float minTemperatureToConsiderCollision = 0;
  unsigned int maxResultPerResSet = 0;


  std::string outFileNamePrefix;
  double constraintsRatio = 0.7;
  double clusterRMSD = 5.0;

  // positional
  std::string suFileName;
  std::string prefix;
  int transNumToRead;
  int bestK;
  std::string constraintsFileName;

  po::options_description desc("Usage: <subunitsFileList> <transFilesPrefix> \
<transNumToRead> <bestKeachStep> <constraintsFile>");

  // optional
  desc.add_options()
    ("help,h","CombDock help")
    ("version", "CombDock 2.0 2019")
    ("penetrationThr,p", po::value<double>(&penetrationThr)->default_value(-5.0),
     "maximum allowed penetration between subunit surfaces (default = -5.0)")
    ("constraintsRatio,r", po::value<double>(&constraintsRatio)->default_value(0.7),
     "constraints ratio (default = 0.7)")
    ("clusterRMSD,c", po::value<double>(&clusterRMSD)->default_value(5.0),
     "final clustering RMSD (default = 5.0)")

      ("maxBackboneCollisionPerChain,b", po::value<float>(&maxBackboneCollisionPerChain)->default_value(0.1),
       "Max percentage(0 to 1) of backbone atoms of a chain that can collide with another chain(default=0.1)")
      ("minTemperatureToConsiderCollision,t", po::value<float>(&minTemperatureToConsiderCollision)->default_value(0),
       "Minimal Bfactor required for atom to be considered when calculating collisions(default=0)")
      ("maxResultPerResSet,j", po::value<unsigned int>(&maxResultPerResSet)->default_value(0),
       "number of results saved for each calculated combination of subunits (default=k)")

    ("outputFileNamePrefix,o", po::value<std::string>(&outFileNamePrefix)->default_value("combdock"),
     "output file name, default name combdock.res");

  // required options: currently 5
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("SUlist", po::value<std::string>(&suFileName)->required(),"SU list file name")
    ("transFilesPrefix", po::value<std::string>(&prefix)->required(),"Trans files prefix")
    ("transNumToRead", po::value<int>(&transNumToRead)->required(), "# of tranformations")
    ("bestK", po::value<int>(&bestK)->required(), "bestK")
    ("constraintsFile", po::value<std::string>(&constraintsFileName)->required(),
     "constraints file name");

  po::options_description cmdline_options;
  cmdline_options.add(desc).add(hidden);

  po::positional_options_description p;
  p.add("SUlist", 1);
  p.add("transFilesPrefix", 1);
  p.add("transNumToRead", 1);
  p.add("bestK", 1);
  p.add("constraintsFile", 1);

  po::variables_map vm;

  try {
    po::store(po::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);

    if(vm.count("help")) {
      std::cout << desc << "\n";
      return 0;
    }
    po::notify(vm);
  }
  catch(po::required_option& e) {
    std::cout << desc << "\n";
      return 0;
  }
  catch(po::error& e) {
    std::cout << desc << "\n";
    return 0;
  }

  if(maxResultPerResSet == 0)
      maxResultPerResSet = bestK;

  // done parsing

  auto start = std::chrono::high_resolution_clock::now();
  HierarchicalFold::timerAll_.reset();

  std::vector<std::string> pdbs;
  std::vector<TransformationAndScore *> trans1, trans2;

  std::cerr << "Before process input" << std::endl;
  std::string argv_str(argv[0]);
  std::string base = argv_str.substr(0, argv_str.find_last_of("/"));
  std::string chemLibFileName = base + "/chem.lib";
  BBContainer bbContainer(suFileName, penetrationThr, chemLibFileName);
  bbContainer.readTransformationFiles(prefix, transNumToRead);

  std::cout << "Starting HierarchicalFold" << std::endl;
  HierarchicalFold hierarchalFold(bbContainer.getBBs(), bestK, maxResultPerResSet, minTemperatureToConsiderCollision,
                                  maxBackboneCollisionPerChain);
  for (unsigned int i = 0; i < bbContainer.getBBsNumber(); i++) {
    std::shared_ptr<SuperBB> sbb =
      std::make_shared<SuperBB>(bbContainer.getBB(i), bbContainer.getBBsNumber());
    hierarchalFold.add(sbb, i);
  }

  // read constraints
  hierarchalFold.readConstraints(constraintsFileName, constraintsRatio);

  HierarchicalFold::timer_.reset();

  hierarchalFold.outputConnectivityGraph("graph.txt");
  hierarchalFold.checkConnectivity();
//    ProfilerStart("nameOfProfile.log");
  auto startBeforeFold = std::chrono::high_resolution_clock::now();
  hierarchalFold.fold(outFileNamePrefix);
//    ProfilerStop();
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end-start;
  std::chrono::duration<double> diffFold = end-startBeforeFold;
  std::cout << "Overall time "<< diff.count() << " s\n";
  std::cout << "only fold time "<< diffFold.count() << " s\n";
  std::cerr << "countResults_ " << HierarchicalFold::countResults_ << std::endl;
  std::cerr << "timer_ " << HierarchicalFold::timer_ << std::endl;
  std::cerr << "timerAll_ " << HierarchicalFold::timerAll_ << std::endl;
  //SuperBB::reportCounters();
  std::string outFileName = outFileNamePrefix + ".res";
  clusterResults(outFileName, bbContainer.getBBs(), clusterRMSD);

  return 0;
}

void clusterResults(const std::string resFileName,
                    const std::vector<std::shared_ptr<const BB>>& bbs,
                    double clusterRMSD) {

  std::ifstream resFile(resFileName);
  if (!resFile) {
    std::cerr << "Can't open results file " << resFileName << std::endl;
    return;
  }
  if(resFile.peek() == std::ifstream::traits_type::eof()) {
    std::cerr << "Results file empty " << resFileName << std::endl;
    return;
  }

  std::cout << "Start final clustering " << resFileName << " rmsd " << clusterRMSD << std::endl;

  ClusterSBBs cluster(clusterRMSD);
    char line[10000];
    std::vector<SuperBB *> sbbs;
    while (!resFile.eof()) {
        if (!resFile.getline(line, 5000)) break;
        std::string str(index(line, '[') + 1);
        //std::cerr << str << std::endl;
        std::stringstream s(str, std::stringstream::in);
        SuperBB *sbb = new SuperBB(bbs[0], bbs.size());
        sbbs.push_back(sbb);
        //std::unique_ptr<SuperBB> sbb(new SuperBB(BB::bbs_[0]));
        RigidTrans3 tr;
        char *tScore = strstr(line, "transScore_");
        if (tScore != NULL) {
            sscanf(tScore, "transScore_ %d", &sbb->transScore_);
        }
        char *rmsd = strstr(line, "rmsd");
        if (rmsd != NULL) {
            sscanf(rmsd, "rmsd %f", &sbb->rmsd_);
        }
        char *multPen = strstr(line, "multPen_");
        if (multPen != NULL) {
            sscanf(multPen, "multPen_ %d", &sbb->multPen_);
        }
        char *singlePen = strstr(line, "singlePen_");
        if (singlePen != NULL) {
            sscanf(singlePen, "singlePen_ %d", &sbb->singlePen_);
        }
        char *newScore = strstr(line, "newScore_");
        if (newScore != NULL) {
            sscanf(newScore, "newScore_ %f", &sbb->newScore_);
        }
        char *numOfBuried = strstr(line, "numOfBuried_");
        if (numOfBuried != NULL) {
            sscanf(numOfBuried, "numOfBuried_ %d", &sbb->numOfBuried_);
        }

        char *numOfHydBuried = strstr(line, "numOfHydBuried_");
        if (numOfHydBuried != NULL) {
            sscanf(numOfHydBuried, "numOfHydBuried_ %d", &sbb->numOfHydBuried_);
        }

        char *maxPen = strstr(line, "maxPen_");
        if (maxPen != NULL) {
            sscanf(maxPen, "maxPen_ %f", &sbb->maxPen_);
        }

        char *restraintsRatio =  strstr(line, "restraintsRatio_");
        if (restraintsRatio != NULL) {
          float r;
          sscanf(restraintsRatio, "restraintsRatio_ %f", &r);
          sbb->setRestraintsRatio(r);
        }

        sbb->yetAnother_ = sbb->newScore_ * ((float) sbb->singlePen_ / (sbb->singlePen_ - sbb->multPen_));
        for (unsigned int i = 1; i < bbs.size(); i++) {
            sbb->bbs_.push_back(bbs[i]);
            sbb->trans_.push_back(tr);
        }
        sbb->setSize(bbs.size());
        for (unsigned int i = 0; i < bbs.size(); i++) {
            //std::cerr << " i " << i << std::endl;
            char c;
            int index;
            RigidTrans3 trans;
            s >> index >> c >> trans >> c >> c;
            //std::cerr << "index " << index << " trans " << trans << " c " << c << std::endl;
            sbb->trans_[index] = trans;
        }
        ///
        cluster.addSBB(sbb);
    }
    cluster.merge();
    std::ofstream outFile("clusters.res");
    cluster.reportAll(outFile);
    outFile.close();
    for (unsigned int i = 0; i < sbbs.size(); i++) delete sbbs[i];
    sbbs.clear();
}
