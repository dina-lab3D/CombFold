#include "BBContainer.h"
#include "HierarchicalFold.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <vector>
// #include <memory>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

// #include <gperftools/profiler.h>

int main(int argc, char *argv[]) {
    // output arguments
    for (int i = 0; i < argc; i++)
        std::cout << argv[i] << " ";
    std::cout << std::endl;

    // input parsing
    double penetrationThr;
    float maxBackboneCollisionPerChain;
    float minTemperatureToConsiderCollision;
    unsigned int maxResultPerResSet;

    std::string outFileNamePrefix;
    double restraintsRatio;
    double clusterRMSD;

    // positional
    std::string suFileName;
    std::string transFilesPrefix;
    int transNumToRead;
    int bestK;
    std::string constraintsFileName;

    po::options_description desc("Usage: <subunitsFileList> <transFilesPrefix> \
<transNumToRead> <bestKeachStep> <constraintsFile>");

    // optional
    desc.add_options()("help,h", "Combinatoral Assembly help")("version", "CombFold 1.0 2022")(
        "penetrationThr,p", po::value<double>(&penetrationThr)->default_value(-1.0),
        "maximum allowed penetration between subunit surfaces (default = -1.0)")(
        "restraintsRatio,r", po::value<double>(&restraintsRatio)->default_value(0.1),
        "constraints ratio (default = 0.1)")("clusterRMSD,c", po::value<double>(&clusterRMSD)->default_value(5.0),
                                             "final clustering RMSD (default = 5.0)")

        ("maxBackboneCollisionPerChain,b", po::value<float>(&maxBackboneCollisionPerChain)->default_value(0.1),
         "Max percentage(0 to 1) of backbone atoms of a chain that can collide with another chain(default=0.1)")(
            "minTemperatureToConsiderCollision,t",
            po::value<float>(&minTemperatureToConsiderCollision)->default_value(0),
            "Minimal Bfactor required for atom to be considered when calculating collisions(default=0)")(
            "maxResultPerResSet,j", po::value<unsigned int>(&maxResultPerResSet)->default_value(0),
            "number of results saved for each calculated combination of subunits (default=k)")

            ("outputFileNamePrefix,o", po::value<std::string>(&outFileNamePrefix)->default_value("output"),
             "output file name, default name output.res");

    // required options: currently 5
    po::options_description hidden("Hidden options");
    hidden.add_options()("SUlist", po::value<std::string>(&suFileName)->required(), "SU list file name")(
        "transFilesPrefix", po::value<std::string>(&transFilesPrefix)->required(),
        "Trans files prefix")("transNumToRead", po::value<int>(&transNumToRead)->required(),
                              "# of tranformations")("bestK", po::value<int>(&bestK)->required(), "bestK")(
        "constraintsFile", po::value<std::string>(&constraintsFileName)->required(), "constraints file name");

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
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }
        po::notify(vm);
    } catch (po::required_option &e) {
        std::cout << desc << "\n";
        return 0;
    } catch (po::error &e) {
        std::cout << desc << "\n";
        return 0;
    }

    if (maxResultPerResSet == 0)
        maxResultPerResSet = bestK;

    // done parsing

    auto start = std::chrono::high_resolution_clock::now();
    HierarchicalFold::timerAll_.reset();

    std::cerr << "Before process input" << std::endl;
    std::string argv_str(argv[0]);
    std::string base = argv_str.substr(0, argv_str.find_last_of("/"));
    std::string chemLibFileName = base + "/chem_params.txt";
    BBContainer bbContainer(suFileName, chemLibFileName, minTemperatureToConsiderCollision);
    bbContainer.readTransformationFiles(transFilesPrefix, transNumToRead);

    std::cout << "Starting HierarchicalFold" << std::endl;
    HierarchicalFold hierarchalFold(bbContainer, bestK, maxResultPerResSet, minTemperatureToConsiderCollision,
                                    maxBackboneCollisionPerChain, penetrationThr, restraintsRatio);

    // read constraints
    hierarchalFold.readConstraints(constraintsFileName);

    HierarchicalFold::timer_.reset();

    hierarchalFold.outputConnectivityGraph("graph.txt");
    hierarchalFold.checkConnectivity();
    //    ProfilerStart("nameOfProfile.log");
    auto startBeforeFold = std::chrono::high_resolution_clock::now();
    hierarchalFold.fold(outFileNamePrefix);
    //    ProfilerStop();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::chrono::duration<double> diffFold = end - startBeforeFold;
    std::cout << "Overall time " << diff.count() << " s\n";
    std::cout << "only fold time " << diffFold.count() << " s\n";
    std::cerr << "countResults_ " << HierarchicalFold::countResults_ << std::endl;
    std::cerr << "timer_ " << HierarchicalFold::timer_ << std::endl;
    std::cerr << "timerAll_ " << HierarchicalFold::timerAll_ << std::endl;
    // SuperBB::reportCounters();

    return 0;
}
