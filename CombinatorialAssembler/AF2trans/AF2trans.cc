#include "Atom.h"
#include "Match.h"
#include "Molecule.h"

#include <fstream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

/* float calculateRMSD(Molecule<Atom>& origMol, Molecule<Atom>& transMol) {
  float square_rmsd =0;
  if((origMol.size() != transMol.size()) || (origMol.size() == 0)) {
    std::cerr << "different molecules " << origMol.size() << " " << transMol.size()<< std::endl;
    exit(1);
  }
  for(unsigned int i=0; i<origMol.size();i++) {
    square_rmsd+=origMol(i).dist2(transMol(i));
  }
  return sqrt(square_rmsd/origMol.size());
} */

double getTempPercentile(Molecule<Atom> &mol1, Molecule<Atom> &mol2, float percentile) {
    std::vector<double> temps1, temps2;
    for (unsigned int i = 0; i < mol1.size(); i++) {
        temps1.push_back(mol1[i].getTempFactor());
        temps2.push_back(mol2[i].getTempFactor());
    }

    std::sort(temps1.begin(), temps1.end());
    std::sort(temps2.begin(), temps2.end());
    int index1 = std::max((int)(percentile * mol1.size() - 1), 0);
    int index2 = std::max((int)(percentile * mol1.size() - 1), 0);
    return std::min(temps1[index1], temps2[index2]);
}

Match calculateTrans(Molecule<Atom> &origMol, Molecule<Atom> &transMol) {
    if ((origMol.size() != transMol.size()) || (origMol.size() == 0)) {
        std::cerr << "different molecules " << origMol.size() << " " << transMol.size() << std::endl;
        exit(1);
    }
    Match match;
    float tempThreshold = std::min(80.0, getTempPercentile(origMol, transMol, 0.5));
    for (unsigned int i = 0; i < origMol.size(); i++) {
        if (origMol[i].getTempFactor() < tempThreshold || transMol[i].getTempFactor() < tempThreshold)
            continue;
        match.add(i, i);
    }
    match.calculateBestFit(origMol, transMol);
    // std::cout << "Got rmsd: " << match.rmsd() << std::endl;
    return match;
}

Molecule<Atom> readMolecule(std::string molName, bool all_atoms) {
    Molecule<Atom> mol;
    std::ifstream molFile(molName);
    if (!molFile) {
        std::cerr << "Can't open file " << molName << std::endl;
        exit(0);
    }

    if (!all_atoms) {
        mol.readPDBfile(molFile, PDB::CAlphaSelector());
        // try nucleic acids
        if (mol.size() == 0) {
            molFile.clear();
            molFile.seekg(0, std::ios::beg);
            mol.readPDBfile(molFile, PDB::PSelector());
        }
    }

    // in case no CA were read, we use all atoms
    if (all_atoms || mol.size() == 0) {
        molFile.clear();
        molFile.seekg(0, std::ios::beg);
        mol.readAllPDBfile(molFile);
    }
    molFile.close();

    return mol;
}

int main(int argc, char **argv) {
    // output arguments
    for (int i = 0; i < argc; i++)
        std::cerr << argv[i] << " ";
    std::cerr << std::endl;

    bool all_atoms = false;
    po::options_description desc(
        "Usage: AF2trans <receptorRef> <ligandRef> <receptorAF2_1> <ligandAF2_1> <receptorAF2_2> <ligandAF2_2> ...\n "
        "translates AF2 complexes into transformations of the ligandRef onto the receptorRef\n");
    desc.add_options()("help",
                       "AF2mer2trans - produces ligand onto receptor docking like transformations from AF2 models\n")(
        "input-files", po::value<std::vector<std::string>>(), "input files")("all,a",
                                                                             "all atoms rmsd (default = false)");

    po::positional_options_description p;
    p.add("input-files", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    std::vector<std::string> files;
    if (vm.count("input-files")) {
        files = vm["input-files"].as<std::vector<std::string>>();
    }
    if (vm.count("help") || files.size() < 4) {
        std::cout << desc << "\n";
        return 0;
    }
    if (vm.count("all_atoms")) {
        all_atoms = true;
    }

    Molecule<Atom> receptorRef = readMolecule(files[0], all_atoms);
    Molecule<Atom> ligandRef = readMolecule(files[1], all_atoms);

    std::cout.precision(4);
    for (unsigned int i = 2; i < files.size(); i += 2) {
        Molecule<Atom> receptorAF2 = readMolecule(files[i], all_atoms);
        Molecule<Atom> ligandAF2 = readMolecule(files[i + 1], all_atoms);
        Match m1 = calculateTrans(receptorRef, receptorAF2);
        Match m2 = calculateTrans(ligandRef, ligandAF2);
        RigidTrans3 T = m1.rigidTrans() * (!m2.rigidTrans());
        std::cout << i / 2 << " | " << m1.rmsd() << "_" << files[i] << " | " << m2.rmsd() << "_" << files[i + 1]
                  << " | " << T << std::endl;
        //    std::cout << rms << std::endl;
    }

    return 0;
}
