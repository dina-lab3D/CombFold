#include <stdlib.h>
#include <stdio.h>
#include <cctype>

#include "CIF.h"
#include <string.h>
#include <boost/algorithm/string.hpp>


std::map<std::string, unsigned int> CIF::atomSiteTable_;

bool CIF::readAtomSiteTable(std::istream& ifile)
{
  unsigned int index = 0;
  int place = 0;
  atomSiteTable_.clear();
  while (!ifile.eof()) {
    std::string line;
    getline(ifile, line);
    if(isAtomSiteLine(line)) {
      boost::trim(line);
      std::string key = line.substr(11);
      atomSiteTable_[key] = index;
      index++;
    } else {
      // we are done reading the atom site dictionary, go back to the beg. of the line
      if(index > 0) {
        ifile.seekg(place);
        break;
      }
    }
    place = ifile.tellg();
  }
  if(index > 0) return true;
  return false;
}

void CIF::writeAtomSiteTable(std::ostream& ofile) {
  ofile << "#\n loop_\n";

  ofile << "_atom_site.group_PDB " << std::endl;
  ofile << "_atom_site.id " << std::endl;
  ofile << "_atom_site.type_symbol " << std::endl;
  ofile << "_atom_site.label_atom_id " << std::endl;
  ofile << "_atom_site.label_alt_id " << std::endl;
  ofile << "_atom_site.label_comp_id " << std::endl;
  ofile << "_atom_site.label_asym_id " << std::endl;
  ofile << "_atom_site.label_entity_id " << std::endl;
  ofile << "_atom_site.label_seq_id " << std::endl;
  ofile << "_atom_site.pdbx_PDB_ins_code " << std::endl;
  ofile << "_atom_site.Cartn_x " << std::endl;
  ofile << "_atom_site.Cartn_y " << std::endl;
  ofile << "_atom_site.Cartn_z " << std::endl;
  ofile << "_atom_site.occupancy " << std::endl;
  ofile << "_atom_site.B_iso_or_equiv " << std::endl;
  ofile << "_atom_site.pdbx_formal_charge " << std::endl;
  ofile << "_atom_site.auth_seq_id " << std::endl;
  ofile << "_atom_site.auth_comp_id " << std::endl;
  ofile << "_atom_site.auth_asym_id " << std::endl;
  ofile << "_atom_site.auth_atom_id " << std::endl;
  ofile << "_atom_site.pdbx_PDB_model_num " << std::endl;
}

bool CIF::isAtomSiteLine(const std::string& CIFline)
{
  return CIFline.substr(0,11) == "_atom_site.";
}

std::string CIF::getColumn(const std::string& CIFline, const std::string& key)
{
  std::vector<std::string> result;
  boost::split(result, CIFline, boost::is_any_of("\t "), boost::token_compress_on);
  if(atomSiteTable_.find(key) != atomSiteTable_.end()) {
    unsigned int col = atomSiteTable_[key];
    if(col < result.size())
      return result[col];
  }
  return std::string();
}

bool CIF::isATOMrec(const std::string& CIFline)
{
  return getColumn(CIFline, "group_PDB") == "ATOM";
}

bool CIF::isHETATMrec(const std::string& CIFline)
{
  return getColumn(CIFline, "group_PDB") == "HETATM";
}

float CIF::atomXCoord(const std::string& CIFline)
{
  try {
    return std::stof(getColumn(CIFline, "Cartn_x"));
  } catch(...) {
    return 0.0;
  }
}

float CIF::atomYCoord(const std::string& CIFline)
{
  try {
    return std::stof(getColumn(CIFline, "Cartn_y"));
  } catch(...) {
    return 0.0;
  }
}

float CIF::atomZCoord(const std::string& CIFline)
{
  try {
    return std::stof(getColumn(CIFline, "Cartn_z"));
  } catch(...) {
    return 0.0;
  }
}

int CIF::atomIndex(const std::string& CIFline)
{
  try {
    return std::stof(getColumn(CIFline, "id"));
  } catch(...) {
    return 0;
  }
}

std::string CIF::atomChainId(const std::string& CIFline)
{
  return getColumn(CIFline, "label_asym_id");
}

std::string CIF::atomChainIdAuthor(const std::string& CIFline)
{
  return getColumn(CIFline, "auth_asym_id");
}

int CIF::atomResidueIndex(const std::string& CIFline)
{
  int ret = 0;
  try {
    ret = std::stoi(getColumn(CIFline, "auth_seq_id"));
  } catch(std::exception& err) {
    try {
      ret = std::stoi(getColumn(CIFline, "label_seq_id"));
    }
    catch(std::exception& err) {
      return ret;
    }
  }
  return ret;
}

std::string CIF::atomType(const std::string& CIFline)
{
  return getColumn(CIFline, "label_atom_id");
}

const std::string CIF::getAtomResidueLongName(const std::string& CIFline)
{
  return getColumn(CIFline, "label_comp_id");
}

AtomEntryType CIF::getAtomEntryType(const std::string& CIFline)
{
  if(isHETATMrec(CIFline))
    return(HETATM);
  else if (isATOMrec(CIFline))
    return(ATOM);
  else
    return(UNK);
}

float CIF::getOccupancy(const std::string& CIFline)
{
  return std::stof(getColumn(CIFline,"occupancy"));
}

float CIF::getTempFactor(const std::string& CIFline)
{
  return std::stof(getColumn(CIFline,"B_iso_or_equiv"));
}
