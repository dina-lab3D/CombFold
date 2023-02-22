#include "Interface.h"

Interface::Interface(): adjacencyCounter(0){}

void Interface::addAdjacency(unsigned int index1, unsigned int index2, float dist){
  //If an exact entry exists (both indexes are found) - update distance if needed
  InterfaceAdjacency::iterator it1=interfaceAdjacency.find(index1);
  if(it1!=interfaceAdjacency.end()) {//index1 exists
     ParticleAdjacency& pAdjacency = it1->second;
     ParticleAdjacency::iterator it2=pAdjacency.find(index2);
     if(it2!= pAdjacency.end()) {//index2 exists
       if((it2->second).dist() > dist) { // update distance
	 (it2->second).setDistance(dist);
       }
       return;
     }
     else{//index2 does not exist
       Adjacency adjacency(index1, index2, dist);
       adjacencyCounter++;
       pAdjacency[index2]= adjacency;
     }
  }
  else{//index1 doesn't exist
    Adjacency adjacency(index1, index2, dist);
    adjacencyCounter++;
    ParticleAdjacency particleAdjacency;
    particleAdjacency[index2]= adjacency;
    interfaceAdjacency [index1]=particleAdjacency;
  }
}

void Interface::buildHigherLevel(Interface& highLevelInterface,
				 const std::vector<unsigned int>& group1,
				 const std::vector<unsigned int>& group2) {
  for(iterator iter = begin(); iter!= end(); ++iter) {
    const Adjacency& adj=*iter;
    int first = group1[adj.first()];
    int second = group2[adj.second()];
    highLevelInterface.addAdjacency(first, second, adj.dist());
  }
}

bool Interface::isAdjacent(unsigned int index1, unsigned int index2) const {
  InterfaceAdjacency::const_iterator interfaceIt=interfaceAdjacency.find(index1);
  if(interfaceIt!=interfaceAdjacency.end()) {
    ParticleAdjacency::const_iterator particleIt = interfaceIt->second.find(index2);
    if(particleIt!=interfaceIt->second.end())
      return true;
  }
  return false;
}

bool Interface::isInterface(unsigned int index1) const {
  InterfaceAdjacency::const_iterator it=interfaceAdjacency.find(index1);
  if(it!=interfaceAdjacency.end())
    return true;
  return false;
}

unsigned int Interface::size() const {
  return interfaceAdjacency.size();
}

unsigned int Interface::neighboursNumber(unsigned int index1) const {
  InterfaceAdjacency::const_iterator it=interfaceAdjacency.find(index1);
  if(it==interfaceAdjacency.end())
    return 0;
  return it->second.size();
}

void Interface::adjacencies(std::vector<Adjacency>& vAdj, unsigned int index1) const {
  InterfaceAdjacency::const_iterator interfaceIt=interfaceAdjacency.find(index1);
  if(interfaceIt != interfaceAdjacency.end()) {
    vAdj.reserve(interfaceIt->second.size());
    for(ParticleAdjacency::const_iterator particleIt = interfaceIt->second.begin();
	particleIt!=interfaceIt->second.end(); particleIt++)
      vAdj.push_back(particleIt->second);
  }
}

const Interface::ParticleAdjacency& Interface::adjacencies(unsigned int index1) const {
  InterfaceAdjacency::const_iterator interfaceIt=interfaceAdjacency.find(index1);
  static const ParticleAdjacency empty;
  if(interfaceIt != interfaceAdjacency.end())
    return interfaceIt->second;
  return empty;
}

void Interface::neighboursIntersection(const Interface::ParticleAdjacency& pa1,
				       const Interface::ParticleAdjacency& pa2,
				       std::vector<unsigned int>& intersection) const {
  if(pa1.size() > pa2.size())
    return neighboursIntersection(pa2, pa1, intersection);
  intersection.reserve(pa1.size());
  for(ParticleAdjacency::const_iterator it = pa1.begin(); it!=pa1.end(); it++)
    if(pa2.find(it->first)!= pa2.end())
      intersection.push_back(it->first);
}

void Interface::neighboursIntersection(unsigned int index1, unsigned int index2,
				       std::vector<unsigned int>& intersection) const {
  neighboursIntersection(adjacencies(index1), adjacencies(index2), intersection);
}

void Interface::neighboursUnion(const Interface::ParticleAdjacency& pa1,
				const Interface::ParticleAdjacency& pa2,
				std::vector<unsigned int>& nUnion) const {
  if(pa1.size() > pa2.size())
    return neighboursUnion(pa2, pa1, nUnion);
  nUnion.reserve(pa1.size() + pa2.size());
  for(ParticleAdjacency::const_iterator it = pa1.begin(); it!=pa1.end(); it++)
    nUnion.push_back(it->first);
  for(ParticleAdjacency::const_iterator it = pa2.begin(); it!=pa2.end(); it++)
    if(pa1.find(it->first) == pa1.end())
      nUnion.push_back(it->first);
}

void Interface::neighboursUnion(unsigned int index1, unsigned int index2,
				std::vector<unsigned int>& nUnion) const {
  neighboursUnion(adjacencies(index1), adjacencies(index2), nUnion);
}

Interface::iterator::iterator(InterfaceAdjacency::iterator iInterfaceAdjIt,
			      ParticleAdjacency::iterator iParticleAdjIt,
			      InterfaceAdjacency::iterator iInterfaceAdjacencyEnd) {
  interfaceAdjIt=iInterfaceAdjIt;
  particleAdjIt=iParticleAdjIt;
  interfaceAdjacencyEnd=iInterfaceAdjacencyEnd;
}

Interface::iterator Interface::iterator::operator++() {
  if ( ++particleAdjIt == interfaceAdjIt->second.end() )//promote line index. if end of line than
    if( (++interfaceAdjIt) != interfaceAdjacencyEnd )//promote to next line. If not the last line than
      particleAdjIt=interfaceAdjIt->second.begin();//set line index to beginning of line
  return *this;
}

const Interface::Adjacency& Interface::iterator::operator*() const{
  return particleAdjIt->second;
}

inline bool Interface::iterator::operator==(const Interface::iterator& it1) const{
  if(it1.interfaceAdjIt == it1.interfaceAdjacencyEnd &&
     this->interfaceAdjIt == this->interfaceAdjacencyEnd)
    return true;
  else
    return (it1.interfaceAdjIt == this->interfaceAdjIt &&
	    it1.particleAdjIt == this->particleAdjIt );
}

bool Interface::iterator::operator!=(const Interface::iterator& it1) const{
  return !(it1 == *this);
}


Interface::iterator Interface::begin() {
  InterfaceAdjacency::iterator interfaceStart=interfaceAdjacency.begin();
  if ( interfaceAdjacency.empty() ){//interface hash map empty
    ParticleAdjacency::iterator emptyParticleIt=ParticleAdjacency::iterator();
    return iterator(interfaceStart, emptyParticleIt, interfaceStart );
  }
  return iterator( interfaceStart, (interfaceStart->second).begin(), interfaceAdjacency.end() );
}

Interface::iterator Interface::end() {
  InterfaceAdjacency::iterator interfaceEnd=interfaceAdjacency.end();
  ParticleAdjacency::iterator emptyParticleIt=ParticleAdjacency::iterator();
  return iterator( interfaceEnd, emptyParticleIt, interfaceEnd );
}
