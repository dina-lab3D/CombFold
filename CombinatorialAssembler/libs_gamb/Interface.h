#ifndef INTERFACE_H
#define INTERFACE_H

#include <iostream>
#include <vector>
#include <unordered_map>

/*
CLASS
  Interface

  Defines an interface container and used to store interfaces between
  molecules.


KEYWORDS
  Interface, Molecule, Particle

AUTHORS
  Inbal Halperin (inbalhal@math.tau.ac.il),
  Dina Duhovny (duhovka@math.tau.ac.il) and
  Yuval Inbar (inbaryuv@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI> 10/01/05 Dina
     in addAdjacency function: if the adjacency already exists, we now check
     if we can update the distance. if the new distance is smaller, we update it.
     setDistance was added to Adjacency for this purpose.
     buildHigherLevel function was re-written, to iterate over interface instead of group1
<LI> 25/12/03 Dina
     added union and intersection functions for neighbours of two particles
     added adjacencies function that returns ParticleAdjacency for particle,
     instead of previous option that copied them into vector.
</UL

GOALS
  The Interface class saves as a container of the interacting pairs between
  two molecules. Those can be atoms, residues or any other particles.
  The pairs are stored in a hash table, allowing fast queries. Not that the
  interface is not symmetric, it is one sided. The queries of a type
  all_neighbours_of_particle are fast for particles of first molecule only
  since the hashing is not symmetric. So if quieries of type
  all_neighbours_of_particle for both interface sides are required, it is
  necessary to construct two interface objects, one hashing (mol1_index,
  mol2_index) pairs and the other (mol2_index, mol1_index).

USAGE
  Interface face1;
  face1.addAdjacency(mol_index1, mol_index2, 3.0);
  bool neighbour = face1.isAdjacent(mol_index1, mol_index2);

END

*/
class Interface {
 public:
  // GROUP: construction methods.

  //// Construct a new interface.
  Interface();

  //// add a pair of adjacent particles to the interface
  void addAdjacency(unsigned int index1, unsigned int index2, float dist);

  //// create high level interface: this function is usefull when we are given
  //   atomic interface and want to get residue interface instead
  void buildHigherLevel(Interface& highLevelInterface,
			const std::vector<unsigned int>& group1,
			const std::vector<unsigned int>& group2);

  //// Class that stores a pair of adjacent particles and their distance
  class Adjacency {
  public:
    Adjacency() {}
    Adjacency(unsigned int p1, unsigned int p2, float d):
      particle1(p1), particle2(p2), distance(d) {}
    Adjacency(const Interface::Adjacency& a):
      particle1(a.particle1), particle2(a.particle2), distance(a.distance) {}
    unsigned int first() const { return particle1; }
    unsigned int second() const { return particle2; }
    float dist() const { return distance; }
    void setDistance(float dist) { distance = dist; }
    friend std::ostream& operator<<(std::ostream& s, const Adjacency& adjacency) {
      s << adjacency.first() << ' ' << adjacency.second() << ' ' <<adjacency.dist() << std::endl;
      return s;
    }
  protected:
    unsigned int particle1, particle2;
    float distance;
  };

  // GROUP: queries

  //// Return true if the two particles are neighbours
  bool isAdjacent(unsigned int index1, unsigned int index2) const;

  //// Return true if the particles belongs to interface
  bool isInterface(unsigned int index1) const;

  //// Returns number of particles in interface
  unsigned int size() const;

  //// Returns number of adjacencies (couples of particles in interface).
  unsigned int adjacenciesNumber() const { return adjacencyCounter; }

  //// Returns number of neighbours of a given particle
  unsigned int neighboursNumber(unsigned int index1) const;

  //// return all the neighbours of a given particle
  void adjacencies(std::vector<Adjacency>& vAdj, unsigned int index1) const;

  typedef std::unordered_map<int, Adjacency> ParticleAdjacency;//int = index2
  typedef std::unordered_map<int, ParticleAdjacency> InterfaceAdjacency;//int = index1

  //// return ParticleAdjacency for a given index
  const ParticleAdjacency& adjacencies(unsigned int index1) const;

  //// return the intersection of neighbours of two particle adjacencies
  void neighboursIntersection(const ParticleAdjacency& pa1, const ParticleAdjacency& pa2,
			      std::vector<unsigned int>& intersection) const;

  //// return the intersection of neighbours of two indices
  void neighboursIntersection(unsigned int index1, unsigned int index2,
			      std::vector<unsigned int>& intersection) const;

  //// return the union of neighbours of two particle adjacencies
  void neighboursUnion(const ParticleAdjacency& pa1, const ParticleAdjacency& pa2,
		       std::vector<unsigned int>& nUnion) const;

  //// return the union of neighbours of two indices
  void neighboursUnion(unsigned int index1, unsigned int index2,
		       std::vector<unsigned int>& nUnion) const;

  //// Interface iterator
  class iterator {
  public:
    iterator();
    iterator(InterfaceAdjacency::iterator iInterfaceAdjIt,
	     ParticleAdjacency::iterator iParticleAdjIt,
	     InterfaceAdjacency::iterator iInterfaceAdjacencyEnd);
    iterator operator++();
    const Adjacency& operator*() const;
    bool operator==(const Interface::iterator& it1) const;
    bool operator!=(const Interface::iterator& it1) const;
  protected:
    InterfaceAdjacency::iterator interfaceAdjIt;
    ParticleAdjacency::iterator particleAdjIt;
    InterfaceAdjacency::iterator interfaceAdjacencyEnd;
  };

  iterator begin();
  iterator end();

 protected:
  //// Interface container
  InterfaceAdjacency interfaceAdjacency;
  //// Number of pairs in interface
  unsigned int adjacencyCounter;
};


#endif
