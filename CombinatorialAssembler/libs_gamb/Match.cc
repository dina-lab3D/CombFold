#include <math.h>
#include "Match.h"
#include "macros.h"

typedef Match::ParticleMatch ParticleMatch;
typedef Match::MatchList MatchList;
typedef MatchList::const_iterator const_iterator;
typedef MatchList::iterator iterator;

ParticleMatch::ParticleMatch()
  : model(0), scene(0), score(0.0), priority(0.0) {}
ParticleMatch::ParticleMatch(unsigned int modelParticle,
			     unsigned int sceneParticle,
			     float sc, float pri) :
  model(modelParticle), scene(sceneParticle), score(sc), priority(pri) {}


std::ostream& operator<<(std::ostream& s, const ParticleMatch& pm)
{
  s.setf(std::ios::fixed, std::ios::floatfield);
  s.setf(std::ios::right, std::ios::adjustfield);
  s.precision(3);

  s.width(8);
  s << pm.model;
  s.width(8);
  s << pm.scene;
  s.width(10);
  s << pm.score;
  s.width(10);
  s << pm.priority;
  return s;
}

Match::Match() : rms(0.0), sc(0.0), trans() {}

Match::Match(const RigidTrans3& rt) : rms(0.0), sc(0.0), trans(rt) {}

Match::Match(const Match& m) :
  pairs(m.pairs), rms(m.rms), sc(m.sc), trans(m.trans) {}

// The add method is complicated by the fact that the pair being added may
// clash with a pair that has already been inserted into the match. This
// is handled using the priority parameter. Pairs with higher priorities
// are preferred.
bool Match::add(const unsigned int modelParticle,
                const unsigned int sceneParticle,
                const float score,
                const float priority)
{
  iterator sceneClash = pairs.end();
  iterator modelClash = pairs.end();
  bool canAdd = true;

  // Find if there are any clashes when adding the pair.
  // Developer's Note: I've tried this class out using two binary search
  // trees (STL map.h) instead of a vector. The idea was to speed up the
  // search for clashes using binary trees. I tested the two implementations
  // using random match additions to a set of matches. Besides being at
  // least 4 times more costly in terms of memory the binary tree
  // implementation also ran almost twice as slow! So, if you're going to
  // try it remember I was there first.
  for (iterator it = pairs.begin(); it != pairs.end() && canAdd;
       it++) {
    if ((*it).scene == sceneParticle) {   // check for clashes with scene part.
      sceneClash = it;
      canAdd = (priority > (*sceneClash).priority);
      if (modelClash != pairs.end()) break;  // if both scene and model clashes found break
    }
    else if ((*it).model == modelParticle) { // check for model part. clash
      modelClash = it;
      canAdd = (priority > (*modelClash).priority);
      if (sceneClash != pairs.end()) break;  // if both scene and model clashes found break
    }
  }

  // if pair has higher priority or no clashes were found then add the new
  // pair. if clashes occurred remove old pair or pairs from the data
  // structures
  if (canAdd) {
    if (sceneClash != pairs.end()) {
      // replace the scene clashed pair with the new pair.
      *sceneClash = ParticleMatch(modelParticle, sceneParticle,
				  score, priority);
      // if new pair clashed both with a model particle and a scene particle
      // of different two pairs then then remove one pair.
      if (modelClash != pairs.end()) pairs.erase(modelClash);
    }
    else if (modelClash != pairs.end())
      // only model particle clashed. replace the clashed pair with the
      // new pair.
      *modelClash = ParticleMatch(modelParticle, sceneParticle,
				  score, priority);

    else    // no clashes.
      pairs.push_back(ParticleMatch(modelParticle, sceneParticle,
				       score, priority));
  }
  return canAdd;
}

void Match::add(const unsigned int modelParticle,
                const unsigned int sceneParticle,
                const float score)
{
  pairs.push_back(ParticleMatch(modelParticle, sceneParticle, score, 0.0));
}


void Match::add(const ParticleMatch& pm){
  pairs.push_back(pm);
}

unsigned int Match::size()  const
{
  return pairs.size();
}

float Match::rmsd() const
{
  return rms;
}

float Match::totalScore() const
{
  return sc;
}

const RigidTrans3& Match::rigidTrans() const
{
  return trans;
}

int Match::sceneParticle(const unsigned int modelParticle) const
{
  for (const_iterator i = pairs.begin(); i != pairs.end(); i++)
    if ((*i).model == modelParticle)
      return (*i).scene;
  return -1;
}

int Match::modelParticle(const unsigned int sceneParticle) const
{
  for (const_iterator i = pairs.begin(); i != pairs.end(); i++)
    if ((*i).scene == sceneParticle)
      return (*i).model;
  return -1;
}

float Match::calculateTotalScore()
{
  sc = 0.0;
  for (iterator it = pairs.begin(); it != pairs.end(); it++)
    sc+= (*it).score;
  return sc;
}

std::ostream& operator<<(std::ostream& s, const Match& m)
{
  s << "#model  #scene  #score   #priority" << std::endl;
  for(std::vector<ParticleMatch>::const_iterator it=m.pairs.begin() ;
      it != m.pairs.end(); ++it)
    s << *it << std::endl;
  s << std::endl;
  s << "#The transformation is:"  << std::endl;
  s << m.trans << std::endl;
  s << "#score: " << m.sc  << "#Rms: " << m.rms << std::endl;
  return s;
}
