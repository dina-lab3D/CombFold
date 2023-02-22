#ifndef _Match_h
#define _Match_h

#include "RigidTrans3.h"
#include <vector>
#include <algorithm>

/*
CLASS
  Match

  Defines a match object used to build and refine matches between Molecule
  type objects. Matches are build by matching pairs of Particles from each
  Molecule.

KEYWORDS
  Match, Molecule, RMSD, linear, algebra, rigid

AUTHORS
  Meir Fuchs. (meirfux@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI>
</UL

GOALS
  The Match class was designed for match building, handling and reefinement.
  A match assumes two Molecule objects. Using specified pairs of Particles
  taken from the two Molecules the match is built. A Match has a score, an
  RMS score and a rigid transformation used by the match. This rigid
  transformation may be refined and changed by finding the best rigid
  transformation to match the given pairs of particles stored in the Match.

  Special attention should be given to the add method. This method lets the
  user add a pair to the match. Besides specifying the indices of the Particles
  within the Molcules the add method recieves a score parameter and a
  priority parameter. The score parameter specifies how much the Particle
  pair contributes to the score. The priority tells the add method how
  "important" the pair is. If a different pair is added, one that clashes
  with the first pair (i.e. the particles overlap) add has to decide which
  pair to bump off). The descision is made using the highest priority. An
  example of a priority function would the inverse distance. In this case
  add will prefer the closer pair even if the pair with the longer distance
  scores better.

  Becasue of the difficulty in determinig how a match should be constucted and
  built the non-public members of this class were left protected. This allows
  a user to write Match descendants that best fit his needs.

USAGE
  The following program is a demo utilizing most of the methods supported by
  the Match object. Given two Molecule<Particle> objects, model and scene this
  program will randomly match pairs of particles from the molecules using the
  inverse distance for both score and priority. Every 200 attempts to
  add random pairs to the match the program recalculates the best linear
  transformation. This process is repeated 1000 times for a total of 200,000
  attempted pair additions.
  EXAMPLE
    Match M;
    for (int i=0; i<1000; i++) {
      for (int j=0; j<200; j++) {
        // randomly pick pair
        unsigned int modelParticle = model.size()*drand48();
	unsigned int sceneParticle = scene.size()*drand48();
	// score is pair distance inverse.
	float score = (1 / (1+(model[modelParticle] | scene[sceneParticle])));
	// add pair according to score
	M.add(modelParticle, sceneParticle, score ,score);
      } // j
      M.calculateBestFit(model, scene);  // every 200 random pairs recalc
      cout << M.rigidTrans() << '\n';    // rig. trans. and print it.
    } // i

    int matchingParticle;
    for(int i=0; i<model.size(); i++)
      if ((matchingParticle = M.sceneParticle(i)) >= 0)
        cout << i << '\t' << matchingParticle << '\n';
  END

  Comment: Generally, it is expected (By me at least)
  that if a very large linear match infact exists between model and scene
  (e.g. try model=scene) the calculated rigid transformation will
  probablistically tend to such a match.

*/
class Match
{
public:

  //// ParticleMatch is used by the class implementation to store matching
  // pairs within an STL vector.
  class ParticleMatch
  {
  public:
    ParticleMatch();
    ParticleMatch(unsigned int model, unsigned int scene,
		  float sc, float pri);

    friend std::ostream& operator<<(std::ostream& s, const ParticleMatch& pm);
    unsigned int model;
    unsigned int scene;
    float score;
    float priority;
  };

  // A vector of ParticleMatch is used to store the matching pairs of
  // particles with their score and priority
  typedef std::vector<ParticleMatch> MatchList;


  // GROUP: Contructors

  //// Construct a new match.
  Match();

  //// Contruct a new Match object using the candidate rigid transformation
  // used to match the two sets of points. This rigid transformation may be
  // further enhanced after a set of matching pairs is added using the
  // calculateBestFit method.
  Match(const RigidTrans3& rt);

  //// Match copy contructor. Implicitely declared to make sure that all
  // data referred to by the Match object is actually copied.
  Match(const Match& m);

  // GROUP: Match construction and refinement methods.

  //// Add a match pair. The pair is defined using two integers (possibly
  // indices into some array of objects), a match score and a priority. The
  // match score defines how much matching the model particle to the scene
  // particle contributes to the total match score. This score should
  // reflect "how well the pair of particles fit together". A second
  // parameter sent with the pair is their priority. Two particles may
  // have a high score but matching one of the particles with a different
  // particle might be "more correct". For example one might use inverse
  // distance as a priority - prefering particles that are closer and use
  // some bio-chemical properties for score (e.g. hydrophobicity, residue
  // match using PAM scores etc.).
  // The function returns true if the pair of particles was actually added.
  bool add(const unsigned int modelParticle, const unsigned int sceneParticle,
           const float score,
           const float priority);

  //// This version of add does not use priority and does not check for double
  // references to the same particles in the match list. The user may prefer
  // to overlook doubly referenced particles in many cases. This method is
  // also faster because no search is performed on the existing pairs.
  void add(const unsigned int modelParticle, const unsigned int sceneParticle,
           const float score = 1.0);

  //// Adds ParticleMatch object
  void add(const ParticleMatch& pm);

  //// Calculate the total score of the match using the pair scores given
  // for each pair in the match (See the add method).
  float calculateTotalScore();


  // GROUP: Inspection methods.

  //// Returns the match size.
  unsigned int size() const;

  //// After calculating the best linear transformation using
  // calculateBestFit this method will return the RMSD of the match.
  float rmsd() const;

  //// After calculating the total score of the match using calculateTotalScore
  // this method will return the total score.
  float totalScore() const;

  //// At any point, the rigid transformation currently used by the match will
  // be returned.
  const RigidTrans3& rigidTrans() const;

  //// Given a model particle returns the matching scene particle. If no such
  // particle is found the method will return a negative value.
  int sceneParticle(const unsigned int model) const;

  //// Given a scene particle returns the matching model particle. If no such
  // particle is found the method will return a negative value.
  int modelParticle(const unsigned int scene) const;

  //// Return information about the pair indexed by index in the match in
  // the form of a ParticleMatch object reference.
  const ParticleMatch& operator[](const unsigned int index) const{return pairs[index];}

  //// Return information about the pair indexed by index in the match in
  // the form of a ParticleMatch object reference.
  ParticleMatch& operator[](const unsigned int index){return pairs[index];}

  //// Calculate the best rigid transformation to minimize the inter-point
  // RMSD using the pair scores as weights. The method accepts the two
  // TMolecules from which the particles are taken from and from which it
  // will get the particles positions.
  // TMolecule should support operator[], like vector<>.
  template<class TMolecule>
  void calculateBestFit(const TMolecule& model, const TMolecule& scene);


  ////
  // Calculates the least-squares error between the two given
  // TMolecules. A negative value is returned in case of an error.
  // <br>
  // TMolecule is assumed to be a collection of elements with 3D
  // coordinates. Thus, TMolecule should support operator[] that
  // returns an element with a cast Vector3() opertator. Examples for
  // TMolecule are: vector<Vector3> and Molecule<Atom>. <br>
  // Author: Oranit Dror (oranit@tau.ac.il)
  template<class TMolecule>
  static float calculateLeastSquaresError(const TMolecule& model, const TMolecule& scene);

  ////
  // Calculates the RMSD between the two given TMolecules. A negative
  // is be returned in case of an error. <br>
  // TMolecule is assumed to be a collection of elements with 3D
  // coordinates. Thus, TMolecule should support operator[] that
  // returns an element with a cast Vector3() opertator. Examples for
  // TMolecule are: vector<Vector3> and Molecule<Atom>. <br>
  // Author: Oranit Dror (oranit@tau.ac.il)
  template<class TMolecule>
  static float calculateRMSD(const TMolecule& model, const TMolecule& scene);

  ////
  // Sort the match list according to the given Comperator object. The
  // Comperator object is a Binary Predicate that compares two
  // objects, returning true if the first precedes the second. This
  // predicate must satisfy the standard mathematical definition of a
  // strict weak ordering. The precise requirements are stated in the
  // STL StrictWeakOrdering documentation, but what they roughly mean
  // is that the predicate has to behave the way that "less than"
  // behaves: if a is less than b then b is not less than a, if a is
  // less than b and b is less than c then a is less than c, and so
  // on.
  // Author: Oranit Dror (oranit@tau.ac.il)
  template<class TComperator>
  void sortMatchList(const TComperator& comperator);

  ////
  // Erase the particle pairs of the match list
  void clear() {
    pairs.clear();
  }

  // The rest of this class is left protected in order to allow match to
  // evolve freely.

  friend std::ostream& operator<<(std::ostream& s, const Match& m);

protected:
  MatchList     pairs;

  float         rms;
  float         sc;
  RigidTrans3   trans;

};


template<class TComperator>
void Match::sortMatchList(const TComperator& comperator) {
  sort(pairs.begin(), pairs.end(), comperator);
}

template<class TMolecule>
float Match::calculateLeastSquaresError(const TMolecule& model, const TMolecule& scene) {
  if(model.size() != scene.size())
    return -1;

  float leastSquaresError = 0;
  for(unsigned int i = 0 ; i < model.size() ; i++) {
    leastSquaresError += ((Vector3) model[i]).dist2((Vector3) scene[i]);
  }

  return leastSquaresError;
}


template<class TMolecule>
float Match::calculateRMSD(const TMolecule& model, const TMolecule& scene) {
  float leastSquaresError = calculateLeastSquaresError(model, scene);

  if (leastSquaresError < 0) {
    return leastSquaresError;
  }

  if (model.size() == 0) {
    return 0.0;
  }

  return sqrt(leastSquaresError/model.size());
}

template<class TMolecule>
void Match::calculateBestFit(const TMolecule &model,
			     const TMolecule &scene)
{
    /* Initialized data */

    static double sqrt3 = 1.73205080756888;
    static int ip[9] = { 1,2,4,2,3,5,4,5,6 };
    double* u = new double[9];
    double* t = new double[3];
    int n;

    /* System generated locals */
    int i__1;
    double d__1, d__2, d__3;
    static double equiv_5[6], equiv_11[6], equiv_14[3];


    /* Local variables */
    static double spur, a[9]	/* was [3][3] */, b[9]	/* was [3][3]
	    */, d__;
#define e (equiv_14)
    static double g, h__;
    static int i__, j, k, m, l;
    static double p, r__[9]	/* was [3][3] */;
    static double sigma;
    static double e0;
#define e1 (equiv_14)
#define e2 (equiv_14 + 1)
#define e3 (equiv_14 + 2)
    static int m1, m2, m3;
    static double sqrth, wc, xc[3], yc[3];
#define rr (equiv_5)
#define ss (equiv_11)
#define rr1 (equiv_5)
#define rr2 (equiv_5 + 1)
#define rr3 (equiv_5 + 2)
#define rr4 (equiv_5 + 3)
#define rr5 (equiv_5 + 4)
#define rr6 (equiv_5 + 5)
#define ss1 (equiv_11)
#define ss2 (equiv_11 + 1)
#define ss3 (equiv_11 + 2)
#define ss4 (equiv_11 + 3)
#define ss5 (equiv_11 + 4)
#define ss6 (equiv_11 + 5)
    static double cof, det, cth, sth;

/* **** CALCULATES BEST ROTATION & TRANSLATION BETWEEN TWO VECTOR SETS */
/* **** SUCH THAT U*X+T IS THE BEST APPROXIMATION TO Y. */
/* **** THIS VERSION OF THE ALGORITHM IS OPTIMIZED FOR THREE-DIMENSIONAL
*/
/* **** FLOAT VECTOR SPACE. */
/* **** USE OF THIS ROUTINE IS RESTRICTED TO NON-PROFIT ACADEMIC */
/* **** APPLICATIONS. */
/* **** PLEASE REPORT ERRORS TO */
/* **** PROGRAMMER:  W.KABSCH  MAX-PLANCK-INSTITUTE FOR MEDICAL RESEARCH
*/
/*                            JAHNSTRASSE 29, 6900 HEIDELBERG, FRG. */
/* **** REFERENCES:  W.KABSCH  ACTA CRYST.(1978).A34,827-828 */
/*                  W.KABSCH  ACTA CRYST.(1976).A32,922-923 */

/*  W     - W(M) IS WEIGHT FOR ATOM PAIR  # M                     (GIVEN)
*/
/*  X     - X(I,M) ARE COORDINATES OF ATOM # M IN SET X           (GIVEN)
*/
/*  Y     - Y(I,M) ARE COORDINATES OF ATOM # M IN SET Y           (GIVEN)
*/
/*  N     - N IS NUMBER OF ATOM PAIRS                             (GIVEN)
*/
/*  MODE  - 0:CALCULATE RMS ONLY                                  (GIVEN)
*/
/*          1:CALCULATE RMS,U,T   (TAKES LONGER) */
/*  RMS   - SUM OF W*(UX+T-Y)**2 OVER ALL ATOM PAIRS             (RESULT)
*/
/*  U     - U(I,J) IS   ROTATION  MATRIX FOR BEST SUPERPOSITION  (RESULT)
*/
/*  T     - T(I)   IS TRANSLATION VECTOR FOR BEST SUPERPOSITION  (RESULT)
*/
/*  IER   - 0:NO ERROR; -1:N WAS <2; -2:ILLEGAL WEIGHTS W(M)     (RESULT)
*/

/* -----------------------------------------------------------------------
 */
    /* Parameter adjustments */
    --t;
    u -= 4;
    n = size();

    /* Function Body */
    if (n < 3) {
	return;
    }
    wc = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	xc[i__ - 1] = 0.;
/* L1: */
	yc[i__ - 1] = 0.;
    }
/* DIR$ NEXT SCALAR */
    i__1 = n;
    for (m = 1; m <= i__1; ++m) {
/* L222: */
	if (pairs[m-1].score < (float)0.) {
	    return;
	}
    }
    i__1 = n;
    for (m = 1; m <= i__1; ++m) {
	wc += pairs[m-1].score;
	for (i__ = 1; i__ <= 3; ++i__) {
	    xc[i__ - 1] += pairs[m-1].score * scene[pairs[m-1].scene][i__-1];
/* L2: */
	    yc[i__ - 1] += pairs[m-1].score * model[pairs[m-1].model][i__-1];
	}
    }
    if (wc <= 0.) {
	return;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	xc[i__ - 1] /= wc;
	yc[i__ - 1] /= wc;
	for (j = 1; j <= 3; ++j) {
/* L3: */
	    r__[i__ + j * 3 - 4] = 0.;
	}
    }
    e0 = 0.;
    i__1 = n;
    for (m = 1; m <= i__1; ++m) {
/*         DO 4 I=1,3 */
/* Computing 2nd power */
	d__1 = scene[pairs[m-1].scene][0] - xc[0];
/* Computing 2nd power */
	d__2 = model[pairs[m-1].model][0] - yc[0];
	e0 += pairs[m-1].score * (d__1 * d__1 + d__2 * d__2);
	d__ = pairs[m-1].score * (model[pairs[m-1].model][0] - yc[0]);
/*              DO 4 J=1,3 */
/* 4             R(I,J)=R(I,J)+D*(X(J,M)-XC(J)) */
	r__[0] += d__ * (scene[pairs[m-1].scene][0] - xc[0]);
	r__[3] += d__ * (scene[pairs[m-1].scene][1] - xc[1]);
	r__[6] += d__ * (scene[pairs[m-1].scene][2] - xc[2]);
/* Computing 2nd power */
	d__1 = scene[pairs[m-1].scene][1] - xc[1];
/* Computing 2nd power */
	d__2 = model[pairs[m-1].model][1] - yc[1];
	e0 += pairs[m-1].score * (d__1 * d__1 + d__2 * d__2);
	d__ = pairs[m-1].score * (model[pairs[m-1].model][1] - yc[1]);
	r__[1] += d__ * (scene[pairs[m-1].scene][0] - xc[0]);
	r__[4] += d__ * (scene[pairs[m-1].scene][1] - xc[1]);
	r__[7] += d__ * (scene[pairs[m-1].scene][2] - xc[2]);
/* Computing 2nd power */
	d__1 = scene[pairs[m-1].scene][2] - xc[2];
/* Computing 2nd power */
	d__2 = model[pairs[m-1].model][2] - yc[2];
	e0 += pairs[m-1].score * (d__1 * d__1 + d__2 * d__2);
	d__ = pairs[m-1].score * (model[pairs[m-1].model][2] - yc[2]);
	r__[2] += d__ * (scene[pairs[m-1].scene][0] - xc[0]);
	r__[5] += d__ * (scene[pairs[m-1].scene][1] - xc[1]);
	r__[8] += d__ * (scene[pairs[m-1].scene][2] - xc[2]);
/* L4: */
    }
/* **** CALCULATE DETERMINANT OF R(I,J) */
    det = r__[0] * (r__[4] * r__[8] - r__[7] * r__[5]) - r__[3] * (r__[1] *
	    r__[8] - r__[7] * r__[2]) + r__[6] * (r__[1] * r__[5] - r__[4] *
	    r__[2]);
    sigma = det;
/* **** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R */
    m = 0;
/* DIR$ NEXT SCALAR */
    for (j = 1; j <= 3; ++j) {
	i__1 = j;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ++m;
/* L5: */
	    rr[m - 1] = r__[i__ * 3 - 3] * r__[j * 3 - 3] + r__[i__ * 3 - 2] *
		     r__[j * 3 - 2] + r__[i__ * 3 - 1] * r__[j * 3 - 1];
	}
    }
/* ***************** EIGENVALUES *****************************************
 */
/* **** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0 */
    spur = (*rr1 + *rr3 + *rr6) / 3.;
    cof = (*rr3 * *rr6 - *rr5 * *rr5 + *rr1 * *rr6 - *rr4 * *rr4 + *rr1 * *
	    rr3 - *rr2 * *rr2) / 3.;
    det *= det;
/* **** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y-SPUR */
    d__ = spur * spur;
    h__ = d__ - cof;
    g = spur * (cof * 1.5 - d__) - det * .5;
/* **** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER */
    if (h__ <= d__ * 1e-9) {
	goto L8;
    }
    sqrth = sqrt(h__);
    d__ = -g / (h__ * sqrth);
    if (d__ > .9999999) {
	goto L6;
    }
    if (d__ < -.9999999) {
	goto L7;
    }
/* .....HANDLE CASE OF THREE DISTINCT ROOTS */
    d__ = acos(d__) / 3.;
    cth = sqrth * cos(d__);
    sth = sqrth * sqrt3 * sin(d__);
    *e1 = spur + cth + cth;
    *e2 = spur - cth + sth;
    *e3 = spur - cth - sth;
    if (*e3 < 0.) {
	*e3 = 0.;
    }
    m1 = 3;
    m2 = 1;
    m3 = 2;
    goto L10;
/* .....HANDLE SPECIAL CASE OF TWO IDENTICAL ROOTS */
L6:
    *e1 = spur + sqrth + sqrth;
    *e2 = spur - sqrth;
    *e3 = *e2;
    m = 1;
    m1 = 3;
    m2 = 1;
    m3 = 2;
    goto L20;
L7:
    *e1 = spur + sqrth;
    *e2 = *e1;
    *e3 = spur - sqrth - sqrth;
    if (*e3 < 0.) {
	*e3 = 0.;
    }
    m = 3;
    m1 = 1;
    m2 = 2;
    m3 = 3;
    goto L20;
/* .....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS */
L8:
    *e1 = spur;
    *e2 = spur;
    *e3 = spur;

    goto L26;

/* **************** EIGENVECTORS *****************************************
 */
/* .....EIGENVECTORS IN CASE OF THREE DISTINCT ROOTS */
/* DIR$ NEXT SCALAR */
L10:
    for (l = 1; l <= 2; ++l) {
	d__ = e[l - 1];
	*ss1 = (d__ - *rr3) * (d__ - *rr6) - *rr5 * *rr5;
	*ss2 = (d__ - *rr6) * *rr2 + *rr4 * *rr5;
	*ss3 = (d__ - *rr1) * (d__ - *rr6) - *rr4 * *rr4;
	*ss4 = (d__ - *rr3) * *rr4 + *rr2 * *rr5;
	*ss5 = (d__ - *rr1) * *rr5 + *rr2 * *rr4;
	*ss6 = (d__ - *rr1) * (d__ - *rr3) - *rr2 * *rr2;
	j = 1;
	if (tabs(*ss1) >= tabs(*ss3)) {
	    goto L12;
	}
	j = 2;
	if (tabs(*ss3) >= tabs(*ss6)) {
	    goto L13;
	}
L11:
	j = 3;
	goto L13;
L12:
	if (tabs(*ss1) < tabs(*ss6)) {
	    goto L11;
	}
L13:
	d__ = 0.;
	j = (j - 1) * 3;
	for (i__ = 1; i__ <= 3; ++i__) {
	    k = ip[i__ + j - 1];
	    a[i__ + l * 3 - 4] = ss[k - 1];
/* L14: */
	    d__ += ss[k - 1] * ss[k - 1];
	}
	d__ = sqrt(d__);
	for (i__ = 1; i__ <= 3; ++i__) {
/* L15: */
	    a[i__ + l * 3 - 4] /= d__;
	}
    }
L16:
    a[m1 * 3 - 3] = a[m2 * 3 - 2] * a[m3 * 3 - 1] - a[m3 * 3 - 2] * a[m2 * 3
	    - 1];
    a[m1 * 3 - 2] = a[m2 * 3 - 1] * a[m3 * 3 - 3] - a[m3 * 3 - 1] * a[m2 * 3
	    - 3];
    a[m1 * 3 - 1] = a[m2 * 3 - 3] * a[m3 * 3 - 2] - a[m3 * 3 - 3] * a[m2 * 3
	    - 2];
    goto L30;
/* .....EIGENVECTORS IN CASE OF TWO DISTINCT ROOTS */
L20:
    p = 0.;
/* ccccccccccccccc      H=SPUR+G */
    h__ = *e2;
/* DIR$ NEXT SCALAR */
    for (i__ = 1; i__ <= 3; ++i__) {
	k = (i__ * i__ + i__) / 2;
	d__ = (d__1 = rr[k - 1] - h__, tabs(d__1));
	if (d__ < p) {
	    goto L21;
	}
	j = i__;
	p = d__;
L21:
	;
    }
    p = 0.;
    d__ = 0.;
    l = (j - 1) * 3;
    for (i__ = 1; i__ <= 3; ++i__) {
	k = ip[i__ + l - 1];
	a[i__ + 2] = 1.;
	if (i__ != j) {
	    goto L22;
	}
	a[i__ + m * 3 - 4] = rr[k - 1] - h__;
	goto L23;
L22:
	a[i__ + m * 3 - 4] = rr[k - 1];
	p -= a[i__ + m * 3 - 4];
L23:
/* Computing 2nd power */
	d__1 = a[i__ + m * 3 - 4];
	d__ += d__1 * d__1;
    }
    a[j + 2] = p / a[j + m * 3 - 4];
    d__ = sqrt(d__);
/* Computing 2nd power */
    d__1 = a[3];
/* Computing 2nd power */
    d__2 = a[4];
/* Computing 2nd power */
    d__3 = a[5];
    p = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    for (i__ = 1; i__ <= 3; ++i__) {
	a[i__ + 2] /= p;
/* L24: */
	a[i__ + m * 3 - 4] /= d__;
    }
    goto L16;
/* .....EIGENVECTORS IN CASE OF THREE IDENTICAL ROOTS */
L26:
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    d__ = 0.;
	    if (i__ == j) {
		d__ = 1.;
	    }
/* L27: */
	    a[i__ + j * 3 - 4] = d__;
	}
    }
/* **** CALCULATE ROTATION MATRIX */
/* DIR$ NEXT SCALAR */
L30:
    for (l = 1; l <= 2; ++l) {
	d__ = 0.;
	for (i__ = 1; i__ <= 3; ++i__) {
	    b[i__ + l * 3 - 4] = r__[i__ - 1] * a[l * 3 - 3] + r__[i__ + 2] *
		    a[l * 3 - 2] + r__[i__ + 5] * a[l * 3 - 1];
/* L31: */
/* Computing 2nd power */
	    d__1 = b[i__ + l * 3 - 4];
	    d__ += d__1 * d__1;
	}
	d__ = sqrt(d__);
	for (i__ = 1; i__ <= 3; ++i__) {
/* L32: */
	    b[i__ + l * 3 - 4] /= d__;
	}
    }
    b[6] = b[1] * b[5] - b[4] * b[2];
    b[7] = b[2] * b[3] - b[5] * b[0];
    b[8] = b[0] * b[4] - b[3] * b[1];
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L33: */
	    u[i__ + j * 3] = b[i__ - 1] * a[j - 1] + b[i__ + 2] * a[j + 2] +
		    b[i__ + 5] * a[j + 5];
	}
    }
/* **** CALCULATE TRANSLATION VECTOR */
    for (i__ = 1; i__ <= 3; ++i__) {
/* L34: */
	t[i__] = yc[i__ - 1] - u[i__ + 3] * xc[0] - u[i__ + 6] * xc[1] - u[
		i__ + 9] * xc[2];
    }
/* **** CALCULATE RMS ERROR */
    d__ = sqrt(*e3);
    if (sigma < (float)0.) {
	d__ = -d__;
    }
    d__ = d__ + sqrt(*e2) + sqrt(*e1);
    rms = (float)sqrt((d__1 = e0 - d__ - d__, tabs(d__1))/n);

    ++t;
    u+= 4;
    trans = RigidTrans3(Matrix3(Vector3((float)u[0], (float)u[3], (float)u[6]),
                                Vector3((float)u[1], (float)u[4], (float)u[7]),
                                Vector3((float)u[2], (float)u[5], (float)u[8])),
                        Vector3((float)t[0], (float)t[1], (float)t[2]));
    delete[] u;
    delete[] t;
/* **** NORMAL EXIT */

#undef ss6
#undef ss5
#undef ss4
#undef ss3
#undef ss2
#undef ss1
#undef rr6
#undef rr5
#undef rr4
#undef rr3
#undef rr2
#undef rr1
#undef ss
#undef rr
#undef e3
#undef e2
#undef e1
#undef e
}


#endif
