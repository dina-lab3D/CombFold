
/* Note that more complete documentation of the original Connolly method
   can be found in QCPE, program 429, e.g. at
   http://atlas.physbio.mssm.edu/~mezei/molmod_core/docs/connolly.html
 */

#include "connolly_surface.h"
#include <numerics.h>

#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>

#include <math.h>
#include <vector>
#include <algorithm>


/* Put GridPoint in a non-anonymous namespace, since g++ 4.2 has a bug which
   results in a compilation failure if it's in the anonymous namespace:
   https://svn.boost.org/trac/boost/ticket/3729
 */
namespace detail {
struct GridPoint {
  GridPoint() : icube(-1), scube(false), sscube(false) {}

  int icube;
  bool scube, sscube;
};
}

namespace {

static const int MAXSPH = 1000;
static const int MAXCIR = 1000;

struct AtomTypeInfo {
  // unit vectors
  std::vector<Vector3> ua;

  // extended vectors
  std::vector<Vector3> eva;
};

struct AtomInfo {
  AtomInfo() : ico(3), skip(false), icuptr(-1) {}

  // Integer cube coordinate
  std::vector<int> ico;

  // Skip from grid
  bool skip;

  // Index of next atom in same cube as this one
  int icuptr;
};

class Cube {
 public:
  void grid_coordinates(const std::vector<Vector3> &CO, float radmax, float rp) {
    // calculate width of cube from maximum atom radius and probe radius
    width_ = 2. * (radmax + rp);

    // minimum atomic coordinates (cube corner)
    Vector3 comin(1000000.0, 1000000.0, 1000000.0);
    for (unsigned n = 0; n < CO.size(); ++n) {
      comin.updateX(std::min(comin[0], CO[n][0]));
      comin.updateY(std::min(comin[1], CO[n][1]));
      comin.updateZ(std::min(comin[2], CO[n][2]));
    }

    dim_ = 0;
    atom_info_.resize(CO.size());
    for (unsigned n = 0; n < CO.size(); ++n) {
      for (int k = 0; k < 3; ++k) {
        int ico = static_cast<int>((CO[n][k] - comin[k]) / width_);
        /*if(ico >= 0) {
          std::cerr << "cube coordinate out of range" << std::endl;
          }*/
        dim_ = std::max(dim_, ico);
        atom_info_[n].ico[k] = ico;
      }
    }
    dim_++;

    cube_.resize(boost::extents[dim_][dim_][dim_]);

    for (unsigned n = 0; n < CO.size(); ++n) {
      AtomInfo &ai = atom_info_[n];
      if (!ai.skip) {
        add_atom_to_cube(CO, n);
      }
    }

    update_adjacent();
  }

  bool get_neighbors(int n, const std::vector<Vector3> &CO, float rp,
                     const std::vector<int> &IAT,
                     const std::vector<float> &rtype,
                     std::vector<int> &neighbors) {
    AtomInfo &ai = atom_info_[n];
    neighbors.resize(0);
    if (ai.skip) {
      return false;
    }
    int ici = ai.ico[0];
    int icj = ai.ico[1];
    int ick = ai.ico[2];
    if (!cube_[ici][icj][ick].sscube) {
      return false;
    }
    float sumi = 2 * rp + rtype[IAT[n]];
    for (int jck = ick - 1; jck <= ick + 1; ++jck) {
      if (jck >= 0 && jck < dim_) {
        for (int jcj = icj - 1; jcj <= icj + 1; ++jcj) {
          if (jcj >= 0 && jcj < dim_) {
            for (int jci = ici - 1; jci <= ici + 1; ++jci) {
              if (jci >= 0 && jci < dim_) {
                for (int jatom = cube_[jci][jcj][jck].icube; jatom >= 0;
                     jatom = atom_info_[jatom].icuptr) {
                  float sum = sumi + rtype[IAT[jatom]];
                  if (n != jatom && CO[n].dist2(CO[jatom]) < sum * sum) {
                    neighbors.push_back(jatom);
                  }
                }
              }
            }
          }
        }
      }
    }
    return true;
  }

 private:
  void add_atom_to_cube(const std::vector<Vector3> &CO, int n) {
    AtomInfo &ai = atom_info_[n];
    int existing_atom = cube_[ai.ico[0]][ai.ico[1]][ai.ico[2]].icube;
    if (existing_atom < 0) {
      // first atom in this cube
      cube_[ai.ico[0]][ai.ico[1]][ai.ico[2]].icube = n;
      cube_[ai.ico[0]][ai.ico[1]][ai.ico[2]].scube = true;
    } else {
      while (true) {
        double dist = CO[n].dist2(CO[existing_atom]);
        if (dist <= 0.) {
          ai.skip = true;
          std::cerr << "Skipped atom " << n << " with same coordinates as "
                    << existing_atom << std::endl;
          return;
        }
        if (atom_info_[existing_atom].icuptr == -1) {
          atom_info_[existing_atom].icuptr = n;
          return;
        }
        existing_atom = atom_info_[existing_atom].icuptr;
      }
    }
  }

  void update_adjacent() {
    for (int i = 0; i < dim_; ++i) {
      for (int j = 0; j < dim_; ++j) {
        for (int k = 0; k < dim_; ++k) {
          update_cube(i, j, k);
        }
      }
    }
  }

  void update_cube(int i, int j, int k) {
    if (cube_[i][j][k].icube >= 0) {
      for (int i1 = i - 1; i1 <= i + 1; ++i1)
        if (i1 >= 0 && i1 < dim_) {
          for (int j1 = j - 1; j1 <= j + 1; ++j1)
            if (j1 >= 0 && j1 < dim_) {
              for (int k1 = k - 1; k1 <= k + 1; ++k1) {
                if (k1 >= 0 && k1 < dim_ && cube_[i1][j1][k1].scube) {
                  cube_[i][j][k].sscube = true;
                  return;
                }
              }
            }
        }
    }
  }

  int dim_;
  float width_;
  std::vector<AtomInfo> atom_info_;
  boost::multi_array<detail::GridPoint, 3> cube_;
};

struct YonProbe {
  YonProbe(Vector3 py, Vector3 ay)
      : center(py), altitude(ay) {}
  Vector3 center;
  Vector3 altitude;
};

struct SurfPoint {
  bool yon;
  Vector3 s;
  float area;
  int n1, n2, n3;
};

enum PointType {
  VICTIM,
  YON,
  OTHER
};

struct ProbePoint {
  ProbePoint() : type(OTHER) {}

  char ishape;
  Vector3 position;
  Vector3 to_center;
  PointType type;
  std::vector<SurfPoint> points;
};

class YonCube {
 public:
  YonCube(const std::vector<YonProbe> &yon_probes, float rp, float dp,
          float radmax)
      : comin_(1000000.0, 1000000.0, 1000000.0) {
    width_ = 2. * (radmax + rp);
    dp2_ = dp * dp;
    for (unsigned n = 0; n < yon_probes.size(); ++n) {
      comin_.updateX(std::min(comin_[0], yon_probes[n].center[0]));
      comin_.updateY(std::min(comin_[1], yon_probes[n].center[1]));
      comin_.updateZ(std::min(comin_[2], yon_probes[n].center[2]));
      // for (int k = 0; k < 3; ++k) {
      //   comin_[k] = std::min(comin_[k], yon_probes[n].center[k]);
      // }
    }

    dim_ = 0;
    atom_info_.resize(yon_probes.size());
    for (unsigned n = 0; n < yon_probes.size(); ++n) {
      atom_info_[n].ico = get_cube_coordinates(yon_probes[n].center);
      for (int k = 0; k < 3; ++k) {
        dim_ = std::max(dim_, atom_info_[n].ico[k]);
      }
    }
    dim_++;
    cube_.resize(boost::extents[dim_][dim_][dim_]);
    for (int i = 0; i < dim_; ++i) {
      for (int j = 0; j < dim_; ++j) {
        for (int k = 0; k < dim_; ++k) {
          cube_[i][j][k] = -1;
        }
      }
    }
    for (unsigned n = 0; n < yon_probes.size(); ++n) {
      add_probe_to_cube(n);
    }
  }

  bool probe_overlap(const ProbePoint &probe,
                     const std::vector<YonProbe> &yon_probes) {
    std::vector<int> ic = get_cube_coordinates(probe.position);
    for (int k = 0; k < 3; ++k) {
      ic[k] = std::max(ic[k], 0);
      ic[k] = std::min(ic[k], dim_ - 1);
    }
    // Check adjoining cubes
    for (int jci = ic[0] - 1; jci <= ic[0] + 1; ++jci) {
      if (jci >= 0 && jci < dim_) {
        for (int jcj = ic[1] - 1; jcj <= ic[1] + 1; ++jcj) {
          if (jcj >= 0 && jcj < dim_) {
            for (int jck = ic[2] - 1; jck <= ic[2] + 1; ++jck) {
              if (jck >= 0 && jck < dim_) {
                for (int jp = cube_[jci][jcj][jck]; jp >= 0;
                     jp = atom_info_[jp].icuptr) {
                  const YonProbe &yp = yon_probes[jp];
                  if (yp.center.dist2(probe.position) < dp2_ &&
                      yp.altitude * probe.to_center < 0) {
                    return true;
                  }
                }
              }
            }
          }
        }
      }
    }
    return false;
  }

 private:
  void add_probe_to_cube(int n) {
    AtomInfo &ai = atom_info_[n];
    int existing_atom = cube_[ai.ico[0]][ai.ico[1]][ai.ico[2]];
    if (existing_atom < 0) {
      cube_[ai.ico[0]][ai.ico[1]][ai.ico[2]] = n;
    } else {
      while (atom_info_[existing_atom].icuptr != -1) {
        existing_atom = atom_info_[existing_atom].icuptr;
      }
      atom_info_[existing_atom].icuptr = n;
    }
  }

  std::vector<int> get_cube_coordinates(const Vector3 &c) {
    std::vector<int> ret;
    for (int k = 0; k < 3; ++k) {
      ret.push_back(static_cast<int>((c[k] - comin_[k]) / width_));
    }
    return ret;
  }

  Vector3 comin_;
  int dim_;
  float width_;
  float dp2_;
  std::vector<AtomInfo> atom_info_;
  boost::multi_array<int, 3> cube_;
};

// Generate unit vectors over sphere
void genun(std::vector<Vector3> &vec, unsigned n) {
  vec.reserve(n);
  vec.resize(0);

  int nequat = static_cast<int>(std::sqrt(n * 3.14159));
  int nvert = std::max(1, nequat / 2);

  for (int i = 0; i <= nvert; ++i) {
    float fi = (3.14159 * i) / nvert;
    float z = std::cos(fi);
    float xy = std::sin(fi);
    int nhor = std::max(1, static_cast<int>(nequat * xy));
    for (int j = 0; j < nhor; ++j) {
      float fj = (2. * 3.14159 * j) / nhor;
      float x = std::cos(fj) * xy;
      float y = std::sin(fj) * xy;
      if (vec.size() >= n) return;
      vec.push_back(Vector3(x, y, z));
    }
  }
}

// Triple product of three vectors
float det(const Vector3 &a, const Vector3 &b,
          const Vector3 &c) {
  Vector3 ab = a&b; //algebra::get_vector_product(a, b);
  return ab * c;
}

// Make identity matrix
void imatx(std::vector<Vector3> &ghgt) {
  ghgt.resize(3);
  ghgt[0] = Vector3(1., 0., 0.);
  ghgt[1] = Vector3(0., 1., 0.);
  ghgt[2] = Vector3(0., 0., 1.);
}

// Concatenate matrix b into matrix a
void cat(std::vector<Vector3> &a, const std::vector<Vector3> &b) {
  std::vector<Vector3> temp(3);
  /*
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      temp[j][i] = a[0][i] * b[j][0] + a[1][i] * b[j][1] + a[2][i] * b[j][2];
    }
  }
  */

  for (int j = 0; j < 3; ++j) {
    temp[j].updateX(a[0][0] * b[j][0] + a[1][0] * b[j][1] + a[2][0] * b[j][2]);
    temp[j].updateY(a[0][1] * b[j][0] + a[1][1] * b[j][1] + a[2][1] * b[j][2]);
    temp[j].updateZ(a[0][2] * b[j][0] + a[1][2] * b[j][1] + a[2][2] * b[j][2]);
  }

  a = temp;
}

// Conjugate matrix g with matrix h giving ghgt
void conj(const std::vector<Vector3> &h, const std::vector<Vector3> &g,
          std::vector<Vector3> &ghgt) {
  // INITIALIZE GHGT MATRIX TO IDENTITY
  imatx(ghgt);
  // CONCATENATE G H GT
  cat(ghgt, g);
  cat(ghgt, h);
  // CALCULATE GT
  std::vector<Vector3> gt(3);
  /*for (int k = 0; k < 3; ++k) {
    for (int l = 0; l < 3; ++l) {
      gt[l][k] = g[k][l];
    }
    }*/

  for (int l = 0; l < 3; ++l) {
    gt[l].updateX(g[0][l]);
    gt[l].updateY(g[1][l]);
    gt[l].updateZ(g[2][l]);
  }
  cat(ghgt, gt);
}

// Multiply v by a giving w
void multv(const Vector3 &v, const std::vector<Vector3> &a,
           Vector3 &w) {
  //for (int i = 0; i < 3; ++i) {
  //w[i] = a[0][i] * v[0] + a[1][i] * v[1] + a[2][i] * v[2];
    //}
  w.updateX(a[0][0] * v[0] + a[1][0] * v[1] + a[2][0] * v[2]);
  w.updateY(a[0][1] * v[0] + a[1][1] * v[1] + a[2][1] * v[2]);
  w.updateZ(a[0][2] * v[0] + a[1][2] * v[1] + a[2][2] * v[2]);
}

// Return b perpendicular to a
Vector3 vperp(const Vector3 &a) {
  Vector3 b(0., 0., 0.), p;

  // FIND SMALLEST COMPONENT
  float small = 10000.0;
  int m = -1;
  for (int k = 0; k < 3; ++k) {
    if (std::abs(a[k]) < small) {
      small = std::abs(a[k]);
      m = k;
    }
  }
  if(m == 0) b.updateX(1.0);
  if(m == 1) b.updateY(1.0);
  if(m == 2) b.updateZ(1.0);
  //b[m] = 1.0;

  // take projection along a and subtract from b
  float dt = a[m] / a.norm2(); //a.get_squared_magnitude();
  b = b - dt * a;
  // renormalize b
  return b.getUnitVector();
}

// Collision check of probe with neighboring atoms
bool collid(const Vector3 &p, const std::vector<Vector3> &cnbr,
            const std::vector<float> &ernbr, int jnbr, int knbr, int lkf,
            const std::vector<int> &lknbr) {
  for (int i = lkf; i >= 0; i = lknbr[i]) {
    if (i == jnbr || i == knbr) continue;
    float dist2 = p.dist2(cnbr[i]);
    if (dist2 < ernbr[i] * ernbr[i]) {
      return true;
    }
  }
  return false;
}

struct SurfaceInfo {
  SurfaceInfo() : area(0.), npoints(0), nlost_saddle(0), nlost_concave(0) {}

  float area;
  int npoints;
  int nlost_saddle;
  int nlost_concave;
};

void handle_atom(Surface &surface_points, int iatom, float d,
                 const std::vector<int> inbr, const std::vector<Vector3> &CO,
                 float rp, std::vector<bool> &srs, const std::vector<Vector3> &up,
                 std::vector<YonProbe> &yon_probes,
                 const std::vector<Vector3> &circle, const std::vector<int> &IAT,
                 const std::vector<float> &rtype,
                 const std::vector<AtomTypeInfo> &attyp_info,
                 std::vector<ProbePoint> &beforept, SurfaceInfo &surface) {
  float ri = rtype[IAT[iatom]];
  Vector3 ci = CO[iatom];

  /* transfer data from main arrays to neighbors */
  std::vector<Vector3> cnbr;
  std::vector<float> rnbr;
  std::vector<float> ernbr;
  std::vector<float> disnbr;
  std::vector<int> lknbr;
  for (unsigned iuse = 0; iuse < inbr.size(); ++iuse) {
    int jatom = inbr[iuse];
    cnbr.push_back(CO[jatom]);
    rnbr.push_back(rtype[IAT[jatom]]);
    ernbr.push_back(rtype[IAT[jatom]] + rp);
    disnbr.push_back(ci.dist2(CO[jatom]));
    // initialize link to next farthest out neighbor
    lknbr.push_back(-1);
  }

  // Set up a linked list of neighbors in order of
  // increasing distance from iatom
  int lkf = 0;
  // put remaining neighbors in linked list at proper position
  for (int l = lkf + 1; l < static_cast<int>(inbr.size()); ++l) {
    int before = -1;
    int after = lkf;
    // step through the list until before:after bracket the current neighbor
    while (after >= 0 && disnbr[l] < disnbr[after]) {
      before = after;
      after = lknbr[after];
    }
    if (before < 0) {
      lkf = l;
    } else {
      lknbr[before] = l;
    }
    lknbr[l] = after;
  }
  // Handle no-neighbors case
  if (inbr.size() == 0) {
    lkf = -1;
  }

  // medium loop for each neighbor of iatom
  for (unsigned jnbr = 0; jnbr < inbr.size(); ++jnbr) {
    int jatom = inbr[jnbr];
    if (jatom <= iatom) continue;

    float rj = rnbr[jnbr];
    Vector3 cj = cnbr[jnbr];

    /* HERE FOLLOW GEOMETRIC CALCULATIONS OF POINTS, VECTORS AND
       DISTANCES USED FOR PROBE PLACEMENT IN BOTH SADDLE AND
       CONCAVE REENTRANT SURFACE GENERATION */

    /* CALCULATE THE INTERSECTION
       OF THE EXPANDED SPHERES OF IATOM AND JATOM
       THIS CIRCLE IS CALLED THE SADDLE CIRCLE
       THE PLANE IT LIES IN IS CALLED THE SADDLE PLANE */

    Vector3 vij = cj - ci;

    /* CREATE AN ORTHONORMAL FRAME
       WITH UIJ POINTING ALONG THE INTER-ATOMIC AXIS
       AND Q AND T DEFINING THE SADDLE PLANE */
    float dij = vij.norm();
    if (dij <= 0.) {
      std::cerr << "Atoms " << iatom << " and " << jatom
                << " have the same center" << std::endl;
      continue;
    }
    Vector3 uij = vij.getUnitVector();
    Vector3 q = vperp(uij);
    Vector3 t = uij&q;//algebra::get_vector_product(uij, q);

    // CALCULATE THE SADDLE CIRCLE CENTER AND RADIUS
    float f = 0.5 * (1.0 + ((ri + rp) * (ri + rp) - (rj + rp) * (rj + rp)) /
                               (dij * dij));
    // BASE POINT
    Vector3 bij = ci + f * vij;
    float f1 = ri + rj + 2. * rp;
    f1 = f1 * f1 - dij * dij;
    // SKIP TO BOTTOM OF MIDDLE LOOP IF ATOMS ARE TOO FAR APART
    if (f1 <= 0.0) continue;

    float f2 = dij * dij - (ri - rj) * (ri - rj);
    // SKIP TO BOTTOM OF MIDDLE LOOP IF ONE ATOM INSIDE THE OTHER
    if (f2 <= 0.0) continue;

    // HEIGHT (RADIUS OF SADDLE CIRCLE)
    float hij = std::sqrt(f1 * f2) / (2. * dij);
    // A STARTING ALTITUDE
    Vector3 aij = hij * q;

    // CONCAVE REENTRANT SURFACE

    // GATHER MUTUAL NEIGHBORS OF IATOM AND JATOM
    int mutual = 0;
    std::vector<bool> mnbr(inbr.size());
    for (unsigned knbr = 0; knbr < inbr.size(); ++knbr) {
      float d2 = cj.dist2(cnbr[knbr]);
      float radsum = 2. * rp + rj + rnbr[knbr];
      mnbr[knbr] = d2 < radsum * radsum && knbr != jnbr;
      if (mnbr[knbr]) {
        ++mutual;
      }
    }

    // INNER LOOP FOR EACH MUTUAL NEIGHBOR OF IATOM AND JATOM
    char ishape = 3;
    for (unsigned knbr = 0; knbr < inbr.size(); ++knbr) {
      if (!mnbr[knbr]) continue;

      int katom = inbr[knbr];
      if (katom <= jatom) continue;

      // transfer from neighbor array to katom variables
      float rk = rnbr[knbr];
      Vector3 ck = cnbr[knbr];

      /* CALCULATE INTERSECTION OF EXPANDED SPHERE OF KATOM
         WITH SADDLE PLANE. WE WILL CALL THIS THE KATOM CIRCLE. */

      /* PROJECTION OF VECTOR,
         FROM KATOM TO A POINT ON THE SADDLE PLANE,
         ONTO IATOM-JATOM AXIS,
         IN ORDER TO GET DISTANCE KATOM IS FROM SADDLE PLANE */
      float dk = uij[0] * (bij[0] - ck[0]) + uij[1] * (bij[1] - ck[1]) +
                 uij[2] * (bij[2] - ck[2]);

      // CALCULATE RADIUS OF KATOM CIRCLE
      float rijk = (rk + rp) * (rk + rp) - dk * dk;
      // SKIP CONCAVE CALCULATION IF NO INTERSECTION
      if (rijk <= 0.0) continue;
      rijk = std::sqrt(rijk);

      // CALCULATE CENTER OF KATOM CIRCLE
      Vector3 cijk = ck + dk * uij;

      // CALCULATE INTERSECTION OF THE KATOM CIRCLE WITH THE SADDLE CIRCLE
      Vector3 vijk = cijk - bij;
      float dijk = vijk.norm();

      if (dijk <= 0.0) {
        std::cerr << "Atoms " << iatom << ", " << jatom << ", and " << katom
                  << " have concentric circles" << std::endl;
        continue;
      }
      float f = 0.5 * (1.0 + (hij * hij - rijk * rijk) / (dijk * dijk));
      // BASE POINT BIJK IS ON SYMMETRY PLANE AND SADDLE PLANE
      Vector3 bijk = bij + f * vijk;

      float f1 = (hij + rijk) * (hij + rijk) - dijk * dijk;
      // SKIP TO BOTTOM OF INNER LOOP IF KATOM TOO FAR AWAY
      if (f1 <= 0.0) continue;

      float f2 = dijk * dijk - (hij - rijk) * (hij - rijk);
      // SKIP TO BOTTOM OF INNER LOOP IF KATOM CIRCLE INSIDE SADDLE CIRCLE
      // OR VICE-VERSA
      if (f2 <= 0.0) continue;

      float hijk = std::sqrt(f1 * f2) / (2. * dijk);
      Vector3 uijk = vijk.getUnitVector();

      // UIJ AND UIJK LIE IN THE SYMMETRY PLANE PASSING THROUGH THE ATOMS
      // SO THEIR CROSS PRODUCT IS PERPENDICULAR TO THIS PLANE
      Vector3 aijk0 = uij&uijk; //algebra::get_vector_product(uij, uijk);

      std::vector<Vector3> aijk;
      aijk.push_back(aijk0 * hijk);
      aijk.push_back(-aijk0 * hijk);

      // PROBE PLACEMENT AT ENDS OF ALTITUDE VECTORS
      std::vector<Vector3> pijk(2);
      std::vector<bool> pair(2);
      for (int ip = 0; ip < 2; ++ip) {
        pijk[ip] = bijk + aijk[ip];
        // COLLISION CHECK WITH MUTUAL NEIGHBORS
        pair[ip] = !collid(pijk[ip], cnbr, ernbr, jnbr, knbr, lkf, lknbr);
      }
      // IF NEITHER PROBE POSITION IS ALLOWED, SKIP TO BOTTOM OF INNER LOOP
      if (!pair[0] && !pair[1]) continue;
      bool both = pair[0] && pair[1];
      // SOME REENTRANT SURFACE FOR ALL THREE ATOMS
      srs[iatom] = srs[jatom] = srs[katom] = true;

      // GENERATE SURFACE POINTS
      float area = (4. * pi * rp * rp) / up.size();
      for (int ip = 0; ip < 2; ++ip) {
        if (!pair[ip]) continue;

        // DETERMINE WHETHER PROBE HAS SURFACE ON FAR SIDE OF PLANE
        bool yonprb = hijk < rp && !both;
        // CALCULATE VECTORS DEFINING SPHERICAL TRIANGLE
        // THE VECTORS ARE GIVEN THE PROBE RADIUS AS A LENGTH
        // ONLY FOR THE PURPOSE OF MAKING THE GEOMETRY MORE CLEAR
        Vector3 vpi = (ci - pijk[ip]) * rp / (ri + rp);
        Vector3 vpj = (cj - pijk[ip]) * rp / (rj + rp);
        Vector3 vpk = (ck - pijk[ip]) * rp / (rk + rp);
        float sign = det(vpi, vpj, vpk);

        ProbePoint probe_point;
        // GATHER POINTS ON PROBE SPHERE LYING WITHIN TRIANGLE
        for (unsigned i = 0; i < up.size(); ++i) {
          SurfPoint sp;
          // IF THE UNIT VECTOR IS POINTING AWAY FROM THE SYMMETRY PLANE
          // THE SURFACE POINT CANNOT LIE WITHIN THE INWARD-FACING TRIANGLE
          if (up[i] * aijk[ip] > 0.) continue;
          if (sign * det(up[i], vpj, vpk) < 0.) continue;
          if (sign * det(vpi, up[i], vpk) < 0.) continue;
          if (sign * det(vpi, vpj, up[i]) < 0.) continue;
          // CALCULATED WHETHER POINT IS ON YON SIDE OF PLANE
          sp.yon = aijk[ip] * (aijk[ip] + up[i]) < 0.;
          // OVERLAPPING REENTRANT SURFACE REMOVAL
          // FOR SYMMETRY-RELATED PROBE POSITIONS
          if (sp.yon && both) continue;

          // CALCULATE COORDINATES OF SURFACE POINT
          sp.s = pijk[ip] + up[i] * rp;

          // FIND THE CLOSEST ATOM AND PUT THE THREE ATOM NUMBERS
          // IN THE PROPER ORDER
          // N1 IS CLOSEST, N2 < N3
          float dsi = sp.s.dist(ci) - ri;
          float dsj = sp.s.dist(cj) - rj;
          float dsk = sp.s.dist(ck) - rk;
          if (dsi <= dsj && dsi <= dsk) {
            sp.n1 = iatom;
            sp.n2 = jatom;
            sp.n3 = katom;
          } else if (dsj <= dsi && dsj <= dsk) {
            sp.n1 = jatom;
            sp.n2 = iatom;
            sp.n3 = katom;
          } else {
            sp.n1 = katom;
            sp.n2 = iatom;
            sp.n3 = jatom;
          }
          sp.area = area;
          probe_point.points.push_back(sp);
        }
        if (probe_point.points.size() > 0) {
          probe_point.ishape = ishape;
          probe_point.position = pijk[ip];
          probe_point.to_center = aijk[ip];
          probe_point.type = yonprb ? YON : OTHER;
          beforept.push_back(probe_point);

          // SAVE PROBE IN YON PROBE ARRAYS
          if (yonprb) {
            yon_probes.push_back(YonProbe(pijk[ip], aijk[ip]));
          }
        }
      }
    }

    // SADDLE-SHAPED REENTRANT
    ishape = 2;

    // SPECIAL CHECK FOR BURIED TORI

    /* IF NEITHER ATOM HAS ANY REENTRANT SURFACE SO FAR
       (AFTER TRIANGLES WITH ALL KATOMS HAVE BEEN CHECKED)
       AND IF THERE IS SOME MUTUAL NEIGHBOR IN THE SAME MOLECULE
       CLOSE ENOUGH SO THAT THE TORUS CANNOT BE FREE,
       THEN WE KNOW THAT THIS MUST BE A BURIED TORUS */
    if (!srs[iatom] && !srs[jatom] && mutual > 0) {
      bool buried_torus = false;
      for (unsigned knbr = 0; knbr < inbr.size() && !buried_torus; ++knbr) {
        if (!mnbr[knbr]) continue;
        float d2 = bij.dist2(cnbr[knbr]);
        float rk2 = ernbr[knbr] * ernbr[knbr] - hij * hij;
        if (d2 < rk2) {
          buried_torus = true;
        }
      }
      if (buried_torus) continue;
    }

    // CALCULATE NUMBER OF ROTATIONS OF PROBE PAIR,
    // ROTATION ANGLE AND ROTATION MATRIX
    float rij = ri / (ri + rp) + rj / (rj + rp);
    float avh = (std::abs(hij - rp) + hij * rij) / 3.;
    int nrot = std::max(static_cast<int>(std::sqrt(d) * pi * avh), 1);
    float angle = pi / nrot;

    // SET UP ROTATION MATRIX AROUND X-AXIS
    std::vector<Vector3> h(3, Vector3(0., 0., 0.));
    h[0].updateX(1.0);
    h[1].updateY(std::cos(angle));
    h[2].updateZ(std::cos(angle));
    h[1].updateZ(std::sin(angle));
    h[2].updateY(-h[1][2]);
    // CALCULATE MATRIX TO ROTATE X-AXIS ONTO IATOM-JATOM AXIS
    std::vector<Vector3> g;
    g.push_back(uij);
    g.push_back(q);
    g.push_back(t);

    // MAKE THE PROBE PAIR ROTATION MATRIX BE ABOUT THE IATOM-JATOM AXIS
    std::vector<Vector3> ghgt;
    conj(h, g, ghgt);

    // ARC GENERATION
    Vector3 pij = bij + aij;
    Vector3 vpi = (ci - pij) * rp / (ri + rp);
    Vector3 vpj = (cj - pij) * rp / (rj + rp);

    // ROTATE CIRCLE ONTO IATOM-JATOM-PROBE PLANE
    // AND SELECT POINTS BETWEEN PROBE-IATOM AND
    // PROBE-JATOM VECTOR TO FORM THE ARC
    int narc = 0;
    std::vector<bool> ayon;
    std::vector<float> arca;
    std::vector<Vector3> vbs0[2];
    for (unsigned i = 0; i < circle.size(); ++i) {
      Vector3 vps0;
      // ROTATION
      multv(circle[i], g, vps0);
      // IF THE VECTOR IS POINTING AWAY FROM THE SYMMETRY LINE
      // THE SURFACE POINT CANNOT LIE ON THE INWARD-FACING ARC
      if (vps0 * aij > 0.) continue;
      Vector3 vector = vpi&vps0;//algebra::get_vector_product(vpi, vps0);
      if (g[2] * vector < 0.) continue;
      vector = vps0 & vpj; //algebra::get_vector_product(vps0, vpj);
      if (g[2] * vector < 0.) continue;

      /* MAKE ARC POINT VECTORS ORIGINATE WITH SADDLE CIRCLE CENTER BIJ
         RATHER THAN PROBE CENTER BECAUSE THEY WILL BE
         ROTATED AROUND THE IATOM-JATOM AXIS */
      vbs0[0].push_back(vps0 + aij);
      // INVERT ARC THROUGH LINE OF SYMMETRY
      float duij = uij * vbs0[0][narc];
      vbs0[1].push_back(-vbs0[0][narc] + 2 * duij * uij);

      /* CHECK WHETHER THE ARC POINT CROSSES THE IATOM-JATOM AXIS
         AND CALCULATE THE AREA ASSOCIATED WITH THE POINT */
      float ht = aij * vbs0[0][narc] / hij;
      ayon.push_back(ht < 0.);
      arca.push_back((2. * pi * pi * rp * std::abs(ht)) /
                     (circle.size() * nrot));
      narc++;
    }

    // INITIALIZE POWER MATRIX TO IDENTITY
    std::vector<Vector3> pow;
    imatx(pow);

    // ROTATE THE PROBE PAIR AROUND THE PAIR OF ATOMS
    for (int irot = 0; irot < nrot; ++irot, cat(pow, ghgt)) {
      std::vector<Vector3> aijp(2);
      // MULTIPLY ALTITUDE VECTOR BY POWER MATRIX
      multv(aij, pow, aijp[0]);
      // SET UP OPPOSING ALTITUDE
      aijp[1] = -aijp[0];

      // SET UP PROBE SPHERE POSITIONS
      std::vector<Vector3> pijp(2);
      bool pair[2];
      for (int ip = 0; ip < 2; ++ip) {
        pijp[ip] = bij + aijp[ip];
        // CHECK FOR COLLISIONS WITH NEIGHBORING ATOMS
        pair[ip] = !collid(pijp[ip], cnbr, ernbr, jnbr, -1, lkf, lknbr);
      }

      // NO SURFACE GENERATION IF NEITHER PROBE POSITION IS ALLOWED
      if (!pair[0] && !pair[1]) continue;
      bool both = pair[0] && pair[1];
      // SOME REENTRANT SURFACE FOR BOTH ATOMS
      srs[iatom] = srs[jatom] = true;

      /* SKIP TO BOTTOM OF MIDDLE LOOP IF IATOM AND JATOM
         ARE CLOSE ENOUGH AND THE SURFACE POINT DENSITY IS
         LOW ENOUGH SO THAT THE ARC HAS NO POINTS */
      if (narc <= 0) continue;

      // SURFACE GENERATION
      for (int ip = 0; ip < 2; ++ip) {
        if (!pair[ip]) continue;

        // DETERMINE WHETHER PROBE HAS SURFACE ON FAR SIDE OF LINE
        bool yonprb = hij < rp && !both;
        ProbePoint probe_point;
        // THE SADDLE-SHAPED REENTRANT SURFACE POINTS COME FROM THE ARC
        for (int i = 0; i < narc; ++i) {
          SurfPoint sp;
          /* OVERLAPPING REENTRANT SURFACE REMOVAL
             FOR SYMMETRY-RELATED PROBE POSITIONS */
          if (both && ayon[i]) continue;
          // ROTATE THE ARC FROM THE XY PLANE ONTO THE IATOM-JATOM-PROBE PLANE
          Vector3 vbs;
          multv(vbs0[ip][i], pow, vbs);
          // MAKE COORDINATES RELATIVE TO ORIGIN
          sp.s = bij + vbs;
          // FIND THE CLOSEST ATOM AND SET UP THE ATOM NUMBERS FOR THE POINT
          float dsi = sp.s.dist(ci) - ri;
          float dsj = sp.s.dist(cj) - rj;
          if (dsi <= dsj) {
            sp.n1 = iatom;
            sp.n2 = jatom;
          } else {
            sp.n1 = jatom;
            sp.n2 = iatom;
          }
          sp.n3 = -1;

          // WE'VE GOT A SURFACE POINT
          sp.yon = ayon[i];
          sp.area = arca[i];
          probe_point.points.push_back(sp);
        }
        if (probe_point.points.size() > 0) {
          probe_point.ishape = ishape;
          probe_point.position = pijp[ip];
          probe_point.to_center = aijp[ip];
          probe_point.type = yonprb ? YON : OTHER;
          beforept.push_back(probe_point);
          // SAVE PROBE IN YON PROBE ARRAYS
          if (yonprb) {
            yon_probes.push_back(YonProbe(pijp[ip], aijp[ip]));
          }
        }
      }
    }
  }

  /* IF THE PROBE RADIUS IS GREATER THAN ZERO
     AND IATOM HAS AT LEAST ONE NEIGHBOR, BUT NO REENTRANT SURFACE,
     THEN IATOM MUST BE COMPLETELY INACCESSIBLE TO THE PROBE */
  if (rp > 0. && inbr.size() > 0 && !srs[iatom]) return;

  const AtomTypeInfo &attyp = attyp_info[IAT[iatom]];
  float area = (4. * pi * ri * ri) / attyp.ua.size();

  // CONTACT PROBE PLACEMENT LOOP
  for (unsigned i = 0; i < attyp.ua.size(); ++i) {
    // SET UP PROBE COORDINATES
    Vector3 pipt = ci + attyp.eva[i];
    // CHECK FOR COLLISION WITH NEIGHBORING ATOMS
    if (collid(pipt, cnbr, ernbr, -1, -1, lkf, lknbr)) continue;

    // INCREMENT SURFACE POINT COUNTER FOR CONVEX SURFACE
    surface.npoints++;

    // ADD SURFACE POINT AREA TO CONTACT AREA
    surface.area += area;
    Vector3 outco = ci + ri * attyp.ua[i];
    Vector3 outvec = attyp.ua[i];

    // CONTACT
    //SurfacePoint sp(iatom, -1, -1, outco, area, outvec);
    SurfacePoint sp(outco, outvec, area, iatom);
    surface_points.push_back(sp);
  }
  return;
}

SurfaceInfo generate_contact_surface(
    Surface &surface_points, const std::vector<Vector3> &CO,
    float radmax, float rp, float d, std::vector<YonProbe> &yon_probes,
    const std::vector<int> &IAT, const std::vector<float> &rtype,
    const std::vector<AtomTypeInfo> &attyp_info,
    std::vector<ProbePoint> &beforept) {
  Cube cube;
  cube.grid_coordinates(CO, radmax, rp);

  // set up probe sphere and circle
  std::vector<Vector3> up;
  int nup = static_cast<int>(4. * pi * rp * rp * d);
  nup = std::max(1, nup);
  nup = std::min(MAXSPH, nup);
  genun(up, nup);

  int ncirc = static_cast<int>(2. * pi * rp * std::sqrt(d));
  ncirc = std::max(1, ncirc);
  ncirc = std::min(MAXCIR, ncirc);
  std::vector<Vector3> circle;
  for (int i = 0; i < ncirc; ++i) {
    float fi = (2. * pi * i) / ncirc;
    circle.push_back(
        Vector3(rp * std::cos(fi), rp * std::sin(fi), 0.));
  }

  SurfaceInfo surface;
  std::vector<bool> srs(CO.size(), false);
  for (unsigned i = 0; i < CO.size(); ++i) {
    std::vector<int> itnl;
    if (cube.get_neighbors(i, CO, rp, IAT, rtype, itnl)) {
      std::sort(itnl.begin(), itnl.end());
      handle_atom(surface_points, i, d, itnl, CO, rp, srs, up, yon_probes,
                  circle, IAT, rtype, attyp_info, beforept, surface);
    }
  }
  return surface;
}

void get_victim_probes(const std::vector<YonProbe> &yon_probes,
                       std::vector<ProbePoint> &beforept, float rp,
                       float radmax, std::vector<int> &victims) {
  // NO VICTIM PROBES IF NO YON PROBES
  if (yon_probes.size() == 0) return;

  // Probe diameter
  float dp = 2. * rp;

  YonCube cube(yon_probes, rp, dp, radmax);

  int ivic = 0;
  for (std::vector<ProbePoint>::iterator pit = beforept.begin();
       pit != beforept.end(); ++pit, ++ivic) {
    if (pit->type == YON) continue;

    // CHECK IF PROBE TOO FAR FROM SYMMETRY ELEMENT FOR POSSIBLE OVERLAP
    if (pit->to_center.norm2() > dp * dp) continue;

    // LOOK FOR OVERLAP WITH ANY YON PROBE IN THE SAME MOLECULE
    if (cube.probe_overlap(*pit, yon_probes)) {
      victims.push_back(ivic);
      pit->type = VICTIM;
    }
  }
}

void get_eaten_points(const std::vector<YonProbe> &yon_probes,
                      const std::vector<ProbePoint> &beforept, float dp2,
                      const ProbePoint &probe, unsigned &neat, unsigned &nyeat,
                      const std::vector<int> &victims, std::vector<Vector3> &eat,
                      int &pi) {
  neat = nyeat = 0;
  eat.resize(0);
  if (yon_probes.size() == 0) return;

  pi++;
  // DETERMINE IF PROBE IS A YON OR VICTIM PROBE
  if (probe.type == OTHER) {
    return;
  }

  // CHECK THIS VICTIM OR YON PROBE AGAINST ALL YON PROBES
  for (unsigned j = 0; j < yon_probes.size(); ++j) {
    if (probe.position.dist2(yon_probes[j].center) >= dp2) continue;
    if (probe.to_center * yon_probes[j].altitude >= 0.) continue;

    // THIS YON PROBE COULD EAT SOME OF THE PROBE'S POINTS
    neat++;
    nyeat++;
    eat.push_back(yon_probes[j].center);
  }

  // ONLY YON PROBES CAN HAVE THEIR POINTS EATEN BY VICTIMS
  if (probe.type != YON) return;

  // CHECK THIS YON PROBE AGAINST ALL VICTIM PROBES
  for (unsigned j = 0; j < victims.size(); ++j) {
    const ProbePoint &victim = beforept[victims[j]];
    if (probe.position.dist2(victim.position) >= dp2) continue;
    if (probe.to_center * victim.to_center >= 0.) continue;
    // THIS VICTIM PROBE COULD EAT SOME OF THE PROBE'S POINTS
    neat++;
    eat.push_back(victim.position);
  }
}

void check_eaten_points(Surface &surface_points,
                        const std::vector<YonProbe> &yon_probes,
                        const std::vector<ProbePoint> &beforept, float rp,
                        const std::vector<int> &victims, SurfaceInfo &surface) {
  float rp2 = rp * rp;
  float dp = rp * 2.;
  float dp2 = dp * dp;
  int pi = 0;
  for (std::vector<ProbePoint>::const_iterator pit = beforept.begin();
       pit != beforept.end(); ++pit) {
    unsigned neat, nyeat;
    std::vector<Vector3> eat;
    get_eaten_points(yon_probes, beforept, dp2, *pit, neat, nyeat, victims, eat,
                     pi);

    // READ THE SURFACE POINTS BELONGING TO THE PROBE
    for (std::vector<SurfPoint>::const_iterator sit = pit->points.begin();
         sit != pit->points.end(); ++sit) {
      // CHECK SURFACE POINT AGAINST ALL EATERS OF THIS PROBE
      bool point_eaten = false;
      for (unsigned k = 0; k < eat.size(); ++k) {
        // VICTIM PROBES CANNOT EAT NON-YON POINTS OF YON PROBES
        if (!(pit->type == YON && !sit->yon && k >= nyeat) &&
            eat[k].dist2(sit->s) < rp2) {
          point_eaten = true;
          break;
        }
      }

      if (!point_eaten) {
        Vector3 outvec = (pit->position - sit->s) / rp;
        // REENTRANT SURFACE POINT
        surface.npoints++;
        surface.area += sit->area;

        //SurfacePoint sp(sit->n1, sit->n2, sit->n3, sit->s, sit->area, outvec);
        SurfacePoint sp(sit->s, outvec, sit->area, sit->n1, sit->n2, sit->n3);
        surface_points.push_back(sp);
      } else {
        if (pit->ishape == 2) {
          surface.nlost_saddle++;
        } else {
          surface.nlost_concave++;
        }
      }
    }
  }
}

SurfaceInfo generate_reentrant_surface(Surface &surface_points,
                                       const std::vector<YonProbe> &yon_probes,
                                       std::vector<ProbePoint> &beforept,
                                       float rp, float radmax) {
  SurfaceInfo surface;
  std::vector<int> victims;
  get_victim_probes(yon_probes, beforept, rp, radmax, victims);
  std::cerr << yon_probes.size() << " yon and " << victims.size()
            << " victim probes" << std::endl;
  check_eaten_points(surface_points, yon_probes, beforept, rp, victims,
                     surface);
  return surface;
}

void msdots(Surface &surface_points, float d, float rp,
            const std::vector<float> &rtype, const std::vector<Vector3> &CO,
            const std::vector<int> &IAT) {
  if(rp <= 0) {
    std::cerr << "Negative probe radius: " << rp << std::endl; return;
  }
  float radmax = 0.;

  std::vector<AtomTypeInfo> attyp_info(rtype.size());

  // Read atom types
  for (unsigned n = 0; n < rtype.size(); ++n) {
    if(rtype[n] <= 0) {
      std::cerr << "Negative atom radius: " << n << " " << rtype[n] << std::endl;
      return;
    }

    radmax = std::max(radmax, rtype[n]);
    // number of unit vectors depends on sphere area and input density
    int nua = std::max(
        1, static_cast<int>((4. * pi * rtype[n] * rtype[n]) * d));
    nua = std::min(nua, MAXSPH);

    // create unit vector arrays
    genun(attyp_info[n].ua, nua);
    // compute extended vectors for later probe placement
    attyp_info[n].eva.resize(attyp_info[n].ua.size());
    for (unsigned isph = 0; isph < attyp_info[n].ua.size(); ++isph) {
      attyp_info[n].eva[isph] = (rtype[n] + rp) * attyp_info[n].ua[isph];
    }
  }
  std::cerr << "Calculation of surface of "
            << CO.size() << " atoms, with surface point density " << d
            << " and probe radius " << rp << std::endl;

  std::vector<YonProbe> yon_probes;
  std::vector<ProbePoint> beforept;
  SurfaceInfo contact_surface =
      generate_contact_surface(surface_points, CO, radmax, rp, d, yon_probes,
                               IAT, rtype, attyp_info, beforept);

  SurfaceInfo reentrant_surface = generate_reentrant_surface(
      surface_points, yon_probes, beforept, rp, radmax);


  std::cerr << reentrant_surface.nlost_saddle
            << " saddle and " << reentrant_surface.nlost_concave
            << " concave surface points removed during non-symmetry "
            << "overlapping reentrant surface removal" << std::endl;
  std::cerr << contact_surface.npoints
            << " contact, " << reentrant_surface.npoints
            << " reentrant, and "
            << contact_surface.npoints + reentrant_surface.npoints
            << " total surface points" << std::endl;
  std::cerr << "Contact area: " << contact_surface.area << "; reentrant area: "
            << reentrant_surface.area << "; total area: "
            << contact_surface.area + reentrant_surface.area << std::endl;
}

}  // namespace

Surface get_connolly_surface(const ChemMolecule& molecule, float d, float rp) {
  typedef boost::unordered_map<float, int> M;
  M radii2type;

  std::vector<int> IAT(molecule.size());
  std::vector<Vector3> CO(molecule.size());
  std::vector<float> rvdw;
  for (unsigned int i = 0; i < molecule.size(); ++i) {
    float r = molecule[i].getRadius();
    M::const_iterator it = radii2type.find(r);
    int type;
    if (it == radii2type.end()) {
      type = radii2type.size();
      radii2type[r] = type;
      rvdw.push_back(r);
    } else {
      type = it->second;
    }
    IAT[i] = type;
    CO[i] = molecule[i].position();
  }
  std::cerr << "Number of vdW radius types: " << rvdw.size() << std::endl;

  std::cerr << "Total number of atoms: " << CO.size() << std::endl;

  /* --------- RUN CONNOLLY'S MOLECULAR SURFACE PROGRAM -------- */
  Surface surface_points;
  msdots(surface_points, d, rp, rvdw, CO, IAT);
  return surface_points;
}
