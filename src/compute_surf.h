/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(surf,ComputeSurf)

#else

#ifndef SPARTA_COMPUTE_SURF_H
#define SPARTA_COMPUTE_SURF_H

#include "compute.h"
#include "surf.h"
#include "hash3.h"

namespace SPARTA_NS {

class ComputeSurf : public Compute {
 public:
  ComputeSurf(class SPARTA *, int, char **);
  ComputeSurf(class SPARTA* sparta) : Compute(sparta) {} // needed for Kokkos
  ~ComputeSurf();
  virtual void init();
  void compute_per_surf();
  virtual void clear();
  virtual void surf_tally(int, Particle::OnePart *, 
                          Particle::OnePart *, Particle::OnePart *);
  virtual int tallyinfo(surfint *&);
  virtual void tallysum(int);
  void reallocate();
  bigint memory_usage();

 protected:
  int groupbit,imix,nvalue,ngroup,ntotal;
  int *which;
  bigint last_tallysum;      // last timestep tallysum was called
  double nfactor_inverse;

  int ntally;              // # of surfs I have tallied for
  int maxtally;            // # of tallies currently allocated
  double **array;          // tally values, maxtally in length
  surfint *tally2surf;     // tally2surf[I] = surf ID of Ith tally

  // hash for surf IDs

#ifdef SPARTA_MAP
  typedef std::map<surfint,int> MyHash;
#elif defined SPARTA_UNORDERED_MAP
  typedef std::unordered_map<surfint,int> MyHash;
#else
  typedef std::tr1::unordered_map<surfint,int> MyHash;
#endif

  MyHash *hash;

  int dim;                 // local copies
  Surf::Line *lines;
  Surf::Tri *tris;

  int weightflag;          // 1 if cell weighting is enabled
  double weight;           // particle weight, based on initial cell
  double *normflux;        // normalization factor for each surf element

  virtual void init_normflux();
  void grow_tally();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute surf mixture ID does not exist

Self-explanatory.

E: Number of groups in compute surf mixture has changed

This mixture property cannot be changed after this compute command is
issued.

*/
