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

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(radeq,SurfCollideRadeq)

#else

#ifndef SPARTA_SURF_COLLIDE_RADEQ_H
#define SPARTA_SURF_COLLIDE_RADEQ_H

#include "surf_collide.h"
#include "fix.h"

namespace SPARTA_NS {

class SurfCollideRadeq : public SurfCollide {
 public:
  SurfCollideRadeq(class SPARTA *, int, char **);
  SurfCollideRadeq(class SPARTA *sparta) : SurfCollide(sparta) {}
  ~SurfCollideRadeq();
  void init();
  Particle::OnePart *collide(Particle::OnePart *&, double *, double &,
                             int, int &, int);

 protected:
  int m;
  bigint nvalid;
  double twall,prefactor,emi;
  char *id_qw;
  int qwindex;
  class Fix *fqw;
  double *qw;

  int distributed,implicit;  // Surf settings
  int firstflag;

  int groupbit;              // mask for surface group
  int dimension;
  int nown;                  // # of surf elements owned by this proc
  bigint nsurf;              // total # of surf elements, lines or tris
  int nchoose;               // # of surf elements output by this proc
  int *cglobal;              // indices of global elements for nchoose

  double vstream[3];
  class RanPark *random;     // RNG for particle reflection

  double radeq(Particle::OnePart *, double *, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot use surf_collide radeq with implicit surfs

Self-explanatory.

E: Surf_collide radeq twall <= 0.0

The starting temperature needs to be strictly positive.

E: Invalid qw in Surf_collide radeq command

Self-explanatory.

E: Could not find Surf_collide radeq fix ID

Self-explanatory.

E: Surf_collide radeq fix does not compute per-surf info

Self-explanatory.

E: Surf_collide radeq fix does not compute per-surf vector

Self-explanatory.

E: Surf_collide radeq fix does not compute per-surf array

Self-explanatory.

E: Surf_collide radeq fix array is accessed out-of-range

Self-explanatory.

E: Emissivity of surface has to be between 0 and 1

Self-explanatory.

E: Surf_collide radeq group ID does not exist

Self-explanatory.
*/
