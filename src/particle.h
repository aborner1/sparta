/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#ifndef DSMC_PARTICLE_H
#define DSMC_PARTICLE_H

#include "pointers.h"

namespace DSMC_NS {

class Particle : protected Pointers {
 public:
  bigint nmove;             // dummy statistic for now
  bigint ncellcross;        // dummy statistic for now

  int nmigrate;             // # of particles to migrate to new procs
  int *mlist;               // indices of particles to migrate

  struct OnePart {
    int id,type;            // particle ID, type
    int icell;              // grid cell the particle is in (0 to N-1)
    double x[3];            // coords of particle
    double v[3];            // velocity of particle
  };

  bigint nglobal;           // global # of particles
  int nlocal;               // # of particles I own
  OnePart *particles;       // list of particles I own

  Particle(class DSMC *);
  ~Particle();
  void init() {}
  void create(int, char **);
  void move();
  void compress();
  void grow(int);
  bigint memory_usage();

 private:
  int seed;

  int maxlocal;             // max # particles list can hold
  int maxmigrate;           // max # migrate list can hold

  void create_all(bigint);
  void create_local(bigint);
};

}

#endif
