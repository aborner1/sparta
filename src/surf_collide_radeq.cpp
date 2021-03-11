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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "surf_collide_radeq.h"
#include "grid.h"
#include "surf.h"
#include "surf_react.h"
#include "input.h"
#include "variable.h"
#include "particle.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "mpi.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

SurfCollideRadeq::SurfCollideRadeq(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal surf_collide radeq command");

  if (surf->implicit)
    error->all(FLERR,"Cannot use surf_collide radeq with implicit surfs");

  id_qw = NULL;

  twall = input->numeric(FLERR,arg[2]);
  if (twall <= 0.0) error->all(FLERR,"Surf_collide radeq twall <= 0.0");

  if (strncmp(arg[3],"f_",2) == 0) {
    int n = strlen(arg[3]);
    id_qw = new char[n];
    strcpy(id_qw,&arg[3][2]);

    char *ptr = strchr(id_qw,'[');
    if (ptr) {
      if (id_qw[strlen(id_qw)-1] != ']')
        error->all(FLERR,"Invalid qw in Surf_collide radeq command");
      qwindex = atoi(ptr+1);
      *ptr = '\0';
    } else qwindex = 0;

  m = modify->find_fix(id_qw);
  if (m < 0) error->all(FLERR,"Could not find Surf_collide radeq fix ID");
  if (modify->fix[m]->per_surf_flag == 0)
    error->all(FLERR,"Surf_collide radeq fix does not "
               "compute per-surf info");
  if (qwindex == 0 && modify->fix[m]->size_per_surf_cols > 0)
    error->all(FLERR,"Surf_collide radeq fix does not "
               "compute per-surf vector");
  if (qwindex > 0 && modify->fix[m]->size_per_surf_cols == 0)
    error->all(FLERR,"Surf_collide radeq fix does not "
               "compute per-surf array");
  if (qwindex > 0 && qwindex > modify->fix[m]->size_per_surf_cols)
    error->all(FLERR,"Surf_collide radeq fix array is "
               "accessed out-of-range");
  }

  emi = input->numeric(FLERR,arg[4]);
  if (emi <= 0.0 || emi > 1.0)
    error->all(FLERR,"Emissivity of surface has to be between 0 and 1");

  prefactor = 1.0 / (emi *  5.670374419e-8);

  // initialize data structures
  // trigger setup of list of owned surf elements belonging to surf group

  firstflag = 1;
  qw = NULL;
  cglobal = NULL;

  // initialize RNG

  random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideRadeq::~SurfCollideRadeq()
{
  if (copy) return;

  delete [] id_qw;
  memory->destroy(qw);
  memory->destroy(cglobal);

  delete random;
}

/* ---------------------------------------------------------------------- */

void SurfCollideRadeq::init()
{
  SurfCollide::init();

  if (!firstflag) return;
  firstflag = 0;

  int ifix = modify->find_fix(id_qw);
  if (ifix < 0)
    error->all(FLERR,"Could not find surf_collide radeq fix ID");
  fqw = modify->fix[ifix];

  nvalid = fqw->per_surf_freq + 1;

  // one-time setup of lists of owned elements contributing to radeq
  // NOTE: will need to recalculate, if allow addition of surf elements
  // nown = # of surf elements I own
  // nchoose = # of nown surf elements in surface group
  // cglobal[] = global indices for nchoose elements
  //             used to access lines/tris in Surf
  // clocal[] = local indices for nchoose elements
  //            used to access nown data from per-surf computes,fixes,variables

  dimension = domain->dimension;
  distributed = surf->distributed;
  implicit = surf->implicit;

  int igroup = 0;
  if (igroup < 0) error->all(FLERR,"Surf_collide radeq group ID does not exist");
  groupbit = surf->bitmask[igroup];

  Surf::Line *lines;
  Surf::Tri *tris;

  if (distributed && !implicit) lines = surf->mylines;
  else lines = surf->lines;
  if (distributed && !implicit) tris = surf->mytris;
  else tris = surf->tris;

  nown = surf->nown;
  nsurf = surf->nsurf;
  int m;
  int me = comm->me;
  int nprocs = comm->nprocs;

  nchoose = 0;
  for (int i = 0; i < nown; i++) {
    if (dimension == 2) {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (lines[m].mask & groupbit) nchoose++;
    } else {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (tris[m].mask & groupbit) nchoose++;
    }
  }

  memory->create(cglobal,nchoose,"surf_collide radeq:cglobal");

  nchoose = 0;
  for (int i = 0; i < nown; i++) {
    if (dimension == 2) {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (lines[m].mask & groupbit) {
        cglobal[nchoose++] = m;
      }
    } else {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (tris[m].mask & groupbit) {
        cglobal[nchoose++] = m;
      }
    }
  }

  memory->create(qw,nsurf,"surf_collide radeq:init");
  for (int i = 0; i < nsurf; i++) qw[i] = 0.0;
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   ip = set to NULL if destroyed by chemsitry
   return jp = new particle if created by chemistry
   return reaction = index of reaction (1 to N) that took place, 0 = no reaction
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideRadeq::

collide(Particle::OnePart *&ip, double *norm, double &, int isr, int &reaction, int isurf)

{
  nsingle++;

  // if surface chemistry defined, attempt reaction
  // reaction > 0 if reaction took place

  Particle::OnePart iorig;
  Particle::OnePart *jp = NULL;
  reaction = 0;

  if (isr >= 0) {
    if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    reaction = surf->sr[isr]->react(ip,norm,jp);
    if (reaction) surf->nreact_one++;
  }

  // diffuse reflection for each particle
  // resets v, roteng, vibeng
  // if new particle J created, also need to trigger any fixes

  if (ip) {
    twall = radeq(ip,norm,isurf);
    if (modify->n_add_particle) {
    int i = ip - particle->particles;
    modify->add_particle(i,twall,twall,twall,vstream);
    }
  }
  if (jp) {
    twall = radeq(jp,norm,isurf);
    if (modify->n_add_particle) {
      int j = jp - particle->particles;
      modify->add_particle(j,twall,twall,twall,vstream);
    }
  }

  // call any fixes with a surf_react() method
  // they may reset j to -1, e.g. fix ambipolar
  //   in which case newly created j is deleted

  if (reaction && modify->n_surf_react) {
    int i = -1;
    if (ip) i = ip - particle->particles;
    int j = -1;
    if (jp) j = jp - particle->particles;
    modify->surf_react(&iorig,i,j);
    if (jp && j < 0) {
      jp = NULL;
      particle->nlocal--;
    }
  }

  return jp;
}

/* ----------------------------------------------------------------------
   diffusive particle collision with surface
   p = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

double SurfCollideRadeq::radeq(Particle::OnePart *p, double *norm, int jsurf)
{
  // diffuse reflection
  // vrm = most probable speed of species, eqns (4.1) and (4.7)
  // vperp = velocity component perpendicular to surface along norm, eqn (12.3)
  // vtan12 = 2 velocity components tangential to surface
  // tangent1 = component of particle v tangential to surface,
  //   check if tangent1 = 0 (normal collision), set randomly
  // tangent2 = norm x tangent1 = orthogonal tangential direction
  // tangent12 are both unit vectors

    double tangent1[3],tangent2[3];
    Particle::Species *species = particle->species;
    int ispecies = p->ispecies;
    double twall_new = twall;

    if (update->ntimestep >= nvalid) {
      nvalid = (update->ntimestep/fqw->per_surf_freq)*fqw->per_surf_freq +
        fqw->per_surf_freq;
      nvalid -= (fqw->per_surf_freq-1)*fqw->nevery;
      if (nvalid <= update->ntimestep) nvalid += fqw->per_surf_freq;

      if (qwindex == 0) {
        double *vector = fqw->vector_surf;
        for (int i = 0; i < nchoose; i++)
          qw[cglobal[i]] = vector[i];
      }
      else {
        double **array = fqw->array_surf;
        int index = qwindex-1;
        for (int i = 0; i < nchoose; i++)
          qw[cglobal[i]] = array[i][index];
      }
    }
     
     if ((qw[jsurf] > 300.0) && (qw[jsurf] < 1.e7)) {
        twall_new = pow((prefactor * qw[jsurf]),0.25);
     }

    double vrm = sqrt(2.0*update->boltz * twall_new / species[ispecies].mass);
    double vperp = vrm * sqrt(-log(random->uniform()));

    double theta = MY_2PI * random->uniform();
    double vtangent = vrm * sqrt(-log(random->uniform()));
    double vtan1 = vtangent * sin(theta);
    double vtan2 = vtangent * cos(theta);

    double *v = p->v;
    double dot = MathExtra::dot3(v,norm);

    tangent1[0] = v[0] - dot*norm[0];
    tangent1[1] = v[1] - dot*norm[1];
    tangent1[2] = v[2] - dot*norm[2];

    if (MathExtra::lensq3(tangent1) == 0.0) {
      tangent2[0] = random->uniform();
      tangent2[1] = random->uniform();
      tangent2[2] = random->uniform();
      MathExtra::cross3(norm,tangent2,tangent1);
    }

    MathExtra::norm3(tangent1);
    MathExtra::cross3(norm,tangent1,tangent2);

    v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
    v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
    v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

    // initialize rot/vib energy

    p->erot = particle->erot(ispecies,twall_new,random);
    p->evib = particle->evib(ispecies,twall_new,random);

    return twall_new;
}
