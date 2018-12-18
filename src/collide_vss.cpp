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
#include "string.h"
#include "stdlib.h"
#include "collide_vss.h"
#include "grid.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "collide.h"
#include "react.h"
#include "comm.h"
#include "fix_vibmode.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NONE,DISCRETE,SMOOTH};            // several files
enum{CONSTANT,VARIABLE};

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

CollideVSS::CollideVSS(SPARTA *sparta, int narg, char **arg) :
  Collide(sparta, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal collide command");

  // optional args

  relaxflag = CONSTANT;
  relaxtypeflag = SERIAL;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"relax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal collide command");
      if (strcmp(arg[iarg+1],"constant") == 0) relaxflag = CONSTANT;
      else if (strcmp(arg[iarg+1],"variable") == 0) relaxflag = VARIABLE;
      else error->all(FLERR,"Illegal collide command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"relaxtype") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal collide command");
      if (strcmp(arg[iarg+1],"serial") == 0) relaxtypeflag = SERIAL;
      else if (strcmp(arg[iarg+1],"prohibdouble") == 0) relaxtypeflag = PROHIBDOUBLE;
      else error->all(FLERR,"Illegal collide command");
      iarg += 2;	  
    } else error->all(FLERR,"Illegal collide command");
  }

  // proc 0 reads file to extract params for current species
  // broadcasts params to all procs

  nparams = particle->nspecies;
  if (nparams == 0) 
    error->all(FLERR,"Cannot use collide command with no species defined");

  params = new Params[nparams];
  if (comm->me == 0) read_param_file(arg[2]);
  MPI_Bcast(params,nparams*sizeof(Params),MPI_BYTE,0,world);

  // check that params were read for all species

  for (int i = 0; i < nparams; i++)
    if (params[i].diam < 0.0) {
      char str[128];
      sprintf(str,"Species %s did not appear in VSS parameter file",
	      particle->species[i].id);
      error->one(FLERR,str);
    }

  // allocate per-species prefactor array

  memory->create(prefactor,nparams,nparams,"collide:prefactor");
}

/* ---------------------------------------------------------------------- */

CollideVSS::~CollideVSS()
{
  if (copymode) return;

  delete [] params;
  memory->destroy(prefactor);
}

/* ---------------------------------------------------------------------- */

void CollideVSS::init()
{
  // initially read-in per-species params must match current species list

  if (nparams != particle->nspecies)
    error->all(FLERR,"VSS parameters do not match current species");

  Collide::init();
}

/* ----------------------------------------------------------------------
   estimate a good value for vremax for a group pair in any grid cell
   called by Collide parent in init()
------------------------------------------------------------------------- */

double CollideVSS::vremax_init(int igroup, int jgroup)
{
  // parent has set mixture ptr

  Particle::Species *species = particle->species;
  double *vscale = mixture->vscale;
  int *mix2group = mixture->mix2group;
  int nspecies = particle->nspecies;

  double vrmgroup = 0.0;

  for (int isp = 0; isp < nspecies; isp++) {
    if (mix2group[isp] != igroup) continue;
    for (int jsp = 0; jsp < nspecies; jsp++) {
      if (mix2group[jsp] != jgroup) continue;

      double diam = 0.5 * (params[isp].diam + params[jsp].diam);
      double omega = 0.5 * (params[isp].omega + params[jsp].omega);
      double tref = 0.5 * (params[isp].tref + params[jsp].tref);
      double mr = species[isp].mass * species[jsp].mass /
	(species[isp].mass + species[jsp].mass);
      double cxs = diam*diam*MY_PI;
      prefactor[isp][jsp] = cxs *
	pow(2.0*update->boltz*tref/mr,omega-0.5)/tgamma(2.5-omega);
      double beta = MAX(vscale[isp],vscale[jsp]);
      double vrm = 2.0 * cxs * beta;
      vrmgroup = MAX(vrmgroup,vrm);
    }
  }

  return vrmgroup;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int icell, int np, double volume)
{
 double fnum = update->fnum;
 double dt = update->dt;

 double nattempt;

 if (remainflag) {
   nattempt = 0.5 * np * (np-1) *
     vremax[icell][0][0] * dt * fnum / volume + remain[icell][0][0];
   remain[icell][0][0] = nattempt - static_cast<int> (nattempt);
 } else 
   nattempt = 0.5 * np * (np-1) *
     vremax[icell][0][0] * dt * fnum / volume + random->uniform();

 // DEBUG
 //nattempt = 10;

  return nattempt;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int icell, int igroup, int jgroup, 
				     double volume)
{
 double fnum = update->fnum;
 double dt = update->dt;

 double nattempt;

 // return 2x the value for igroup != jgroup, since no J,I pairing

 double npairs;
 if (igroup == jgroup) npairs = 0.5 * ngroup[igroup] * (ngroup[igroup]-1);
 else npairs = ngroup[igroup] * (ngroup[jgroup]);
 //else npairs = 0.5 * ngroup[igroup] * (ngroup[jgroup]);

 nattempt = npairs * vremax[icell][igroup][jgroup] * dt * fnum / volume;

 if (remainflag) {
   nattempt += remain[icell][igroup][jgroup];
   remain[icell][igroup][jgroup] = nattempt - static_cast<int> (nattempt);
 } else nattempt += random->uniform();

 return nattempt;
}

/* ----------------------------------------------------------------------
   determine if collision actually occurs
   1 = yes, 0 = no
   update vremax either way
------------------------------------------------------------------------- */

int CollideVSS::test_collision(int icell, int igroup, int jgroup,
			       Particle::OnePart *ip, Particle::OnePart *jp)
{
  double *vi = ip->v;
  double *vj = jp->v;
  int ispecies = ip->ispecies;
  int jspecies = jp->ispecies;
  double du  = vi[0] - vj[0];
  double dv  = vi[1] - vj[1];
  double dw  = vi[2] - vj[2];
  double vr2 = du*du + dv*dv + dw*dw;
  double omega1 = params[ispecies].omega;
  double omega2 = params[jspecies].omega;
  double omega = 0.5 * (omega1+omega2);
  double vro  = pow(vr2,1.0-omega);

  // although the vremax is calcualted for the group,
  // the individual collisions calculated species dependent vre

  double vre = vro*prefactor[ispecies][jspecies];
  vremax[icell][igroup][jgroup] = MAX(vre,vremax[icell][igroup][jgroup]);
  if (vre/vremax[icell][igroup][jgroup] < random->uniform()) return 0;
  precoln.vr2 = vr2;
  return 1;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::setup_collision(Particle::OnePart *ip, Particle::OnePart *jp)
{
  Particle::Species *species = particle->species;

  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  precoln.vr = sqrt(precoln.vr2);

  precoln.ave_rotdof = 0.5 * (species[isp].rotdof + species[jsp].rotdof);
  precoln.ave_vibdof = 0.5 * (species[isp].vibdof + species[jsp].vibdof);
  precoln.ave_dof = (precoln.ave_rotdof  + precoln.ave_vibdof)/2.;

  double imass = precoln.imass = species[isp].mass;
  double jmass = precoln.jmass = species[jsp].mass;
  precoln.mr = imass * jmass / (imass+jmass);

  precoln.etrans = 0.5 * precoln.mr * precoln.vr2;
  precoln.erot = ip->erot + jp->erot;
  precoln.evib = ip->evib + jp->evib;

  precoln.eint   = precoln.erot + precoln.evib;
  precoln.etotal = precoln.etrans + precoln.eint;

  // COM velocity calculated using reactant masses

  double divisor = 1.0 / (imass+jmass);
  double *vi = ip->v;
  double *vj = jp->v;
  precoln.ucmf = ((imass*vi[0])+(jmass*vj[0])) * divisor;
  precoln.vcmf = ((imass*vi[1])+(jmass*vj[1])) * divisor;
  precoln.wcmf = ((imass*vi[2])+(jmass*vj[2])) * divisor;

  postcoln.etrans = precoln.etrans;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  postcoln.eint = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */

int CollideVSS::perform_collision(Particle::OnePart *&ip, 
                                  Particle::OnePart *&jp, 
                                  Particle::OnePart *&kp)
{
  int reactflag,kspecies;
  double x[3],v[3];
  Particle::OnePart *p3;

  // if gas-phase chemistry defined, attempt and perform reaction
  // if a 3rd particle is created, its kspecies >= 0 is returned
  // if 2nd particle is removed, its jspecies is set to -1

  if (react) 
    reactflag = react->attempt(ip,jp,
                               precoln.etrans,precoln.erot,
                               precoln.evib,postcoln.etotal,kspecies);
  else reactflag = 0;
  
  // repartition energy and perform velocity scattering for I,J,K particles
  // reaction may have changed species of I,J particles
  // J,K particles may have been removed or created by reaction

  kp = NULL;

  if (reactflag) {

    // add 3rd K particle if reaction created it
    // index of new K particle = nlocal-1
    // if add_particle() performs a realloc:
    //   make copy of x,v, then repoint ip,jp to new particles data struct

    if (kspecies >= 0) {
      int id = MAXSMALLINT*random->uniform();

      Particle::OnePart *particles = particle->particles;
      memcpy(x,ip->x,3*sizeof(double));
      memcpy(v,ip->v,3*sizeof(double));
      int reallocflag = 
        particle->add_particle(id,kspecies,ip->icell,x,v,0.0,0.0);
      if (reallocflag) {
        ip = particle->particles + (ip - particles);
        jp = particle->particles + (jp - particles);
      }

      kp = &particle->particles[particle->nlocal-1];
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_ThreeBodyScattering(ip,jp,kp);

    // remove 2nd J particle if recombination reaction removed it
    // p3 is 3rd particle participating in energy exchange

    } else if (jp->ispecies < 0) {
      double *vi = ip->v;
      double *vj = jp->v;

      double divisor = 1.0 / (precoln.imass + precoln.jmass);
      double ucmf = ((precoln.imass*vi[0]) + (precoln.jmass*vj[0])) * divisor;
      double vcmf = ((precoln.imass*vi[1]) + (precoln.jmass*vj[1])) * divisor;
      double wcmf = ((precoln.imass*vi[2]) + (precoln.jmass*vj[2])) * divisor;

      vi[0] = ucmf;
      vi[1] = vcmf;
      vi[2] = wcmf;

      jp = NULL;
      p3 = react->recomb_part3;
      setup_collision(ip,p3);
      if (precoln.ave_dof > 0.0) EEXCHANGE_ReactingEDisposal(ip,p3,jp);
      SCATTER_TwoBodyScattering(ip,p3);

    } else {
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_TwoBodyScattering(ip,jp);
    }

  } else { 
    if (precoln.ave_dof > 0.0) {
      if (relaxtypeflag == PROHIBDOUBLE} EEXCHANGE_NonReactingEDisposal_ProhibDouble(ip,jp);
      else EEXCHANGE_NonReactingEDisposal_Serial(ip,jp);
    }
    SCATTER_TwoBodyScattering(ip,jp);
  }

  return reactflag;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::SCATTER_TwoBodyScattering(Particle::OnePart *ip, 
					   Particle::OnePart *jp)
{
  double ua,vb,wc;
  double vrc[3];

  Particle::Species *species = particle->species;
  double *vi = ip->v;
  double *vj = jp->v;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  double mass_i = species[isp].mass;
  double mass_j = species[jsp].mass;

  double alpha_r = 2.0 / (params[isp].alpha + params[jsp].alpha);
  double mr = species[isp].mass * species[jsp].mass /
    (species[isp].mass + species[jsp].mass);

  double eps = random->uniform() * 2*MY_PI;
  if (fabs(alpha_r - 1.0) < 0.001) { 
    double vr = sqrt(2.0 * postcoln.etrans / mr);
    double cosX = 2.0*random->uniform() - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps); 
  } else {
    double scale = sqrt((2.0 * postcoln.etrans) / (mr * precoln.vr2));
    double cosX = 2.0*pow(random->uniform(),alpha_r) - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.0e-6) { 
      ua = scale * ( cosX*vrc[0] + sinX*d*sin(eps) );
      vb = scale * ( cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) - 
                                         vrc[0]*vrc[1]*sin(eps))/d );
      wc = scale * ( cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) + 
                                         vrc[0]*vrc[2]*sin(eps))/d );
    } else {
      ua = scale * ( cosX*vrc[0] ); 
      vb = scale * ( sinX*vrc[1]*cos(eps) ); 
      wc = scale * ( sinX*vrc[2]*sin(eps) );
    }
  }
  
  // new velocities for the products

  double divisor = 1.0 / (mass_i + mass_j);
  // Changed signed of second term to match DS1V
  vi[0] = precoln.ucmf + (mass_j*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_j*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_j*divisor)*wc;
  vj[0] = precoln.ucmf - (mass_i*divisor)*ua;
  vj[1] = precoln.vcmf - (mass_i*divisor)*vb;
  vj[2] = precoln.wcmf - (mass_i*divisor)*wc;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_NonReactingEDisposal_Serial(Particle::OnePart *ip, 
						Particle::OnePart *jp)
{

  double State_prob,Fraction_Rot,Fraction_Vib,E_Dispose;
  int i,rotdof,vibdof,max_level,ivib,irot;

  Particle::OnePart *p;
  Particle::Species *species = particle->species;

  double AdjustFactor = 0.99999999;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  double pevib = 0.0;   

  // handle each kind of energy disposal for non-reacting reactants

  if (precoln.ave_dof == 0) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;

  } else {
    E_Dispose = precoln.etrans;

    double aveomega = 0.5*(params[ip->ispecies].omega + 
                           params[jp->ispecies].omega); 

    for (i = 0; i < 2; i++) {
      if (i == 0) p = ip; 
      else p = jp;

      int sp = p->ispecies;
      rotdof = species[sp].rotdof;
      double rotn_phi = 1.0/species[sp].rotrel; 

      if (rotdof) {
        if (relaxflag == VARIABLE) rotn_phi = rotrel_serial(sp,E_Dispose);
        if (rotn_phi >= random->uniform()) {
          if (rotstyle == NONE) {
            p->erot = 0.0; 
          } else if (rotstyle != NONE && rotdof == 2) {
            E_Dispose += p->erot;
            Fraction_Rot = 
              1- pow(random->uniform(),(1/(2.5-aveomega)));
            p->erot = Fraction_Rot * E_Dispose;
            E_Dispose -= p->erot;
          } else {
            E_Dispose += p->erot;
            p->erot = E_Dispose * 
              sample_bl(random,0.5*species[sp].rotdof-1.0,
                        1.5-params[sp].omega);
            E_Dispose -= p->erot;
          }
        }
      }
      postcoln.erot += p->erot;
 
      vibdof = species[sp].vibdof;
      double vibn_phi = 1.0/species[sp].vibrel[0]; 

      if (vibdof) {
        if (relaxflag == VARIABLE) vibn_phi = vibrel_serial(sp,E_Dispose+p->evib);
        if (vibn_phi >= random->uniform()) {
          if (vibstyle == NONE) {
            p->evib = 0.0; 

          } else if (vibdof == 2) {
            if (vibstyle == SMOOTH) {
              E_Dispose += p->evib;
              Fraction_Vib = 
                1.0 - pow(random->uniform(),(1.0/(2.5-aveomega)));
              p->evib= Fraction_Vib * E_Dispose;
              E_Dispose -= p->evib;

            } else if (vibstyle == DISCRETE) {
              E_Dispose += p->evib;
              max_level = static_cast<int>
                (E_Dispose / (update->boltz * species[sp].vibtemp[0]));
              do {
                ivib = static_cast<int> 
                  (random->uniform()*(max_level+AdjustFactor));
                p->evib = ivib * update->boltz * species[sp].vibtemp[0];
                State_prob = pow((1.0 - p->evib / E_Dispose),
                                 (1.5 - params[sp].omega));
              } while (State_prob < random->uniform());
              E_Dispose -= p->evib;
            }
            
          } else if (vibdof > 2) {
            if (vibstyle == SMOOTH) {
              E_Dispose += p->evib;
              p->evib = E_Dispose * 
                sample_bl(random,0.5*species[sp].vibdof-1.0,
                          1.5-aveomega);
              E_Dispose -= p->evib;
              
            } else if (vibstyle == DISCRETE) {
              p->evib = 0.0;
              
              int nmode = particle->species[sp].nvibmode;
              int **vibmode = 
                particle->eiarray[particle->ewhich[index_vibmode]];
              int pindex = p - particle->particles; 

              for (int imode = 0; imode < nmode; imode++) {
                ivib = vibmode[pindex][imode];
                E_Dispose += ivib * update->boltz * 
                  particle->species[sp].vibtemp[imode];
                max_level = static_cast<int>
                  (E_Dispose / (update->boltz * species[sp].vibtemp[imode]));

                do {
                  ivib = static_cast<int> 
                    (random->uniform()*(max_level+AdjustFactor));
                  pevib = ivib * update->boltz * species[sp].vibtemp[imode];
                  State_prob = pow((1.0 - pevib / E_Dispose),
                                   (1.5 - aveomega));
                } while (State_prob < random->uniform());
                
                vibmode[pindex][imode] = ivib;
                p->evib += pevib;
                E_Dispose -= pevib;
              }
            }
          } // end of vibstyle/vibdof if
        } 
        postcoln.evib += p->evib;
      } // end of vibdof if
    }
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_NonReactingEDisposal_ProhibDouble(Particle::OnePart *ip, 
						Particle::OnePart *jp)
{

  double State_prob,Fraction_Vib,E_Dispose;
  double phi,factor,transdof,pevib,Tt,vibdof_dis,A;
  int i,sp,rotdof,vibdof,max_level,ivib,imode,relaxflag1;

  Particle::OnePart *p1,*p2,*p;
  Particle::Species *species = particle->species;

  double AdjustFactor = 0.99999999;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  relaxflag1 = 0;
  phi = 0.0;
  factor= 1.0;

  // handle each kind of energy disposal for non-reacting reactants

  if (precoln.ave_dof == 0) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;

  } else {
    E_Dispose = precoln.etrans;

    do {
      if (0.5 < random->uniform()) {
          p1 = ip;
          p2 = jp;
      }   else {
          p1 = jp;
          p2 = ip;
      }

      for (i = 0; i < 2; i++) {
        if (i == 0) p = p1;
        else if (i == 1) p = p2;

          sp = p->ispecies;
          vibdof = species[sp].vibdof;
          transdof = 5.0-2.0*params[sp].omega;

          Tt = (E_Dispose+p->evib) / (update->boltz * (transdof+vibdof));

          if (vibdof) {
              if (vibstyle == NONE) {
                p->evib = 0.0;

              } else if (vibstyle == SMOOTH) {

                  factor *= 1/(1-phi);
                  if (relaxflag == VARIABLE) phi = factor*(1.0 + vibdof/transdof)/vibrel_prohibdouble(sp,E_Dispose+p->evib,vibdof+transdof);
                  else phi = factor*(1.0 + vibdof/transdof)/species[sp].vibrel[0];

                  if (phi >= random->uniform()) {
                      E_Dispose += p->evib;
                      if (vibdof == 2) {
                          Fraction_Vib = 1.0 - pow(random->uniform(),(1.0/(2.5-params[sp].omega)));
                          p->evib= Fraction_Vib * E_Dispose;
                      }
                      else if (vibdof > 2) p->evib = E_Dispose * sample_bl(random,0.5*vibdof-1.0,1.5-params[sp].omega);
                      E_Dispose -= p->evib;
                      postcoln.evib += p->evib;
                      relaxflag1 = 1;
                      break;
                  }

              } else if (vibstyle == DISCRETE) {
                  if (vibdof == 2) {

                      vibdof_dis = 2.0*(species[sp].vibtemp[0]/Tt) / (exp(species[sp].vibtemp[0]/Tt)-1);
                      A = pow(vibdof_dis,2)*exp(species[sp].vibtemp[0]/Tt)/2.0;

                      factor *= 1/(1-phi);
                      if (relaxflag == VARIABLE) phi = factor*(1.0 + A/transdof)/vibrel_prohibdouble(sp,E_Dispose+p->evib,vibdof_dis+transdof);
                      else phi = factor*(1.0 + A/transdof)/species[sp].vibrel[0];

                      if (phi >= random->uniform()) {
                          E_Dispose += p->evib;
                          max_level = static_cast<int>
                            (E_Dispose / (update->boltz * species[sp].vibtemp[0]));
                          do {
                            ivib = static_cast<int>
                              (random->uniform()*(max_level+AdjustFactor));
                            p->evib = ivib * update->boltz * species[sp].vibtemp[0];
                            State_prob = pow((1.0 - p->evib / E_Dispose),
                                             (1.5 - params[sp].omega));
                          } while (State_prob < random->uniform());
                          E_Dispose -= p->evib;
                          postcoln.evib += p->evib;
                          relaxflag1 = 1;
                          break;
                      }

                  } else if (vibdof > 2) {

                      int nmode = particle->species[sp].nvibmode;
                      int **vibmode =
                        particle->eiarray[particle->ewhich[index_vibmode]];
                      int pindex = p - particle->particles;

                      imode = 0;
                      while (imode < nmode) {

                        vibdof_dis = 2.0*(species[sp].vibtemp[imode]/Tt) / (exp(species[sp].vibtemp[imode]/Tt)-1);
                        A = pow(vibdof_dis,2)*exp(species[sp].vibtemp[imode]/Tt)/2.0;

                        factor *= 1/(1-phi);
                        if (relaxflag == VARIABLE) phi = factor*(1.0 + A/transdof)/vibrel_prohibdouble(sp,E_Dispose+p->evib,vibdof_dis+transdof);
                        else phi = factor*(1.0 + A/transdof)/species[sp].vibrel[imode];

                        if (phi >= random->uniform()) {
                            ivib = vibmode[pindex][imode];
                            E_Dispose += ivib * update->boltz *
                                    particle->species[sp].vibtemp[imode];
                            p->evib -= ivib * update->boltz *
                                    particle->species[sp].vibtemp[imode];
                            max_level = static_cast<int>
                              (E_Dispose / (update->boltz * species[sp].vibtemp[imode]));

                            do {
                              ivib = static_cast<int>
                                (random->uniform()*(max_level+AdjustFactor));
                              pevib = ivib * update->boltz * species[sp].vibtemp[imode];
                              State_prob = pow((1.0 - pevib / E_Dispose),
                                               (1.5 - params[sp].omega));
                            } while (State_prob < random->uniform());

                            vibmode[pindex][imode] = ivib;
                            p->evib += pevib;
                            E_Dispose -= pevib;
                            postcoln.evib += pevib;
                            relaxflag1 = 1;
                            break;
                        }
                        imode++;
                      }
                      if (relaxflag1) break;
                    }
                } // end of vibstyle if
            } // end of vibdof if
          if (relaxflag1) break;
      }

      if (relaxflag1) break;

      for (i = 0; i < 2; i++) {
        if (i == 0) p = p1;
        else if (i == 1) p = p2;

          sp = p->ispecies;
          rotdof = species[sp].rotdof;
          transdof = 5.0-2.0*params[sp].omega;

          if (rotdof) {
            if (rotstyle == NONE) {
                p->erot = 0.0;
            } else {
                factor *= 1/(1-phi);
                if (relaxflag == VARIABLE) phi = factor*(1.0 + rotdof/transdof)/rotrel_prohibdouble(sp,E_Dispose+p->erot);
                else phi = factor*(1.0 + rotdof/transdof)/species[sp].rotrel;

                if (phi >= random->uniform()) {
                    E_Dispose += p->erot;
                    if (rotdof == 2) {
                        Fraction_Vib = 1.0 - pow(random->uniform(),(1.0/(2.5-params[sp].omega)));
                        p->erot= Fraction_Vib * E_Dispose;
                    }
                    else if (rotdof > 2) p->erot = E_Dispose * sample_bl(random,0.5*rotdof-1.0,1.5-params[sp].omega);
                    E_Dispose -= p->erot;
                    postcoln.erot += p->erot;
                    break;
                }
            }
          }
      }
    } while (0 > 1);
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::SCATTER_ThreeBodyScattering(Particle::OnePart *ip, 
			  		     Particle::OnePart *jp,
			  		     Particle::OnePart *kp)
{
  double vrc[3],ua,vb,wc;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  int ksp = kp->ispecies;
  double mass_i = species[isp].mass;
  double mass_j = species[jsp].mass;
  double mass_k = species[ksp].mass;
  double mass_ij = mass_i + mass_j;
  double *vi = ip->v;
  double *vj = jp->v;
  double *vk = kp->v;

  double alpha_r = 2.0 / (params[isp].alpha + params[jsp].alpha);
  double mr = mass_ij * mass_k / (mass_ij + mass_k);
  postcoln.eint = ip->erot + jp->erot + ip->evib + jp->evib 
                + kp->erot + kp->evib;

  double cosX = 2.0*pow(random->uniform(), alpha_r) - 1.0;
  double sinX = sqrt(1.0 - cosX*cosX);
  double eps = random->uniform() * 2*MY_PI;

  if (fabs(alpha_r - 1.0) < 0.001) { 
    double vr = sqrt(2*postcoln.etrans/mr);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    double scale = sqrt((2.0*postcoln.etrans) / (mr*precoln.vr2));
    // Should this be vij-vk?? (Arnaud 12/17/18)
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.E-6 ) { 
      ua = scale * (cosX*vrc[0] + sinX*d*sin(eps));
      vb = scale * (cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) - 
                                        vrc[0]*vrc[1]*sin(eps))/d);
      wc = scale * (cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) + 
                                        vrc[0]*vrc[2]*sin(eps))/d);
    } else {
      ua = scale * cosX*vrc[0]; 
      vb = scale * sinX*vrc[1]*cos(eps); 
      wc = scale * sinX*vrc[2]*sin(eps);
    }
  }

  // new velocities for the products

  double divisor = 1.0 / (mass_ij + mass_k);
  // Changed signed of second term to match DS1V
  vi[0] = precoln.ucmf + (mass_ij*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_ij*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_ij*divisor)*wc;
  vk[0] = precoln.ucmf - (mass_k*divisor)*ua;
  vk[1] = precoln.vcmf - (mass_k*divisor)*vb;
  vk[2] = precoln.wcmf - (mass_k*divisor)*wc;
  vj[0] = vi[0];
  vj[1] = vi[1];
  vj[2] = vi[2];
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_ReactingEDisposal(Particle::OnePart *ip, 
                                             Particle::OnePart *jp,
                                             Particle::OnePart *kp)
{
  double State_prob,Fraction_Rot,Fraction_Vib;
  int i,numspecies,rotdof,vibdof,max_level,ivib,irot;
  double aveomega,pevib;

  Particle::OnePart *p;
  Particle::Species *species = particle->species;
  double AdjustFactor = 0.99999999;

  if (!kp) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    numspecies = 2;
    aveomega = 0.5*(params[ip->ispecies].omega + params[jp->ispecies].omega); 
  } else {
    ip->erot = 0.0;
    jp->erot = 0.0;
    kp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    kp->evib = 0.0;
    numspecies = 3;
    aveomega = (params[ip->ispecies].omega + params[jp->ispecies].omega + 
                params[kp->ispecies].omega)/3;
  }

  // handle each kind of energy disposal for non-reacting reactants
  // clean up memory for the products
  
  double E_Dispose = postcoln.etotal;

  for (i = 0; i < numspecies; i++) {
    if (i == 0) p = ip; 
    else if (i == 1) p = jp; 
    else p = kp;

    int sp = p->ispecies;
    rotdof = species[sp].rotdof;

    if (rotdof) {
      if (rotstyle == NONE) {
        p->erot = 0.0;
      } else if (rotdof == 2) {
        Fraction_Rot =
          1- pow(random->uniform(),(1/(2.5-aveomega)));
        p->erot = Fraction_Rot * E_Dispose;
        E_Dispose -= p->erot;
        
      } else if (rotdof > 2) {
        p->erot = E_Dispose * 
          sample_bl(random,0.5*species[sp].rotdof-1.0,
                    1.5-aveomega);
        E_Dispose -= p->erot;
      }
    }
    
    vibdof = species[sp].vibdof;

    if (vibdof) {
      if (vibstyle == NONE) {
        p->evib = 0.0;
      } else if (vibdof == 2 && vibstyle == DISCRETE) {
        max_level = static_cast<int> 
          (E_Dispose / (update->boltz * species[sp].vibtemp[0]));
        do {
          ivib = static_cast<int> 
            (random->uniform()*(max_level+AdjustFactor));
          p->evib = (double)
            (ivib * update->boltz * species[sp].vibtemp[0]);
          State_prob = pow((1.0 - p->evib / E_Dispose),
                           (1.5 - aveomega));
        } while (State_prob < random->uniform());
        E_Dispose -= p->evib;
        
      } else if (vibdof == 2 && vibstyle == SMOOTH) {
        Fraction_Vib =
          1.0 - pow(random->uniform(),(1.0 / (2.5-aveomega)));
        p->evib = Fraction_Vib * E_Dispose;
        E_Dispose -= p->evib;

      } else if (vibdof > 2 && vibstyle == SMOOTH) {
          p->evib = E_Dispose * 
          sample_bl(random,0.5*species[sp].vibdof-1.0,
                   1.5-aveomega);
          E_Dispose -= p->evib;
      } else if (vibdof > 2 && vibstyle == DISCRETE) {
          p->evib = 0.0;

          int nmode = particle->species[sp].nvibmode;
          int **vibmode = particle->eiarray[particle->ewhich[index_vibmode]];
          int pindex = p - particle->particles;

          for (int imode = 0; imode < nmode; imode++) {
            ivib = vibmode[pindex][imode];
            E_Dispose += ivib * update->boltz * 
            particle->species[sp].vibtemp[imode];
            max_level = static_cast<int>
            (E_Dispose / (update->boltz * species[sp].vibtemp[imode]));
            do {
              ivib = static_cast<int> 
              (random->uniform()*(max_level+AdjustFactor));
              pevib = ivib * update->boltz * species[sp].vibtemp[imode];
              State_prob = pow((1.0 - pevib / E_Dispose),
                               (1.5 - aveomega));
            } while (State_prob < random->uniform());

            vibmode[pindex][imode] = ivib;
            p->evib += pevib;
            E_Dispose -= pevib;
          }
        }
      }
    }
  
  // compute post-collision internal energies
  
  postcoln.erot = ip->erot + jp->erot;
  postcoln.evib = ip->evib + jp->evib;
  
  if (kp) {
    postcoln.erot += kp->erot;
    postcoln.evib += kp->evib;
  }
  
  // compute portion of energy left over for scattering
  
  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::sample_bl(RanPark *random, double Exp_1, double Exp_2)
{
  double Exp_s = Exp_1 + Exp_2;
  double x,y;
  do {
    x = random->uniform();
    y = pow(x*Exp_s/Exp_1, Exp_1)*pow((1.0-x)*Exp_s/Exp_2, Exp_2);
  } while (y < random->uniform());
  return x;
}

/* ----------------------------------------------------------------------
   compute a variable rotational relaxation parameter
------------------------------------------------------------------------- */

double CollideVSS::rotrel_serial(int isp, double Ec)
{
  double Tr = Ec /(update->boltz * (2.5-params[isp].omega));
  double rotphi = (1.0+params[isp].rotc2/sqrt(Tr) + params[isp].rotc3/Tr)
                / params[isp].rotc1; 
  return rotphi;
}

/* ----------------------------------------------------------------------
   compute a variable vibrational relaxation parameter
------------------------------------------------------------------------- */

double CollideVSS::vibrel_serial(int isp, double Ec)
{
  double Tr = Ec /(update->boltz * (3.5-params[isp].omega));
  double omega = params[isp].omega;
  double vibphi = 1.0 / (params[isp].vibc1/pow(Tr,omega) * 
                         exp(params[isp].vibc2/pow(Tr,1.0/3.0)));
  return vibphi;
}

/* ----------------------------------------------------------------------
   compute a variable rotational relaxation parameter
------------------------------------------------------------------------- */

double CollideVSS::rotrel_prohibdouble(int isp, double Ec)
{
  Particle::Species *species = particle->species;
  double rotphi = params[isp].rotc1 / (1.0 + tgamma(species[isp].rotdof
                       + 2.5-params[isp].omega)/tgamma(species[isp].rotdof+1.5-params[isp].omega)
                       * params[isp].rotc2*sqrt(update->boltz/Ec) + tgamma(species[isp].rotdof+2.5-params[isp].omega)
                       /tgamma(species[isp].rotdof+2.0-params[isp].omega)*params[isp].rotc3*update->boltz/Ec);
  return rotphi;
}

/* ----------------------------------------------------------------------
   compute a variable vibrational relaxation parameter
------------------------------------------------------------------------- */

double CollideVSS::vibrel_prohibdouble(int isp, double Ec, double doftot)
{
  Particle::Species *species = particle->species;
  double Tr = 2.0 * Ec / (update->boltz * doftot);
  double Zmw = 4.0*MY_PIS * pow(params[isp].diam,2.0)*pow(params[isp].tref,params[isp].omega-0.5)
                  * pow(Tr,-params[isp].omega)*101325.0*exp(params[isp].vibc1
                  * (pow(Tr,-1.0/3.0)-params[isp].vibc2)-18.42)/sqrt(species[isp].mass*update->boltz);
  double Zpark = 4.0*MY_PI*pow(params[isp].diam,2.0)*pow(params[isp].tref,params[isp].omega-0.5)
                          *pow(Tr,2.5-params[isp].omega)/(2.5e9*params[isp].park);
  double vibphi = Zmw + Zpark;
  return vibphi;
}

/* ----------------------------------------------------------------------
   read list of species defined in species file
   store info in filespecies and nfilespecies
   only invoked by proc 0
------------------------------------------------------------------------- */

void CollideVSS::read_param_file(char *fname)
{
  FILE *fp = fopen(fname,"r");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open VSS parameter file %s",fname);
    error->one(FLERR,str);
  }

  // set all diameters to -1, so can detect if not read

  for (int i = 0; i < nparams; i++) params[i].diam = -1.0;

  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines must have at least NWORDS 

  int NWORDS = 5;
  if (relaxflag == VARIABLE) {
    if (relaxtypeflag == PROHIBDOUBLE) NWORDS = 10;
    else NWORDS = 9;
  }
  char **words = new char*[NWORDS];
  char line[MAXLINE],copy[MAXLINE];
  int isp;

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    strcpy(copy,line);
    int nwords = wordcount(copy);
    if (nwords < NWORDS)
      error->one(FLERR,"Incorrect line format in VSS parameter file");
    wordparse(NWORDS,line,words);

    isp = particle->find_species(words[0]);
    if (isp < 0) continue;

    params[isp].diam = atof(words[1]);
    params[isp].omega = atof(words[2]);
    params[isp].tref = atof(words[3]);
    params[isp].alpha = atof(words[4]);
    if (relaxflag == VARIABLE) {
      params[isp].rotc1 = atof(words[5]);
      params[isp].rotc2 = atof(words[6]);
      params[isp].rotc3 =  (MY_PI+MY_PI2*MY_PI2)*params[isp].rotc2;
      params[isp].rotc2 =  (MY_PI*MY_PIS/2.)*sqrt(params[isp].rotc2);
      params[isp].vibc1 = atof(words[7]);
      params[isp].vibc2 = atof(words[8]);
      if (relaxtypeflag == PROHIBDOUBLE) params[isp].park = atof(words[9]);
    }
  }

  delete [] words;
  fclose(fp);
}

/* ----------------------------------------------------------------------
   count whitespace-delimited words in line
------------------------------------------------------------------------- */

int CollideVSS::wordcount(char *line)
{
  int nwords = 0;
  char *word = strtok(line," \t");
  while (word) {
    nwords++;
    word = strtok(NULL," \t");
  }
  return nwords;
}

/* ----------------------------------------------------------------------
   parse first N whitespace-delimited words in line
   store ptr to each word in words
------------------------------------------------------------------------- */

void CollideVSS::wordparse(int n, char *line, char **words)
{
  for (int i = 0; i < n; i++) {
    if (i == 0) words[i] = strtok(line," \t");
    else words[i] = strtok(NULL," \t");
  }
}

/* ----------------------------------------------------------------------
   return a per-species parameter to caller
------------------------------------------------------------------------- */

double CollideVSS::extract(int isp, const char *name)
{
  if (strcmp(name,"diam") == 0) return params[isp].diam;
  else if (strcmp(name,"omega") == 0) return params[isp].omega;
  else if (strcmp(name,"tref") == 0) return params[isp].tref;
  else error->all(FLERR,"Request for unknown parameter from collide");
  return 0.0;
}
