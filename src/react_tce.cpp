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
#include "react_tce.h"
#include "particle.h"
#include "collide.h"
#include "update.h"
#include "random_knuth.h"
#include "error.h"
#include "modify.h"
#include "compute.h"

using namespace SPARTA_NS;

enum{NONE,DISCRETE,SMOOTH};
enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // other files

/* ---------------------------------------------------------------------- */

ReactTCE::ReactTCE(SPARTA *sparta, int narg, char **arg) :
  ReactBird(sparta, narg, arg) {}

/* ---------------------------------------------------------------------- */

void ReactTCE::init()
{
  if (!collide || strcmp(collide->style,"vss") != 0)
    error->all(FLERR,"React tce can only be used with collide vss");

  ReactBird::init();
}

/* ---------------------------------------------------------------------- */

int ReactTCE::attempt(Particle::OnePart *ip, Particle::OnePart *jp,
                      double pre_etrans, double pre_erot, double pre_evib,
                      double &post_etotal, int &kspecies)
{
  double pre_etotal,ecc,e_excess,z;
  int imode;
  OneReaction *r;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  int icell = ip->icell;

  double pre_ave_rotdof = (species[isp].rotdof + species[jsp].rotdof)/2.0;

  int n = reactions[isp][jsp].n;
  if (n == 0) return 0;
  int *list = reactions[isp][jsp].list;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  double random_prob = random->uniform();
  double avei = 0.0;
  double z1 = 0.0;
  int nmode;

  // loop over possible reactions for these 2 species

  for (int i = 0; i < n; i++) {
    r = &rlist[list[i]];

    // ignore energetically impossible reactions

    pre_etotal = pre_etrans + pre_erot + pre_evib;

    // two options for total energy in TCE model
    // 0: partialEnergy = true: rDOF model
    // 1: partialEnergy = false: TCE: Rotation + Vibration

    // average DOFs participating in the reaction

    if (partialEnergy) {
       ecc = pre_etrans;
       z = r->coeff[0];
       if (pre_ave_rotdof > 0.1) ecc += pre_erot*z/pre_ave_rotdof;
    }
    else {
       ecc = pre_etotal;
       if (pre_etotal+r->coeff[4] <= 0.0) continue; // Cover cases where coeff[1].neq.coeff[4]
       z = pre_ave_rotdof;
       z1 = pre_ave_rotdof;
//       if (collide->vibstyle == SMOOTH)
       z += (species[isp].vibdof + species[jsp].vibdof)/2.0;
//       else if (collide->vibstyle == DISCRETE) {
//           if (species[isp].vibdof == 2) {
//               avei = ip->evib / (update->boltz * species[isp].vibtemp[0]);
//               if (avei > 0.0) z += avei * log(1 + 1.0/avei);
//               z1 += (species[isp].vibtemp[0]/temp[icell]) / (exp(species[isp].vibtemp[0]/temp[icell])-1);
//           }
//           else if (species[isp].vibdof > 2) {
//               nmode = species[isp].nvibmode;
//               imode = 0;
//               while (imode < nmode) {
//                   z += (species[isp].vibtemp[imode]/temp[icell]) / (exp(species[isp].vibtemp[imode]/temp[icell])-1);
//                   imode++;
//               }
//           }
//           if (species[jsp].vibdof == 2) {
//               avei = jp->evib / (update->boltz * species[jsp].vibtemp[0]);
//               if (avei > 0.0) z += avei * log(1 + 1.0/avei);
//               z1 += (species[jsp].vibtemp[0]/temp[icell]) / (exp(species[jsp].vibtemp[0]/temp[icell])-1);
//           }
//           else if (species[jsp].vibdof > 2) {
//               imode = 0;
//               while (imode < 4) {
//                   z += (species[jsp].vibtemp[imode]/temp[icell]) / (exp(species[jsp].vibtemp[imode]/temp[icell])-1);
//                   imode++;
//               }
//           }
//       }
//       printf("z_vib %f z_vib_temp %f evib %e ecc %e\n",z-2,z1-2,ip->evib+jp->evib,ecc);
    }

    e_excess = ecc - r->coeff[1];
    if (e_excess <= 0.0) continue;

    // compute probability of reaction

    switch (r->type) {
    case DISSOCIATION:
    case IONIZATION:
    case EXCHANGE:
      {
        react_prob += r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
          pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) *
          pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]);
        break;
      }

    case RECOMBINATION:
      {
        // skip if no 3rd particle chosen by Collide::collisions()
        //   this includes effect of boost factor to skip recomb reactions
        // check if this recomb reaction is the same one
        //   that the 3rd particle species maps to, else skip it
        // this effectively skips all recombinations reactions
        //   if selected a 3rd particle species that matches none of them
        // scale probability by boost factor to restore correct stats

        if (recomb_species < 0) continue;
        int *sp2recomb = reactions[isp][jsp].sp2recomb;
        if (sp2recomb[recomb_species] != list[i]) continue;

        react_prob += recomb_boost * recomb_density * r->coeff[2] *
          tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
          pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) *  // extended to general recombination case with non-zero activation energy
          pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]);
        break;
      }

    default:
      error->one(FLERR,"Unknown outcome in reaction");
      break;
    }

    // test against random number to see if this reaction occurs
    // if it does, reset species of I,J and optional K to product species
    // J particle is destroyed in recombination reaction, set species = -1
    // K particle can be created in a dissociation or ionization reaction,
    //   set its kspecies, parent will create it
    // important NOTE:
    //   does not matter what order I,J reactants are in compared
    //     to order the reactants are listed in the reaction file
    //   for two reasons:
    //   a) list of N possible reactions above includes all reactions
    //      that I,J species are in, regardless of order
    //   b) properties of pre-reaction state are stored in precoln:
    //      computed by setup_collision()
    //      used by perform_collision() after reaction has taken place
    //      precoln only stores combined properties of I,J
    //      nothing that is I-specific or J-specific

    if (react_prob > random_prob) {
      tally_reactions[list[i]]++;

      if (!computeChemRates) {
          ip->ispecies = r->products[0];

          switch (r->type) {
          case DISSOCIATION:
          case IONIZATION:
          case EXCHANGE:
            {
              jp->ispecies = r->products[1];
              break;
            }
          case RECOMBINATION:
            {
              // always destroy 2nd reactant species

              jp->ispecies = -1;
              break;
            }
          }

          if (r->nproduct > 2) kspecies = r->products[2];
          else kspecies = -1;

          post_etotal = pre_etotal + r->coeff[4];
      }

      return 1;
    }
  }

  return 0;
}
