/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "observables.h"
#include "update.h"
#include "suN.h"
#include "linear_algebra.h"
#include "error.h"
#include "representation.h"

extern rhmc_par _update_par;

/*
 * compute the local action at every site for the HMC
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */
void local_hmc_action(local_action_type type,
                      scalar_field *loc_action,
                      suNg_av_field *momenta,
                      spinor_field *phi1,
                      spinor_field *phi2) {

  _DECLARE_INT_ITERATOR(i);
  int j;
  double a,tmp;

  /* check input types */
#ifndef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(u_gauge,loc_action); /* check that action is defined on the global lattice */
  _TWO_SPINORS_MATCHING(loc_action,momenta);
  _TWO_SPINORS_MATCHING(loc_action,phi1);
  _TWO_SPINORS_MATCHING(loc_action,phi2);
#endif

  _MASTER_FOR(&glattice,i) {

    a=0.;

    /* Momenta */
    for (j=0;j<4;++j) {
      suNg_algebra_vector *cmom=momenta->ptr+coord_to_index(i,j);
      _algebra_vector_sqnorm_g(tmp,*cmom); 
      a+=tmp; /* this must be positive */
    }
    a*=0.5*_FUND_NORM2;

    /* Gauge action */
    a -= (_update_par.beta/((double)NG))*local_plaq(i);

    switch(type) {
    case NEW:
      *_FIELD_AT(loc_action,i)=a;
      break;
    case DELTA:
      a-=*_FIELD_AT(loc_action,i);
      *_FIELD_AT(loc_action,i)=a;
      break;
    default:
      error(1,1,"local_hmc_action","Invalid type");
    }

  }

  /* pseudofermion fields can be defined only on even sites is the preconditioning is used */
  _MASTER_FOR(phi1->type,i) {
    a=0.;
  /* Fermions */
    for (j=0;j<_update_par.n_pf;++j) {
      _spinor_prod_re_f(tmp,*_FIELD_AT(&phi1[j],i),*_FIELD_AT(&phi2[j],i));
      a += tmp;
    }

    *_FIELD_AT(loc_action,i)+=a;
  }

   
}
