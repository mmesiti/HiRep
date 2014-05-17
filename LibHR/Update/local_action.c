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
#include "logger.h"

/*
 * compute the local action at every site for the HMC
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */
void local_hmc_action(local_action_type type,
                      action_par *par,
                      scalar_field *loc_action,
                      suNg_av_field *momenta,
                      spinor_field *phi1,
                      spinor_field *phi2) {

  
  /* check input types */
  _TWO_SPINORS_MATCHING(u_gauge,loc_action); /* check that action is defined on the global lattice */
  _TWO_SPINORS_MATCHING(loc_action,momenta);
  _TWO_SPINORS_MATCHING(loc_action,phi1);
  _TWO_SPINORS_MATCHING(loc_action,phi2);


  switch(type) {
  case NEW:
    _MASTER_FOR(&glattice,i) {
      *_FIELD_AT(loc_action,i)=0.;
    }
    break;
  case DELTA:
    _MASTER_FOR(&glattice,i) {
      *_FIELD_AT(loc_action,i)= -*_FIELD_AT(loc_action,i);
    }
    break;
  default:
    error(1,1,"local_hmc_action","Invalid type");
  }

  _MASTER_FOR(&glattice,i) {
    double a=0., tmp;
    /* Momenta */
    for (int j=0;j<4;++j) {
      suNg_algebra_vector *cmom=momenta->ptr+coord_to_index(i,j);
      _algebra_vector_sqnorm_g(tmp,*cmom); 
      a+=tmp; /* this must be positive */
    }
    a*=0.5*_FUND_NORM2;
    *_FIELD_AT(loc_action,i)+=a;
//  }
//  _MASTER_FOR(&glattice,i) {
    /* Gauge action */
    *_FIELD_AT(loc_action,i) += -(par->beta/((double)NG))*local_plaq(i);
  }

  /* pseudofermion fields can be defined only on even sites if EO preconditioning is used */
  if(par->n_pf > 0) {
    _MASTER_FOR(phi1->type,i) {
      double a=0., tmp;
    /* Fermions */
      for (int j=0;j<par->n_pf;++j) {
        _spinor_prod_re_f(tmp,*_FIELD_AT(&phi1[j],i),*_FIELD_AT(&phi2[j],i));
        a += tmp;
      }
  
      *_FIELD_AT(loc_action,i)+=a;
    }
  }
}
