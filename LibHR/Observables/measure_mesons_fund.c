/***************************************************************************\
* Copyright (c) 2013 Rudy Arthur, Ari Hietanen                              *
*                                                                           *
*                                                                           *
\***************************************************************************/

#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "logger.h"
#include "io.h"
#include "random.h"
#include "communications.h"
#include "ranlux.h"
#include "gamma_spinor.h"
#include "spin_matrix_fund.h"
#include "propagator_fund.h"
#include <string.h>
#include "meson_observables.h"
#define PI 3.141592653589793238462643383279502884197


/*For fundamental flavors*/
static void spinmatrix_g_op(suNg_spin_matrix *out, suNg_spin_matrix *in, gamma_ind i){
  switch (i){
  case _g5: _spinmatrix_g_g5(*out,*in); break;
  case _id: *out=*in; break;
  case _g0: _spinmatrix_g_g0(*out,*in); break;
  case _g1: _spinmatrix_g_g1(*out,*in); break;
  case _g2: _spinmatrix_g_g2(*out,*in); break;
  case _g3: _spinmatrix_g_g3(*out,*in); break;
  case _g0g1: _spinmatrix_g_g0g1(*out,*in); break;
  case _g0g2: _spinmatrix_g_g0g2(*out,*in); break;
  case _g0g3: _spinmatrix_g_g0g3(*out,*in); break;
  case _g0g5: _spinmatrix_g_g0g5(*out,*in); break;
  case _g5g1: _spinmatrix_g_g5g1(*out,*in); break;
  case _g5g2: _spinmatrix_g_g5g2(*out,*in); break;
  case _g5g3: _spinmatrix_g_g5g3(*out,*in); break;
  case _g0g5g1: _spinmatrix_g_g5g0g1(*out,*in); break;
  case _g0g5g2: _spinmatrix_g_g5g0g2(*out,*in); break;
  case _g0g5g3: _spinmatrix_g_g5g0g3(*out,*in); break;
  default: break;
  }
}

static void op_spinmatrix_g(suNg_spin_matrix* out, suNg_spin_matrix* in, gamma_ind i){
  switch (i){
  case _g5: _g5_spinmatrix_g(*out,*in); break;
  case _id: *out=*in; break;
  case _g0: _g0_spinmatrix_g(*out,*in); break;
  case _g1: _g1_spinmatrix_g(*out,*in); break;
  case _g2: _g2_spinmatrix_g(*out,*in); break;
  case _g3: _g3_spinmatrix_g(*out,*in); break;
  case _g0g1: _g0g1_spinmatrix_g(*out,*in); break;
  case _g0g2: _g0g2_spinmatrix_g(*out,*in); break;
  case _g0g3: _g0g3_spinmatrix_g(*out,*in); break;
  case _g0g5: _g0g5_spinmatrix_g(*out,*in); break;
  case _g5g1: _g5g1_spinmatrix_g(*out,*in); break;
  case _g5g2: _g5g2_spinmatrix_g(*out,*in); break;
  case _g5g3: _g5g3_spinmatrix_g(*out,*in); break;
  case _g0g5g1: _g5g0g1_spinmatrix_g(*out,*in); break;
  case _g0g5g2: _g5g0g2_spinmatrix_g(*out,*in); break;
  case _g0g5g3: _g5g0g3_spinmatrix_g(*out,*in); break;
  default: break;
  }
}


static void op_propagator_fund(suNf_propagator_fund* out, suNf_propagator_fund* in, gamma_ind i){
  switch (i){
  case _g5: _g5_propagator_fund( (*out),(*in) ); break;
  case _id: *out=*in; break;
  case _g0: _g0_propagator_fund( (*out),(*in) ); break;
  case _g1: _g1_propagator_fund( (*out),(*in) ); break;
  case _g2: _g2_propagator_fund( (*out),(*in) ); break;
  case _g3: _g3_propagator_fund( (*out),(*in) ); break;
  default: break;
  }
}




#define corr_ind(px,py,pz,n_mom,tc,nm,cm) ((px)*(n_mom)*(n_mom)*(lt)*(nm)+(py)*(n_mom)*(lt)*(nm)+(pz)*(lt)*(nm)+ ((cm)*(lt)) +(tc))

/*For fundamental fermions*/
/* spinor_fields* are 4xnm arrays of spinor_field ordered([color][spinor])*/
void measure_mesons_core_fund(spinor_field_fund* psi0, spinor_field_fund* psi1, spinor_field_fund* eta, meson_observable* mo, int nm, int tau, int n_mom, int offset,int lt){
  int i,ix,t,x,y,z,beta,px,py,pz,tc;
  double pdotx,cpdotx,spdotx;
  complex tr;
  suNg_spin_matrix sma,smb, sm1,sm2,smtmp1,smtmp2,sm_src;
  meson_observable* motmp=mo;
  _spinmatrix_g_zero(smtmp1);
  _spinmatrix_g_zero(smtmp2);
  lprintf("measure_mesons_core",50,"Measuring channels: ");
  while (motmp!=NULL){
    lprintf("",50," %s",motmp->channel_name);
    motmp=motmp->next;
  }
  lprintf("",50,"\n");
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){
	for(i=0; i<nm; i++) {
	  for (t=0; t<T; t++) {	 
	    tc = (zerocoord[0]+t+GLB_T-tau)%GLB_T+offset;
	    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
		  ix=ipt(t,x,y,z);					
		  pdotx = 2.0*PI*(((double) px)*(x+zerocoord[1])/GLB_X + ((double) py)*(y+zerocoord[2])/GLB_Y + ((double) pz)*(z+zerocoord[3])/GLB_Z);
		  cpdotx=cos(pdotx);
		  spdotx=sin(pdotx);
		  for (beta=0;beta<4;beta++){ 
		    _spinmatrix_g_assign_row(sma, *_FIELD_AT(&psi0[beta*nm+i],ix), beta);
		    _spinmatrix_g_assign_row(smb, *_FIELD_AT(&psi1[beta*nm+i],ix), beta);
		    _spinmatrix_g_assign_row(sm_src, *_FIELD_AT(&eta[beta],ix), beta);
		  }
		  motmp=mo;
		  while (motmp!=NULL){
		    if (motmp->ind2!=_NOGAMMA){
		      spinmatrix_g_op(&smtmp1,&sma,motmp->ind2);
		      _spinmatrix_g_g5(sm1,smtmp1);
		      op_spinmatrix_g(&smtmp2,&smb,motmp->ind1);
		      _g5_spinmatrix_g(sm2,smtmp2);
		      _spinmatrix_g_mul_trace(tr,sm1,sm2);
		    }
		    else{
		      spinmatrix_g_op(&smtmp1,&sma,motmp->ind1);
		      _spinmatrix_g_mul_trace(tr,smtmp1,sm_src);
		    }
		    motmp->corr_re[corr_ind(px,py,pz,n_mom,tc,nm,i)] += motmp->sign*(tr.re*cpdotx+tr.im*spdotx);
		    motmp->corr_im[corr_ind(px,py,pz,n_mom,tc,nm,i)] += motmp->sign*(tr.im*cpdotx-tr.re*spdotx);
		    motmp=motmp->next;
		  }
		}
	  }
	}
      }
  lprintf("measure_mesons_core",50,"Measuring DONE! ");
}



static void measure_conserved_core_fund(spinor_field_fund* psi0, spinor_field_fund* psi1, spinor_field_fund* eta, meson_observable* mo, int nm, int tau, int n_mom, int offset,int lt){

  int i,ix,t,x,y,z,beta,px,py,pz,tc,a,ixmu;
  double pdotx,cpdotx,spdotx;
  complex tr;
  suNf_propagator_fund sp0,sp1,Usp,spf,sptmp1,sptmp2,sptmp3,spleft,spdag;
  suNg *u1;
  meson_observable* motmp=mo;

  lprintf("measure_conserved_core",50,"Measuring channels: ");
  while (motmp!=NULL){
    lprintf("",50," %s",motmp->channel_name);
    motmp=motmp->next;
  }
  lprintf("",50,"\n");
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){
	for(i=0; i<nm; i++) {
	  for (t=0; t<T; t++) {	 
	    tc = (zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T+offset;
	    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) {
		  ix=ipt(t,x,y,z);

					
		  pdotx = 2.0*PI*(((double) px)*(x+zerocoord[1])/GLB_X + ((double) py)*(y+zerocoord[2])/GLB_Y + ((double) pz)*(z+zerocoord[3])/GLB_Z);
		  cpdotx=cos(pdotx);
		  spdotx=sin(pdotx);
		  for (a=0;a<NG;++a){
		    for (beta=0;beta<4;beta++){ 
		      _propagator_fund_assign(sp0, *_FIELD_AT(&psi0[a*4*nm+beta*nm+i],ix),a,beta);
		      _propagator_fund_assign(sp1, *_FIELD_AT(&psi0[a*4*nm+beta*nm+i],ix),a,beta);
		    }
		  }
                  _propagator_fund_dagger(spdag,sp1);

		  motmp=mo;
		  while (motmp!=NULL){

			  ixmu = iup(ix,motmp->ind2);
			  for (a=0;a<NG;++a){
			    for (beta=0;beta<4;beta++){ 
			      _propagator_fund_assign(spf, *_FIELD_AT(&psi0[a*4*nm+beta*nm+i],ixmu), a,beta); //S(x+mu, y)
			    }
			  }
			  u1 = _4FIELD_AT(u_gauge_f_fund,ix,motmp->ind2);

			  // Tr [ g5 (1+g_mu) U^(x) S(x,0) g5 g_mu S^(x+mu, y) ]
			  _suNg_inverse_prop_multiply(Usp,*u1,sp0); //U^(x) S(x,0) 
			  sptmp1=Usp;
			  op_propagator_fund(&sptmp2,&sptmp1,motmp->ind1); //g_mu U^(x) S(x,0) 
			  _propagator_fund_add(sptmp1,sptmp1,sptmp2); //(1+g_mu) U^(x) S(x,0) 
			  _g5_propagator_fund(sptmp2,sptmp1); //g5 (1+g_mu) U^(x) S(x,0) 
			  _propagator_fund_dagger(sptmp1,spf); //S^(x+mu, 0)
			  _g5_propagator_fund(sptmp3,sptmp1); //g5 g_mu S^(x+mu, 0)
			  op_propagator_fund(&sptmp1,&sptmp3,motmp->ind1); //g_mu S^(x+mu, 0)
			  _propagator_fund_mul(spleft,sptmp2,sptmp1);	
			    

			  // Tr [ g5 (1+g_mu) U(x) S(x+mu,0) g5 g_mu S^(x, y) ]
			  _suNg_prop_multiply(Usp,*u1,spf); //U(x) S(x+mu,0)
			  sptmp1=Usp;
			  op_propagator_fund(&sptmp2,&sptmp1,motmp->ind1);//g_mu U(x) S(x,0) 
			  _propagator_fund_sub(sptmp1,sptmp1,sptmp2);//(1-g_mu) U(x) S(x,0) 
			  _g5_propagator_fund(sptmp2,sptmp1);//g5(1-g_mu) U(x) S(x,0) 
			  _g5_propagator_fund(sptmp1,spdag); //g5 S^(x, 0)
			  op_propagator_fund(&sptmp3,&sptmp1,motmp->ind1); //g_mu g5 S^(x, 0)
			  _propagator_fund_mul(sptmp1,sptmp2,sptmp3);

			  _propagator_fund_sub(sptmp2,spleft,sptmp1);
			  _propagator_fund_trace(tr,sptmp2);    

			  motmp->corr_re[corr_ind(px,py,pz,n_mom,tc,nm,i)] += 0.5*motmp->sign*(tr.re*cpdotx+tr.im*spdotx);
			  motmp->corr_im[corr_ind(px,py,pz,n_mom,tc,nm,i)] += 0.5*motmp->sign*(tr.im*cpdotx-tr.re*spdotx);
			  motmp=motmp->next;

		  } //END CORRELATOR LOOP
		} //END SPATIAL LOOP
	  } //END T LOOP
	} //END MASS LOOP
      } //END MOMENTUM LOOP
  lprintf("measure_formfactor_core",50,"Measuring DONE! ");
}

static void init_corrs(int nm, int n_mom, meson_observable* mo){
  int i,n_mom_tot,size;
  meson_observable* motmp=mo;
  if (n_mom>=GLB_X){
    n_mom=GLB_X-1;
    lprintf("measure_mesons_with_momenta",0,"Reduced n_mom to %d (no more momenta accessible with this lattice size) ",n_mom);
  }
  n_mom_tot = n_mom*n_mom*n_mom;
  size = GLB_T*nm*n_mom_tot;
  while (motmp!=NULL){
    if (size!=motmp->corr_size){
      free(motmp->corr_re);
      free(motmp->corr_im);
      motmp->corr_re=malloc(sizeof(double)*size);
      motmp->corr_im=malloc(sizeof(double)*size);    
      motmp->corr_size=size;
      for(i=0; i<motmp->corr_size; i++){
	motmp->corr_re[i] = 0.;
	motmp->corr_im[i] = 0.;
      }
    }
    motmp=motmp->next;
  }
}

static void zero_corrs(meson_observable* mo){
  meson_observable* motmp=mo;
  int i;
  while (motmp!=NULL){
    for(i=0; i<motmp->corr_size; i++){
      motmp->corr_re[i] = 0.;
      motmp->corr_im[i] = 0.;
    }
    motmp=motmp->next;
  }
}

/*For fundamadental flavors*/
void measure_mesons_fund(meson_observable* mo,spinor_field_fund* psi0, spinor_field_fund* eta, int nm, int tau){
  init_corrs(nm,1,mo);
  lprintf("measure_mesons",50,"measure default mesons");
  measure_mesons_core_fund(psi0, psi0, eta, mo,nm, tau, 1, 0,GLB_T);
}

void measure_conserved_currents_fund(meson_observable* mo,spinor_field_fund* psi0, spinor_field_fund* eta, int nm, int tau){
  init_corrs(nm,1,mo);
  lprintf("measure_mesons",50,"measure conserved vector currents");
  measure_conserved_core_fund(psi0, psi0, eta, mo,nm, tau, 1, 0,GLB_T);
}

void measure_point_mesons_momenta_fund(meson_observable* mo,spinor_field_fund* psi0, spinor_field_fund* eta, int nm, int tau, int n_mom){
  init_corrs(nm,n_mom,mo);
  lprintf("measure_mesons",50,"measure point mesons with momenta");
  measure_mesons_core_fund(psi0, psi0, eta, mo, nm, tau, n_mom, 0, GLB_T);
}

static void print_corr_core(meson_observable* mo,int lt,int conf, int nm, double* mass, char* label, int n_mom){
  int i,t,px,py,pz;
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){ 
	for(i=0; i<nm; i++) {
	  /*Real*/
	  if (n_mom>1){
	    lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s_re momentum(%d,%d,%d)= ",conf,mass[i],label,mo->channel_type,mo->channel_name,px,py,pz); 
	  }
	  else{
	    if (mo->ind1==mo->ind2){	    
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s= ",conf,mass[i],label,mo->channel_type,mo->channel_name); //To be compatible with the old output
	    }
	    else{
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s_re= ",conf,mass[i],label,mo->channel_type,mo->channel_name); 
	    }
	  }
	  for(t=0;t<lt;++t) lprintf("MAIN",0,"%e ",mo->corr_re[corr_ind(px,py,pz,n_mom,t,nm,i)]);
	  lprintf("MAIN",0,"\n");
	  /*Imaginary */
	  if (n_mom>1){
	    lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s_im momentum(%d,%d,%d)= ",conf,mass[i],label,mo->channel_type,mo->channel_name,px,py,pz); 
	  }
	  else{
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s_im= ",conf,mass[i],label,mo->channel_type,mo->channel_name); 
	  }
	  for(t=0;t<lt;++t) lprintf("MAIN",0,"%e ",mo->corr_im[corr_ind(px,py,pz,n_mom,t,nm,i)]);
	  lprintf("MAIN",0,"\n"); 
	}
      }
}


static void print_corr(meson_observable* mo,int lt,int conf, int nm, double* mass, char* label, int n_mom){
  meson_observable* motmp=mo;
  while (motmp!=NULL){  
    print_corr_core(motmp,lt,conf,nm,mass,label,n_mom);
    motmp=motmp->next;
  }
}

static void do_global_sum(meson_observable* mo, double norm){
  meson_observable* motmp=mo;
  int i;
  while (motmp!=NULL){
      global_sum(motmp->corr_re,motmp->corr_size);
      global_sum(motmp->corr_im,motmp->corr_size);
      for(i=0; i<motmp->corr_size; i++){
	motmp->corr_re[i] *= norm;
	motmp->corr_im[i] *= norm;
      }
    motmp=motmp->next;
  }
}

void print_mesons_fund(meson_observable* mo,double norm, int conf, int nm, double* mass, int lt, int n_mom, char* label){
  do_global_sum(mo,-((1./norm)/GLB_VOL3));
  print_corr(mo,lt,conf,nm,mass,label,n_mom);
  zero_corrs(mo);
}
