/***************************************************************************\
* Copyright (c) 2008, Claudio Pica, Vincent Drach and Ari Hietanen          *
* Copyright (c) 2017, Martin Hansen                                         *
* All rights reserved.                                                      *
\***************************************************************************/

#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "dirac.h"
#include "inverters.h"
#include "rational_functions.h"
#include "representation.h"
#include "logger.h"
#include "linear_algebra.h"
#include "memory.h"
#include "communications.h"
#include "clover_tools.h"
#include "suN_repr_func_fund.h"
#include <string.h>

/* we need to compute  Tr  U(x,mu) g_5*(1-g_mu) chi2 # chi1^+
 * where # indicates the tensor product and Tr is the trace on Lorentz space.
 * the strategy is the following:
 * given the form of g_5(1-g_mu) one can compute only the first two lorentz
 * components of the spinor; so we first apply g_5(1-g_mu) to chi2 to find the first
 * two components; then we multiply these two vectors by U(x,mu) and
 * store the result in p.c[0], p.c[1]; when computing the trace we can factorize p.c[0] and p.c[1]
 * as they both multiply two components of chi1^+; we store these factors in p.c[2] and p.c[3].
 * the tensor product is performed by the macro
 * _suNf_FMAT(u,p): u = p.c[0] # p.c[2]^+ + p.c[1] # p.c[3]^+
 */

/* these macros use the variables ptmp, p (fundamental representation)*/
#ifdef BC_T_THETA
#define _T_theta_mulc(r) _vector_mulc_g(ptmp,eitheta[0],(r)); (r)=ptmp
#else
#define _T_theta_mulc(r)
#endif
#ifdef BC_X_THETA
#define _X_theta_mulc(r) _vector_mulc_g(ptmp,eitheta[1],(r)); (r)=ptmp
#else
#define _X_theta_mulc(r)
#endif
#ifdef BC_Y_THETA
#define _Y_theta_mulc(r) _vector_mulc_g(ptmp,eitheta[2],(r)); (r)=ptmp
#else
#define _Y_theta_mulc(r)
#endif
#ifdef BC_Z_THETA
#define _Z_theta_mulc(r) _vector_mulc_g(ptmp,eitheta[3],(r)); (r)=ptmp
#else
#define _Z_theta_mulc(r)
#endif

#define _F_fund_DIR0(u,chi1,chi2)				      \
  _vector_add_g(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNg_multiply(p.c[0],*(pu_gauge_f_fund(ix,0)),ptmp);		      \
  _T_theta_mulc(p.c[0]);                                      \
  _vector_add_g(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNg_multiply(p.c[1],*(pu_gauge_f_fund(ix,0)),ptmp);		      \
  _T_theta_mulc(p.c[1]);                                      \
  _vector_sub_g(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_sub_g(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNg_FMAT((u),p)

#define _F_fund_DIR1(u,chi1,chi2)				      \
  _vector_i_add_g(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNg_multiply(p.c[0],*(pu_gauge_f_fund(ix,1)),ptmp);		      \
  _X_theta_mulc(p.c[0]);                                      \
  _vector_i_add_g(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNg_multiply(p.c[1],*(pu_gauge_f_fund(ix,1)),ptmp);		      \
  _X_theta_mulc(p.c[1]);                                      \
  _vector_i_sub_g(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_i_sub_g(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNg_FMAT((u),p)

#define _F_fund_DIR2(u,chi1,chi2)				      \
  _vector_add_g(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNg_multiply(p.c[0],*(pu_gauge_f_fund(ix,2)),ptmp);		      \
  _Y_theta_mulc(p.c[0]);                                      \
  _vector_sub_g(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNg_multiply(p.c[1],*(pu_gauge_f_fund(ix,2)),ptmp);		      \
  _Y_theta_mulc(p.c[1]);                                      \
  _vector_sub_g(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_add_g(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNg_FMAT((u),p)

#define _F_fund_DIR3(u,chi1,chi2)				      \
  _vector_i_add_g(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNg_multiply(p.c[0],*(pu_gauge_f_fund(ix,3)),ptmp);		      \
  _Z_theta_mulc(p.c[0]);                                      \
  _vector_i_sub_g(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNg_multiply(p.c[1],*(pu_gauge_f_fund(ix,3)),ptmp);		      \
  _Z_theta_mulc(p.c[1]);                                      \
  _vector_i_sub_g(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_i_add_g(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNg_FMAT((u),p)


static suNg_av_field *force_sum = NULL;

/* ------------------------------------ */
/* CALCULATE FORCE OF THE HOPPING TERM (fundamental representation)  */
/* ------------------------------------ */
void force_fermion_core_fund(spinor_field_fund *Xs, spinor_field_fund *Ys, int auto_fill_odd, double dt, double* forcestat, double residue)
{
	double coeff;
	spinor_field_fund Xtmp, Ytmp;

	coeff = residue * dt * (_FUND_NORM2/_FUND_NORM2);
	Xtmp = *Xs;
	Ytmp = *Ys;
	Xs->type = &glattice;
	Ys->type = &glattice;



#ifdef UPDATE_EO

	if(auto_fill_odd)
	{
		spinor_field_fund Xe, Xo, Ye, Yo;

		Xe = *Xs;
		Xe.type = &glat_even;
		Xo = *Xs;
		Xo.ptr = Xs->ptr + glat_odd.master_shift;
		Xo.type = &glat_odd;

		Ye = *Ys;
		Ye.type = &glat_even;
		Yo = *Ys;
		Yo.type = &glat_odd;
		Yo.ptr = Ys->ptr + glat_odd.master_shift;

		Dphi_fund_(&Xo, &Xe);
		Dphi_fund_(&Yo, &Ye);
	}

	coeff = -coeff;

#endif

	// Communicate spinor field
	start_sf_sendrecv_fund(Xs);
	start_sf_sendrecv_fund(Ys);


	// Loop over lattice
	_PIECE_FOR(&glattice,xp)
	{
		suNg_algebra_vector f;
		suNg_vector ptmp;
		suNg_spinor p;
		suNg_FMAT s1;

		if(xp == glattice.inner_master_pieces)
		{
			_OMP_PRAGMA(master)
			{
				complete_sf_sendrecv_fund(Xs);
				complete_sf_sendrecv_fund(Ys);
			}
			_OMP_PRAGMA(barrier)
		}
#ifdef MEASURE_FORCEHMC
                _SITE_FOR_SUM(&glattice,xp,ix,forcestat[0],forcestat[1])
#else
                _SITE_FOR(&glattice,xp,ix)
#endif
		{
			int iy;
			suNg_spinor *chi1, *chi2;

			// Direction 0
			iy = iup(ix,0);
			_suNg_FMAT_zero(s1);
			chi1 = _FIELD_AT(Xs,ix);
			chi2 = _FIELD_AT(Ys,iy);
			_F_fund_DIR0(s1,chi1,chi2);
			chi1 = _FIELD_AT(Ys,ix);
			chi2 = _FIELD_AT(Xs,iy);
			_F_fund_DIR0(s1,chi1,chi2);

//			lprintf("FORCE_CORE", 10, "fermion force core %1.6f, %1.6f\n",(s1).c[0].im,(s1).c[10].im);
			_fund_algebra_project_FMAT(f,s1);
//			lprintf("FORCE_CORE", 10, "fermion force core %1.6f, %1.6f\n",(f).c[0],(f).c[9]);
//			lprintf("MD_INT", 10, "or or Here? %1.6f, %1.6f\n",(f).c[0],(f).c[9]);
//			lprintf("MD_INT", 10, "or or Here? %1.6f\n",coeff);
			_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force_sum,ix,0),coeff,f);

#ifdef MEASURE_FORCEHMC
                        double nsq;
                        _algebra_vector_sqnorm_g(nsq,f);
                        nsq=sqrt(nsq);
                        forcestat[0]+=nsq;
                        if (nsq>forcestat[1]) forcestat[1]=nsq;
#endif

			// Direction 1
			iy = iup(ix,1);
			_suNg_FMAT_zero(s1);
			chi1 = _FIELD_AT(Xs,ix);
			chi2 = _FIELD_AT(Ys,iy);
			_F_fund_DIR1(s1,chi1,chi2);
			chi1 = _FIELD_AT(Ys,ix);
			chi2 = _FIELD_AT(Xs,iy);
			_F_fund_DIR1(s1,chi1,chi2);

			_fund_algebra_project_FMAT(f,s1);
			_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force_sum,ix,1),coeff,f);


#ifdef MEASURE_FORCEHMC
                        _algebra_vector_sqnorm_g(nsq,f);
                        nsq=sqrt(nsq);
                        forcestat[0]+=nsq;
                        if (nsq>forcestat[1]) forcestat[1]=nsq;
#endif

			// Direction 2
			iy = iup(ix,2);
			_suNg_FMAT_zero(s1);
			chi1 = _FIELD_AT(Xs,ix);
			chi2 = _FIELD_AT(Ys,iy);
			_F_fund_DIR2(s1,chi1,chi2);
			chi1 = _FIELD_AT(Ys,ix);
			chi2 = _FIELD_AT(Xs,iy);
			_F_fund_DIR2(s1,chi1,chi2);

			_fund_algebra_project_FMAT(f,s1);
			_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force_sum,ix,2),coeff,f);

#ifdef MEASURE_FORCEHMC
                        _algebra_vector_sqnorm_g(nsq,f);
                        nsq=sqrt(nsq);
                        forcestat[0]+=nsq;
                        if (nsq>forcestat[1]) forcestat[1]=nsq;
#endif


			// Direction 3
			iy = iup(ix,3);
			_suNg_FMAT_zero(s1);
			chi1 = _FIELD_AT(Xs,ix);
			chi2 = _FIELD_AT(Ys,iy);
			_F_fund_DIR3(s1,chi1,chi2);
			chi1 = _FIELD_AT(Ys,ix);
			chi2 = _FIELD_AT(Xs,iy);
			_F_fund_DIR3(s1,chi1,chi2);

			_fund_algebra_project_FMAT(f,s1);
			_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force_sum,ix,3),coeff,f);

#ifdef MEASURE_FORCEHMC
                        _algebra_vector_sqnorm_g(nsq,f);
                        nsq=sqrt(nsq);
                        forcestat[0]+=nsq;
                        if (nsq>forcestat[1]) forcestat[1]=nsq;
#endif

		} // sites
	} // pieces

	// Reset spinor geometry
	Xs->type = Xtmp.type;
	Ys->type = Ytmp.type;
}

#undef _F_fund_DIR0
#undef _F_fund_DIR1
#undef _F_fund_DIR2
#undef _F_fund_DIR3

void fermion_force_fund_begin()
{
	if(force_sum == NULL)
	{
		force_sum = alloc_avfield(&glattice);
	}

	// Clear force field
	_MASTER_FOR(&glattice,ix)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			_algebra_vector_zero_g(*_4FIELD_AT(force_sum,ix,mu));
		}
	}

/*#ifdef WITH_CLOVER

	// Clear clover force field
	_MASTER_FOR(&glattice,ix)
	{
		for(int mu = 0; mu < 6; mu++)
		{
      #if defined(GAUGE_SPN) && defined(REPR_FUNDAMENTAL)
			_suNf_FMAT_zero(*_6FIELD_AT(cl_force,ix,mu));
      #else
      _suNf_zero(*_6FIELD_AT(cl_force,ix,mu));
      #endif
		}
	}

#endif*/
}

void fermion_force_fund_end(double dt, suNg_av_field *force)
{
/*#ifdef WITH_CLOVER

	// Evaluate derivative of clover term
	force_clover_core(dt);

#endif*/

#ifdef WITH_SMEARING

	// Evaluate smeared force and add to global force field
	smeared_gauge_force(force_sum, force);

#else

	// Add force to global force field
	_MASTER_FOR(&glattice,ix)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			_algebra_vector_add_assign_g(*_4FIELD_AT(force,ix,mu), *_4FIELD_AT(force_sum,ix,mu));
		}
	}

#endif

	// Boundary conditions
	apply_BCs_on_momentum_field(force);
}
