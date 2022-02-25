/***************************************************************************\
* Copyright (c) 2008, Claudio Pica, Martin Hansen                           *
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
#include "clover_tools.h"
#include "communications.h"

spinor_field_fund *Xs_fund=NULL;
spinor_field_fund *Ys_fund=NULL;
spinor_field_fund *eta_fund=NULL;
#ifdef UPDATE_EO
static spinor_field_fund *xi=NULL;
#endif 

void free_force_hmc_fund()
{
	free_spinor_field_f_fund(Xs_fund);
#ifdef UPDATE_EO
	free_spinor_field_f_fund(eta_fund);
	free_spinor_field_f_fund(xi);
#endif
}

void init_force_hmc_fund()
{
	static int init = 0;
	if(init == 0)
	{
#ifndef UPDATE_EO
		Xs_fund = alloc_spinor_field_f_fund(3, &glattice);
		Ys_fund = Xs_fund+1;
		eta_fund = Ys_fund+1;
#else
		Xs_fund = alloc_spinor_field_f_fund(2, &glattice);
		Ys_fund = Xs_fund+1;
		eta_fund = alloc_spinor_field_f_fund(1, &glat_even);
		xi = alloc_spinor_field_f_fund(1, &glat_odd);
		Xs_fund->type = &glat_even;
		Ys_fund->type = &glat_even;
#endif
		init = 1;
		atexit(free_force_hmc_fund);
	}
}

void force_hmc_fund(double dt, void *vpar)
{
	int n_iters = 0;
	force_hmc_fund_par *par = (force_hmc_fund_par*)vpar;
	suNg_av_field *force = *par->momenta;
	spinor_field_fund *pf = par->pf;

//	lprintf("force_hmc_fund pseudofermion", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",*pf,*(pf+1),*(pf+2),*(pf+3),*(pf+4),*(pf+5),*(pf+6),*(pf+7));
//	lprintf("force_hmc_fund force", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",*(force),*(force+1),*(force+2),*(force+3),*(force+4),*(force+5),*(force+6),*(force+7));
	init_force_hmc_fund();
	set_dirac_mass(par->mass);
	set_twisted_mass(par->mu);
	fermion_force_fund_begin();

	double tmp;
	mshift_par mpar;
	mpar.err2 = par->inv_err2;
	mpar.max_iter = 0;
	mpar.n = 1;
	mpar.shift = &tmp;
	mpar.shift[0] = 0;

        #ifdef MEASURE_FORCEHMC
                double forcestat[2] = {0.,0.};
        #else
                double *forcestat = NULL;
        #endif


#ifndef UPDATE_EO

	if(par->mu == 0 || par->hasenbusch != 0)
	{
		/* X = H^{-1} pf = D^{-1} g5 pf */
		spinor_field_zero_f_fund(Xs_fund);
		spinor_field_g5_assign_f_fund(pf);
		n_iters += g5QMR_mshift_fund(&mpar, &D_fund, pf, Xs_fund);
		spinor_field_g5_assign_f_fund(pf);

		if(par->hasenbusch == 0)
		{
			/* Y  D^{-1} ( g5 X ) */
			spinor_field_g5_f_fund(eta_fund, Xs_fund);
		}
		else if(par->hasenbusch == 1)
		{
			/* Y = H^{-1} ( g5 pf[k] + b X ) = D^{-1}( pf + b g5 X ) */
			spinor_field_g5_f_fund(eta_fund, Xs_fund);
			spinor_field_mul_f_fund(eta_fund, par->b, eta_fund);
			spinor_field_add_assign_f_fund(eta_fund, pf);
		}
		else if(par->hasenbusch == 2)
		{
			/* Y= -i D^{-1} ( pf[k] + imu g5 X )*/
			double mu1 = par->mu;
			double mu2 = par->mu + par->b;
			double muS = mu2*mu2 - mu1*mu1;
			spinor_field_g5_f_fund(eta_fund, Xs_fund);
			spinor_field_mul_f_fund(Xs_fund, muS, Xs_fund);
		}

		spinor_field_zero_f(Ys_fund);
		n_iters += g5QMR_mshift_fund(&mpar, &D_fund, eta_fund, Ys_fund);
	}
	else
	{
		n_iters += cg_mshift(&mpar, QpQm_tm_alt, pf, Xs_fund);/*not available*/
		Qtm_p_alt(Ys_fund,Xs_fund);
	}

#else

	if(par->mu == 0)
	{
		/* X_e = H^{-1} pf */
		/* X_o = D_{oe} X_e = D_{oe} H^{-1} pf */
		spinor_field_g5_assign_f_fund(pf);
		mre_fund_guess(&par->mpar, 0, Xs_fund, &D_fund, pf);
		n_iters += g5QMR_mshift_fund(&mpar, &D_fund, pf, Xs_fund);
		mre_fund_store(&par->mpar, 0, Xs_fund);
		spinor_field_g5_assign_f_fund(pf);/*up to here*/
//		lprintf("DEBUG", 10, "Here? %1.6f, %1.6f\n",(s1).c[0],(s1).c[9]);

//		lprintf("DEBUG", 10, "Here? %1.6f\n",*Xs_fund);
		/* Y_e = H^{-1} ( g5 pf + b X_e ) */
		/* Y_o = D_oe H^{-1} ( g5 pf + b X_e ) */
		if(par->hasenbusch != 1)
		{
			spinor_field_copy_f_fund(eta_fund, Xs_fund);
		}
		else
		{
			spinor_field_g5_f_fund(eta_fund, pf);
			spinor_field_mul_add_assign_f_fund(eta_fund, par->b, Xs_fund);
		}

		spinor_field_g5_assign_f_fund(eta_fund);
		mre_fund_guess(&par->mpar, 1, Ys_fund, &D_fund, eta_fund);
		n_iters += g5QMR_mshift_fund(&mpar, &D_fund, eta_fund, Ys_fund);
		mre_fund_store(&par->mpar, 1, Ys_fund);
		spinor_field_g5_assign_f_fund(eta_fund);

//		lprintf("DEBUG", 10, "Here? %1.6f\n",*eta_fund);
		if(par->hasenbusch == 2)
		{
			double muS = par->b*par->b;
			spinor_field_mul_f_fund(Xs_fund, muS, Xs_fund);
		}
	}
	else
	{
		/* Ye = 1/(QpQm+mu^2) \phi */
		mre_fund_guess(&par->mpar, 0, Ys_fund, QpQm_tm_alt, pf);
		n_iters += 2*cg_mshift(&mpar, QpQm_tm_alt, pf, Ys_fund);/*not avialable*/
		mre_fund_store(&par->mpar, 0, Ys_fund);
		Qtm_m_alt(Xs_fund, Ys_fund);

		if(par->hasenbusch == 2)
		{
			double mu1 = par->mu;
			double mu2 = par->mu + par->b;
			double muS = mu2*mu2 - mu1*mu1;
			spinor_field_mul_f_fund(Xs_fund, muS, Xs_fund);
		}
	}

#endif

//	lprintf("DEBUG", 10, "Here? %1.6f, %1.6f\n",*(Xs_fund+1),*Ys_fund);
	if(par->hasenbusch != 1)
	{
		force_fermion_core_fund(Xs_fund, Ys_fund, 1, dt, forcestat, 1.);
	}
	else
	{
		force_fermion_core_fund(Xs_fund, Ys_fund, 1, dt, forcestat, par->b);
	}

/*#ifdef WITH_CLOVER_EO

	if(par->logdet)
	{
		force_clover_logdet(par->mass, 2.); // 2 = # of flavors
	}

#endif*/

#ifdef MEASURE_FORCEHMC
    global_sum(forcestat,1);
    if (par->hasenbusch !=1){
        forcestat[0]*=fabs(dt)*(_FUND_NORM2/_FUND_NORM2)/((double)(GLB_T*GLB_X*GLB_Y*GLB_Z));
        global_max(forcestat+1,1);
        forcestat[1]*=fabs(dt)*(_FUND_NORM2/_FUND_NORM2);
    }
    else{
        forcestat[0]*=fabs(dt*par->b)*(_FUND_NORM2/_FUND_NORM2)/((double)(GLB_T*GLB_X*GLB_Y*GLB_Z));
        global_max(forcestat+1,1);
        forcestat[1]*=fabs(dt*par->b)*(_FUND_NORM2/_FUND_NORM2);
    }
    force_ave[par->id]+=forcestat[0];
    force_max[par->id]+=forcestat[1];
    lprintf("FORCE_HMC_FUND",10,"avr dt |force| = %1.8e dt maxforce = %1.8e, dt = %1.8e \n",forcestat[0],forcestat[1],dt);
    n_inv_iter[par->id-1]+=n_iters;
#endif


	fermion_force_fund_end(dt, force);


//	lprintf("force", 10, "Here? %1.6f, %1.6f, %1.6f\n",*force,*(force+1),*(force+2));
}
