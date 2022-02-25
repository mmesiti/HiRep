/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
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

#include <math.h>

static unsigned int n_ws = 0;
static spinor_field *chi = NULL;
static spinor_field *Hchi = NULL;

void force_rhmc(double dt, void *vpar)
{
	#ifdef TIMING
	struct timeval start, end;
	struct timeval start1, end1;
	struct timeval etime;

	#ifdef TIMING_WITH_BARRIERS
	MPI_Barrier(GLB_COMM);
	#endif
	gettimeofday(&start,0);
	#endif

	int n_iters = 0;
	force_rhmc_par *par = (force_rhmc_par*)vpar;
	suNg_av_field *force = *par->momenta;
	spinor_field *pf = par->pf;
	rational_app *ratio = par->ratio;
	fermion_force_begin();

	mshift_par mpar;
	mpar.n = ratio->order;
	mpar.shift = ratio->b;
	mpar.err2 = par->inv_err2; /* this should be high for reversibility */
	mpar.max_iter = 0; /* no limit */

        #ifdef MEASURE_FORCEHMC
                double forcestat[2] = {0.,0.};
        #else
                double *forcestat = NULL;
        #endif


	if(n_ws < ratio->order+1)
	{
		if(chi != NULL) free_spinor_field_f(chi);
		n_ws = ratio->order + 1;
		chi = alloc_spinor_field_f(n_ws, &glattice);
		Hchi = chi + n_ws - 1;
	}

	#ifdef UPDATE_EO
	for(int n = 0; n < ratio->order; n++)
	{
		chi[n].type = &glat_even;
	}
	Hchi->type = &glat_even;
	#endif

	for(int k = 0; k < par->n_pf; k++)
	{
		/* compute inverse vectors chi[i] = (H^2 - b[i])^1 * pf */
		if (mpar.n == 1) { spinor_field_zero_f(chi); }

		#ifdef TIMING
		#ifdef TIMING_WITH_BARRIERS
		MPI_Barrier(GLB_COMM);
		#endif
		gettimeofday(&start1,0);
		#endif

		set_dirac_mass(par->mass);
		n_iters += cg_mshift(&mpar, &H2, &pf[k], chi);

		#ifdef TIMING
		#ifdef TIMING_WITH_BARRIERS
		MPI_Barrier(GLB_COMM);
		#endif
		gettimeofday(&end1,0);
		timeval_subtract(&etime,&end1,&start1);
		lprintf("TIMING",0,"cg_mshift in force_rhmc %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
		#endif

		for(int n = 0; n < ratio->order; n++)
		{
			H(Hchi, &chi[n]);
			force_fermion_core(Hchi, &chi[n], 1, dt, forcestat, ratio->a[n+1]);
#ifdef MEASURE_FORCEHMC
			global_sum(forcestat,1);
			forcestat[0]*=fabs(dt*ratio->a[n+1])*(_REPR_NORM2/_FUND_NORM2)/((double)(GLB_T*GLB_X*GLB_Y*GLB_Z));
			global_max(forcestat+1,1);
			forcestat[1]*=fabs(dt*ratio->a[n+1])*(_REPR_NORM2/_FUND_NORM2);
			force_ave[par->id]+=forcestat[0];
			force_max[par->id]+=forcestat[1];
			lprintf("FORCE_RHMC",10,"avr dt |force| = %1.8e dt maxforce = %1.8e, dt = %1.8e \n",forcestat[0],forcestat[1],dt);
//			n_inv_iter[par->id-1]+=n_iters;
#endif

		}
	}

	#ifdef TIMING
	#ifdef TIMING_WITH_BARRIERS
	MPI_Barrier(GLB_COMM);
	#endif
	gettimeofday(&end,0);
	timeval_subtract(&etime,&end,&start);
	lprintf("TIMING",0,"force_rhmc %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
	#endif

#ifdef WITH_CLOVER_EO
	double nf = (-2.0*ratio->n)/ratio->d;
	force_clover_logdet(par->mass, nf);
#endif

	fermion_force_end(dt, force);
}
