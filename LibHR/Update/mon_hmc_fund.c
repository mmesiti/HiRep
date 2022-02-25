/***************************************************************************\
 * Copyright (c) 2015, Martin Hansen, Claudio Pica                        *
 * All rights reserved.                                                   *
 \***************************************************************************/

#include "global.h"
#include "update.h"
#include "logger.h"
#include "memory.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "clover_tools.h"
#include <stdlib.h>

static spinor_field_fund *tmp_pf = NULL;
static int mon_init = 1;

void hmc_fund_gaussian_pf(const struct _monomial *m)
{
	mon_hmc_fund_par *par = (mon_hmc_fund_par*)(m->data.par);
	gaussian_spinor_field_fund(par->pf);
//	lprintf("pseudofermion", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",*(par->pf),*(par->pf+1),*(par->pf+2),*(par->pf+3),*(par->pf+4),*(par->pf+5),*(par->pf+6),*(par->pf+7));
}

void hmc_fund_correct_pf(const struct _monomial *m)
{
	mon_hmc_fund_par *par = (mon_hmc_fund_par*)(m->data.par);
   
//	lprintf("pseudofermion", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",*(par->pf),*(par->pf+1),*(par->pf+2),*(par->pf+3),*(par->pf+4),*(par->pf+5),*(par->pf+6),*(par->pf+7));
	/* compute H2^{1/2}*pf = H*pf */
	spinor_field_copy_f_fund(tmp_pf, par->pf);
	set_dirac_mass(par->mass);
	H_fund(par->pf, tmp_pf);
//	lprintf("pseudofermion", 10, "%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n",par->mass,*(par->pf),*(par->pf+1),*(par->pf+2),*(par->pf+3),*(par->pf+4),*(par->pf+5),*(par->pf+6),*(par->pf+7));
}

void hmc_fund_correct_la_pf(const struct _monomial *m)
{
	mon_hmc_fund_par *par = (mon_hmc_fund_par*)(m->data.par);
	double shift;

	mshift_par mpar;
	mpar.err2 = m->data.MT_prec;
	mpar.max_iter = 0;
	mpar.n = 1;
	mpar.shift = &shift;
	mpar.shift[0] = 0;
   
	/* compute H2^{-1/2}*pf = H^{-1}*pf */
	spinor_field_g5_f_fund(tmp_pf, par->pf);
	set_dirac_mass(par->mass);
	spinor_field_zero_f_fund(par->pf); /* mshift inverter uses this as initial guess for 1 shift */
	g5QMR_mshift_fund(&mpar, &D_fund, tmp_pf, par->pf);/*g5QMR_mshift_fund?*/
}

const spinor_field_fund* hmc_fund_pseudofermion(const struct _monomial *m)
{
	mon_hmc_fund_par *par = (mon_hmc_fund_par*)(m->data.par);
	return par->pf;
}

void hmc_fund_add_local_action(const struct _monomial *m, scalar_field *loc_action)
{
	mon_hmc_fund_par *par = (mon_hmc_fund_par*)(m->data.par);
	pf_local_action_fund(loc_action, par->pf);
//#ifdef WITH_CLOVER_EO
//	clover_la_logdet(2., par->mass, loc_action);
//#endif
}

void hmc_fund_free(struct _monomial *m)
{
	mon_hmc_fund_par *par = (mon_hmc_fund_par*)m->data.par;

	if(par->pf != NULL)
	{
		free_spinor_field_f_fund(par->pf);
	}

	free(par);
	free(m);
}


struct _monomial* hmc_fund_create(const monomial_data *data)
{
	monomial *m = malloc(sizeof(*m));
	mon_hmc_fund_par *par = (mon_hmc_fund_par*)data->par;

	// Copy data structure
	m->data = *data;

	// Allocate memory for spinor field
	if(mon_init)
	{
		tmp_pf = alloc_spinor_field_f_fund(1, &glat_default);
		mon_init = 0;
	}
	par->pf = alloc_spinor_field_f_fund(1, &glat_default);

	// Setup force parameters
	par->fpar.id = data->id;
	par->fpar.n_pf = 1;
	par->fpar.pf = par->pf;
	par->fpar.inv_err2 = data->force_prec;
	par->fpar.inv_err2_flt = 1e-6;
	par->fpar.mass = par->mass;
	par->fpar.b = 0;
	par->fpar.hasenbusch = 0;
	par->fpar.mu = 0;
	par->fpar.logdet = 1;
	par->fpar.momenta = &suN_momenta;

	// Setup chronological inverter
	mre_fund_init(&(par->fpar.mpar), par->mre_past, data->force_prec);

	// Setup pointers to update functions
	m->free = &hmc_fund_free;
	m->update_force = &force_hmc_fund;
	m->force_par = &par->fpar;
	m->update_field = 0;
	m->field_par = 0;

	m->fund_pseudofermion = &hmc_fund_pseudofermion;
	m->gaussian_pf = &hmc_fund_gaussian_pf;
	m->correct_pf = &hmc_fund_correct_pf;
	m->correct_la_pf = &hmc_fund_correct_la_pf;
	m->add_local_action = &hmc_fund_add_local_action;

	return m;
}
