/*******************************************************************************
*
* Wrapper functions for different type of measurements 
*
*******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "spectrum.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "spin_matrix_fund.h"
#include "propagator.h"
#include "propagator_fund.h"
#include "gaugefix.h"

void measure_chimera(double* m, double* mf, int conf_num, double precision){
	// declare point sources and props  
	spinor_field* source = alloc_spinor_field_f(4*NF,&glattice); //This isn't glat_even so that the odd sites will be set to zero explicitly
	spinor_field* prop =  alloc_spinor_field_f(4*NF,&glattice);
	spinor_field_fund* source_fund = alloc_spinor_field_f_fund(4*NG,&glattice); //This isn't glat_even so that the odd sites will be set to zero explicitly
	spinor_field_fund* prop_fund =  alloc_spinor_field_f_fund(4*NG,&glattice);	

	int nm=1;
	int tau=0;
	int k,beta;

	init_propagator_eo(nm, m, precision);
	init_propagator_fund_eo(nm, mf, precision);

	// create point source
	create_full_point_source(source,tau);
	create_full_point_source_fund(source_fund,tau);

	// to test gauge invariance of the correlator  : 
	//	random_gauge_transform(u_gauge);
  //  represent_gauge_field();

	// calc invert
	calc_propagator(prop,source,4*NF);//4x3 for QCD 
	calc_propagator_fund(prop_fund,source_fund,4*NG);//4x3 for QCD 

	// perform contraction
	contract_chimera(prop,prop_fund,tau);
  
  // one should run the meson contraction as well here 
	//for (k=0;k<NF;++k){
	//	for (beta=0;beta<4;++beta){
	//		spinor_field_zero_f(&source2[beta]);
	//		spinor_field_zero_f(&prop2[beta]);
	//	}

	//		measure_mesons(prop2, source2, nm, tau);
	//}
	//print_mesons(1.,conf_num,0,m,"DEFAULT_TEST_POINT");

	// free 
	free_propagator_eo(); 
	free_propagator_fund_eo(); 
	free_spinor_field_f(source);
	free_spinor_field_f_fund(source_fund);
	free_spinor_field_f(prop);
	free_spinor_field_f_fund(prop_fund);
}
