/***************************************************************************\
* Copyright (c) 2014 Vincent Drach, Ari Hietanen, Martin Hansen             *
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
#include "spin_matrix.h"
#include "propagator_fund.h"
#include "propagator.h"
#include <string.h>
#include "meson_observables.h"
#define PI 3.141592653589793238462643383279502884197
#define ONORM 0.7071067811865476

//void _propagator_diquark(complex dq)
//{
//	dq.re = 1.0;
//	dq.im = 0.0;
//}
//
//void _propagator_diquark_tmp(complex dq, suNf_propagator_fund SQ, suNf_propagator_fund SQtilde, int i, int ip, int j, int jp)
//{
//	complex tmp;
//	_complex_0(dq);
//	_complex_0(tmp);
//	int a,b,c,alpha,beta; 
//	int j,jp,k,kp;
//
//	double eas[NF][NG][NG];
/*
	for(int a = 0; a < NF; ++a)
	for(int b = 0; b < NG; ++b)
	for(int c = 0; c < NG; ++c)
	{
		eas[a][b][c] = 0;
	}

	eas[0][0][1] = ONORM;
	eas[0][1][0] = -ONORM;
	eas[1][0][3] = ONORM;
	eas[1][3][0] = -ONORM;
	eas[2][1][2] = ONORM;
	eas[2][2][1] = -ONORM;
	eas[3][0][2] = 0.5;
	eas[3][1][3] = -0.5;
	eas[3][2][0] = -0.5;
	eas[3][3][1] = 0.5;
	eas[4][2][3] = ONORM;
	eas[4][3][2] = -ONORM;


*/
// 	for(int j = 0; j < NG; ++j)
//	for(int jp = 0; jp < NG; ++jp)
// 	for(int k = 0; k < NG; ++k)
//	for(int kp = 0; kp < NG; ++kp)
//	{
//		for(int alpha = 0; alpha < 1; ++alpha)
//		for(int gamma = 0; gamma < 1; gamma++)
//		{
//			dq = SQ.c[k].c[alpha].c[gamma].c[kp];
//			tmp.re = SQ.c[j].c[alpha].c[gamma].c[jp].re; //Do I have to compute real and imaginary parts separately?
//			tmp.im = SQ.c[j].c[alpha].c[gamma].c[jp].im; //Do I have to compute real and imaginary parts separately?
//			tmp.re = eas[i][j][kp]*SQ.c[j].c[alpha].c[gamma].c[jp].re*eas[ip][jp][k]; //Do I have to compute real and imaginary parts separately?
//			tmp.im = eas[i][j][kp]*SQ.c[j].c[alpha].c[gamma].c[jp].im*eas[ip][jp][k]; //Do I have to compute real and imaginary parts separately?
//			_complex_mul_assign(dq, tmp, SQtilde.c[kp].c[alpha].c[gamma].c[k]);
//			_complex_mul_assign(dq, tmp, tmp);
//			_complex_mul_assign(dq, SQtilde.c[kp].c[alpha].c[gamma].c[k], Stilde.c[kp].c[alpha].c[gamma].c[k]);
//	dq.re = 1.0;
//	dq.im = 0.0;
//	_complex_mul_assign(dq, SQ.c[j].c[i].c[ip].c[jp], SQ.c[j].c[i].c[ip].c[jp]);
//		}
//	}
//}


void contract_chimera(spinor_field* chi0, spinor_field_fund* psi0, int tau)
{
	int ix, tc, idx;
	suNg Omega, Otr;
	suNf_propagator S, Stmp, Sp, Sm;
	suNf_propagator_fund SQ, SQtmp, SQtilde;
	complex dq, tmp;

//	complex corr_chimera[GLB_T][4][4];
//	complex corr_chimera[GLB_T];
//	complex corr_ps[GLB_T];
	complex corr_sp[GLB_T];
	complex corr_sm[GLB_T];

	lprintf("contract_chimera", 50, "Performing chimera baryon contraction ");
	error(NG != 4, 1, "contract_chimera [measure_chimera.c]", "The code does not work for that number of colors!\n");
	

	_symplectic(Omega);
	_suNg_transpose(Otr,Omega);


	// Zero variables
//	memset(corr_chimera, 0, sizeof(corr_chimera));
//	memset(corr_ps, 0, sizeof(corr_ps));
	memset(corr_sp, 0, sizeof(corr_sp));
	memset(corr_sm, 0, sizeof(corr_sm));

	double eas[NF][NG][NG];

	for(int a = 0; a < NF; ++a)
	for(int b = 0; b < NG; ++b)
	for(int c = 0; c < NG; ++c)
	{
		eas[a][b][c] = 0;
	}

	eas[0][0][1] = 0.7071067811865476;
	eas[0][1][0] = -0.7071067811865476;
	eas[1][0][3] = 0.7071067811865476;
	eas[1][3][0] = -0.7071067811865476;
	eas[2][1][2] = 0.7071067811865476;
	eas[2][2][1] = -0.7071067811865476;
	eas[3][0][2] = 0.5;
	eas[3][1][3] = -0.5;
	eas[3][2][0] = -0.5;
	eas[3][3][1] = 0.5;
	eas[4][2][3] = 0.7071067811865476;
	eas[4][3][2] = -0.7071067811865476;


	// LOOP VOLUME 
	for(int t = 0; t < T; t++)
	{
		tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;

		for(int x = 0; x < X; x++)
		for(int y = 0; y < Y; y++)
		for(int z = 0; z < Z; z++)
		{
			ix = ipt(t, x, y, z);

			// get S^ab_beta_gamma
			for(int a = 0; a < NF; ++a)
			for(int beta = 0; beta < 4; beta++)
			{
				idx = beta + a*4;
				_propagator_assign(S, *_FIELD_AT(&chi0[idx],ix), a, beta);
			}

			_g0_propagator(Stmp, S);

			_propagator_add(Sp, Stmp, S);
			_propagator_mul_assign(Sp, 0.5);


			_propagator_sub(Sm, Stmp, S);
			_propagator_mul_assign(Sm, 0.5);

//			corr_as[ts] = S.c[0].c[0].c[0].c[0];
//			_propagator_dagger(Stmp,S);

			// Why do we need to tranpose here?
//			_propagator_transpose(S, Stmp);

			// Stilde = Gamma S Gamma_tilde
			// C  =  g0 g2
			// C g5 =  g1 g3
			// Nucleon case: Stilde = -Cg5 S Cg5 = - g5g0g2 S g5g0g2
 			// we forgot about the minus sign that will contribute only to a global sign...

			// get S^ab_beta_gamma for fundamental
			for(int a = 0; a < NG; ++a)
			for(int beta = 0; beta < 4; beta++)
			{
				idx = beta + a*4;
				_propagator_fund_assign(SQ, *_FIELD_AT(&psi0[idx],ix), a, beta);
			}

			//multiply the color structure (Omega)
			_suNg_prop_multiply(SQtmp, Otr, SQ);
			_propagator_fund_transpose(SQ, SQtmp);
			_suNg_prop_multiply(SQtmp, Omega, SQ);
			_propagator_fund_transpose(SQ, SQtmp);

			//multiply the gamma structure
//			_propagator_fund_transpose(SQtr, SQ);
//			_g5g0g2_propagator(SQtmp, SQtr);
//			_propagator_g5g0g2(SQtr, SQtmp);

			_g5g0g2_propagator_fund(SQtmp, SQ);
			_propagator_fund_g5g0g2(SQtilde, SQtmp);

//			_complex_0(tmp);
//
//			for(int l = 0; l < NG; ++l)
//			for(int lp = 0; lp < NG; ++lp)
//		 	for(int k = 0; k < NG; ++k)
//			for(int kp = 0; kp < NG; ++kp)
//			{
//				for(int alpha = 0; alpha < 4; ++alpha)
//				for(int gamma = 0; gamma < 4; gamma++)
//				{
//					tmp.re = eas[0][l][kp]*SQ.c[l].c[alpha].c[gamma].c[lp].re*eas[0][lp][k]; //Do I have to compute real and imaginary parts separately?
//					tmp.im = eas[0][l][kp]*SQ.c[l].c[alpha].c[gamma].c[lp].im*eas[0][lp][k]; //Do I have to compute real and imaginary parts separately?
//					_complex_mul_assign(corr_dq[tc], tmp, SQtilde.c[kp].c[alpha].c[gamma].c[k]);
//				}
//			}
//
//			corr_as[tc].re += S.c[0].c[0].c[0].c[0].re;
//			corr_as[tc].im += S.c[0].c[0].c[0].c[0].im;
	
			for(int ip = 0; ip < NF; ++ip)
			for(int jp = 0; jp < NF; ++jp)
			{
//				_propagator_diquark(dq, SQ, SQtilde, ip, jp);

				_complex_0(dq);
				_complex_0(tmp);

			 	for(int l = 0; l < NG; ++l)
				for(int lp = 0; lp < NG; ++lp)
			 	for(int k = 0; k < NG; ++k)
				for(int kp = 0; kp < NG; ++kp)
				{
					for(int alpha = 0; alpha < 4; ++alpha)
					for(int gamma = 0; gamma < 4; gamma++)
					{
						tmp.re = eas[ip][l][kp]*SQ.c[l].c[alpha].c[gamma].c[lp].re*eas[jp][lp][k]; //Do I have to compute real and imaginary parts separately?
						tmp.im = eas[ip][l][kp]*SQ.c[l].c[alpha].c[gamma].c[lp].im*eas[jp][lp][k]; //Do I have to compute real and imaginary parts separately?
//						tmp.re = SQ.c[l].c[alpha].c[gamma].c[lp].re; //Do I have to compute real and imaginary parts separately?
//						tmp.im = SQ.c[l].c[alpha].c[gamma].c[lp].im; //Do I have to compute real and imaginary parts separately?
						_complex_mul_assign(dq, tmp, SQtilde.c[kp].c[alpha].c[gamma].c[k]);
//						_complex_mul_assign(dq, SQ.c[kp].c[alpha].c[gamma].c[k], SQtilde.c[kp].c[alpha].c[gamma].c[k]);
					}
				}

//				_complex_mul_assign(corr_chimera[tc], S.c[ip].c[0].c[0].c[jp], dq);

				
				for(int i = 0; i < 4; i++)
//				for(int j = 0; j < 4; j++)
				{
//					_propagator_diquark(dq);
//					corr_chimera[t][i][j] = dq;
//					_complex_mul_assign(corr_chimera[tc][i][j], dq, dq);
					_complex_mul_assign(corr_sp[tc], Sp.c[ip].c[i].c[i].c[jp], dq); //We use this for chimera.
					_complex_mul_assign(corr_sm[tc], Sm.c[ip].c[i].c[i].c[jp], dq); //We use this for chimera.
//					_complex_mul_assign(corr_ps[tc], S.c[ip].c[i].c[j].c[jp], Stmp.c[jp].c[j].c[i].c[ip]);
//					_complex_mul_assign(corr_chimera[tc][i][j], S.c[ip].c[i].c[j].c[jp], S.c[ip].c[i].c[j].c[jp]);
//					_complex_mul_assign(corr_chimera[tc][i][j], SQ.c[ip].c[i].c[j].c[jp], SQ.c[ip].c[i].c[j].c[jp]);
				}
				
			} // end cp
		} // END SPATIAL LOOP
	} // END TEMPORAL LOOP

//	global_sum((double*)corr_chimera, 2*4*4*GLB_T);
//	global_sum((double*)corr_ps, 2*GLB_T);
	global_sum((double*)corr_sp, 2*GLB_T);
	global_sum((double*)corr_sm, 2*GLB_T);
//	global_sum((double*)corr_chimera, 2*GLB_T);

	// Print chimera correlator
//	for(int t = 0; t < GLB_T; t++)
//	for(int i = 0; i < 4; i++)
//	{
//
//		lprintf("CORR_CHIMERA", 0, "%d %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e \n",
//				t,
//				corr_chimera[t][i][0].re,
//				corr_chimera[t][i][0].im,
//				corr_chimera[t][i][1].re,
//				corr_chimera[t][i][1].im,
//				corr_chimera[t][i][2].re,
//				corr_chimera[t][i][2].im,
//				corr_chimera[t][i][3].re,
//				corr_chimera[t][i][3].im
//		);
//
//	}
	
	for(int t = 0; t < GLB_T; t++)
	{
		lprintf("CORR_CHIMERA", 0, "%d %3.10e %3.10e %3.10e %3.10e \n",
				t,
//				corr_chimera[t].re,
				corr_sp[t].re,
				corr_sp[t].im,
//				corr_chimera[t].im,
				corr_sm[t].re,
				corr_sm[t].im
		);
	}

	lprintf("contract_chimera", 50, "Measuring DONE!\n");
}
