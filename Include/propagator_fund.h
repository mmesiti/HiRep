/*******************************************************************************
*
* File propagator_fund.h
*
* Type definitions and macros for propagator
*
* 2013 Rudy Arthur, Ari Hietanen
*
*******************************************************************************/

#ifndef PROPAGATOR_FUND_H
#define PROPAGATOR_FUND_H
#include "suN_types.h"
#include "suN.h"
#include "spin_matrix_fund.h"
#include "gamma_spinor.h"

typedef struct
{
  suNg_spin_matrix c[NG]; 
} suNf_propagator_fund;


//r propagator, s spinor, i color index, j spin index
#define _propagator_fund_assign(p, s, i, j) \
do {  \
  int ITMP; \
  for(ITMP=0;ITMP<NG;++ITMP){ \
  (p).c[ITMP].c[0].c[j].c[i] = (s).c[0].c[ITMP]; \
  (p).c[ITMP].c[1].c[j].c[i] = (s).c[1].c[ITMP]; \
  (p).c[ITMP].c[2].c[j].c[i] = (s).c[2].c[ITMP]; \
  (p).c[ITMP].c[3].c[j].c[i] = (s).c[3].c[ITMP]; \
 } \
} while(0) 
  

  //(r).c[i].c[j] = (s)
  //_spinmatrix_g_assign_col(r.c[i], s, j)

//r propagator, s spinmatrix_g, i color index
#define _propagator_fund_assign_spin_matrix(r, s, i) \
  (r).c[i] = (s)

//r propagator
#define _PROP_AT(r,a,alpha,beta,b) ( (r).c[(a)].c[(alpha)].c[(beta)].c[(b)] )
    
#define _PROP_IDX(r,i,j) ( (r).c[((i)/4)].c[((i)%4)].c[((j)%4)].c[((j)/4)] )

//  r=0  (r propagator) 
#define _propagator_fund_zero(r) \
do {  \
  int ITMP; \
  for(ITMP=0;ITMP<NG;++ITMP){ _spinmatrix_g_zero((r).c[ITMP]); } \
} while(0) 

//r propagator k result; Tr [ r ]
#define _propagator_fund_one(r) \
   _propagator_fund_zero(r); \
   do { \
      int ITMP; \
      for(ITMP=0;ITMP<NG;ITMP++){ \
	_complex_1( r.c[ITMP].c[0].c[0].c[ITMP] ); \
	_complex_1( r.c[ITMP].c[1].c[1].c[ITMP] ); \
	_complex_1( r.c[ITMP].c[2].c[2].c[ITMP] ); \
	_complex_1( r.c[ITMP].c[3].c[3].c[ITMP] ); \
      } \
   } while(0) 



//r propagator k result; Tr [ r ]
#define _propagator_fund_trace(k, r) \
   do { \
      (k).re=0.;(k).im=0.; \
      int ITMP; \
      for(ITMP=0;ITMP<NG;ITMP++){ \
	(k).re += (r).c[ITMP].c[0].c[0].c[ITMP].re; \
	(k).im += (r).c[ITMP].c[0].c[0].c[ITMP].im; \
	(k).re += (r).c[ITMP].c[1].c[1].c[ITMP].re; \
	(k).im += (r).c[ITMP].c[1].c[1].c[ITMP].im; \
	(k).re += (r).c[ITMP].c[2].c[2].c[ITMP].re; \
	(k).im += (r).c[ITMP].c[2].c[2].c[ITMP].im; \
	(k).re += (r).c[ITMP].c[3].c[3].c[ITMP].re; \
	(k).im += (r).c[ITMP].c[3].c[3].c[ITMP].im; \
      } \
   } while(0) 


#define _propagator_fund_add(p,q,r)		\
   do { \
     int a,beta;			  \
     for (a=0;a<NG;++a) for (beta=0;beta<4;++beta){			\
	 _spinor_add_g(p.c[a].c[beta],q.c[a].c[beta],r.c[a].c[beta]);	\
       }								\
   } while(0) 

#define _propagator_fund_sub(p,q,r)		\
   do { \
     int a,beta;			  \
     for (a=0;a<NG;++a) for (beta=0;beta<4;++beta){			\
	 _spinor_sub_g(p.c[a].c[beta],q.c[a].c[beta],r.c[a].c[beta]);	\
       }								\
   } while(0) 


//S propagator k factor;
#define _propagator_fund_mul_assign(S, k) \
   do { \
     int a,beta;			  \
     for (a=0;a<NG;++a) for (beta=0;beta<4;++beta){			\
	 _spinor_mul_f((S).c[a].c[beta],k,(S).c[a].c[beta]);	\
       }								\
   } while(0) 

#define _propagator_fund_mulc_assign(S, k) \
   do { \
      int ITMP, JTMP; \
      for(ITMP=0;ITMP<4*NG;ITMP++){ \
      for(JTMP=0;JTMP<4*NG;JTMP++){ \
	complex tmp; tmp = _PROP_IDX( S,ITMP,JTMP); \
	_complex_mul( _PROP_IDX( S,ITMP,JTMP), k, tmp ); \
      }} \
   } while(0)  
	

//r propagator = s^T
#define _propagator_fund_transpose(r,s) \
   do { \
      int ITMP, JTMP; \
      for(ITMP=0;ITMP<4*NG;ITMP++){ \
      for(JTMP=0;JTMP<4*NG;JTMP++){ \
	_PROP_IDX( (r),ITMP,JTMP) = _PROP_IDX( (s),JTMP,ITMP); \
      }} \
   } while(0)

//r propagator = s^dagger
#define _propagator_fund_dagger(r,s) \
   do { \
      int ITMP, JTMP; \
      for(ITMP=0;ITMP<4*NG;ITMP++){ \
      for(JTMP=0;JTMP<4*NG;JTMP++){ \
	_PROP_IDX( (r),ITMP,JTMP).re = _PROP_IDX( (s),JTMP,ITMP).re; \
	_PROP_IDX( (r),ITMP,JTMP).im = -_PROP_IDX( (s),JTMP,ITMP).im; \
      }} \
   } while(0) 

//s spinor = r propagator t spinor
#define _propagator_fund_mul_spinor(s,r,t) \
   do { \
      _spinor_zero_f(s); \
      int a,b,alpha,beta; \
      for(alpha=0;alpha<NG;++alpha){ \
      for(a=0;a<NG;++a){ \
      for(beta=0;beta<NG;++beta){ \
      for(b=0;b<NG;++b){ \
	_complex_mul_assign((s).c[(alpha)].c[(a)],(r).c[(a)].c[(alpha)].c[(beta)].c[(b)],(t).c[(beta)].c[(b)]); \
      }}}} \
   } while(0) 

//s spinor =  t^dagger spinor r propagator
#define _propagator_fund_leftmul_spinor(s,t,r) \
   do { \
      _spinor_zero_f(s); \
      int a,b,alpha,beta; \
      for(beta=0;beta<NG;++beta){ \
      for(b=0;b<NG;++b){ \
      for(alpha=0;alpha<NG;++alpha){ \
      for(a=0;a<NG;++a){ \
	_complex_mul_star_assign((s).c[(beta)].c[(b)],(t).c[(alpha)].c[(a)],(r).c[(a)].c[(alpha)].c[(beta)].c[(b)]); \
      }}}} \
   } while(0) 

//Q propagator S propagator R propagator factor; Q = SR
#define _propagator_fund_mul(Q,S,R) \
do { \
      int a,b,c, alpha, beta, gamma; \
      for(a=0;a<NG;++a){ \
      for(b=0;b<NG;++b){ \
      for(alpha=0;alpha<4;++alpha){ \
      for(beta=0;beta<4;++beta){ \
	_complex_0( Q.c[a].c[alpha].c[beta].c[b] ); \
      for(c=0;c<NG;c++){ \
      for(gamma=0;gamma<4;gamma++){ \
	_complex_mul_assign( Q.c[a].c[alpha].c[beta].c[b], S.c[a].c[alpha].c[gamma].c[c], R.c[c].c[gamma].c[beta].c[b] ); \
      }}}}}} \
   } while(0) 

//tr complex S propagator R propagator factor; tr = Trace[ S R ]
#define _propagator_fund_mul_trace(tr,S,R) \
   do { \
     int a,b, alpha, beta; \
     (tr).re=(tr).im=0.0;\
     for(a=0;a<NG;++a){ \
     for(alpha=0;alpha<4;++alpha){ \
     for(b=0;b<NG;++b){ \
     for(beta=0;beta<4;++beta){ \
       _complex_mul_assign( tr, S.c[a].c[alpha].c[beta].c[b], R.c[b].c[beta].c[alpha].c[a] ); \
     }}}} \
   } while(0)

//tr complex S propagator R propagator factor; tr = Trace[ S^dagger R ]
#define _propagator_fund_muldag_trace(tr,S,R) \
   do { \
     int a,b, alpha, beta; \
     (tr).re=(tr).im=0.0;\
     for(a=0;a<NG;++a){ \
     for(alpha=0;alpha<4;++alpha){ \
     for(b=0;b<NG;++b){ \
     for(beta=0;beta<4;++beta){ \
       _complex_prod_assign( tr, S.c[b].c[beta].c[alpha].c[a], R.c[b].c[beta].c[alpha].c[a] ); \
     }}}} \
   } while(0)

//Color matrix U x propagator
#define _suNg_prop_multiply(us,u,s)\
do {				   \
  suNg_vector v1,v2;		   \
  int a,b,alpha,beta; \
  for (beta=0;beta<4;++beta) for (alpha=0;alpha<4;++alpha){	\
    for (b=0;b<NG;++b){							\
      for (a=0;a<NG;++a){							\
        v1.c[a]=(s).c[a].c[alpha].c[beta].c[b];				\
      }									\
      _suNg_multiply(v2,(u),v1);					\
      for (a=0;a<NG;++a){						\
        (us).c[a].c[alpha].c[beta].c[b]=v2.c[a];				\
      }\
    }\
  }\
} while(0)

//Color matrix U^dag x propagator
#define _suNg_inverse_prop_multiply(us,u,s)\
do {				   \
  suNg_vector v1,v2;		   \
  int a,b,alpha,beta; \
  for (beta=0;beta<4;++beta) for (alpha=0;alpha<4;++alpha){	\
    for (b=0;b<NG;++b){							\
      for (a=0;a<NG;++a){							\
        v1.c[a]=(s).c[a].c[alpha].c[beta].c[b];				\
      }									\
      _suNg_inverse_multiply(v2,(u),v1);					\
      for (a=0;a<NG;++a){						\
        (us).c[a].c[alpha].c[beta].c[b]=v2.c[a];				\
      }									\
    }						\
  }						\
} while(0)

#define _id_propagator_fund(p,q)			\
 do { \
      int ITMP, JTMP; \
      for(ITMP=0;ITMP<4*NG;ITMP++){ \
      for(JTMP=0;JTMP<4*NG;JTMP++){ \
	_PROP_IDX( (p),ITMP,JTMP) = _PROP_IDX( (q),ITMP,JTMP); \
      }} \
   } while(0) 


//P = Gamma Q 
#define _g0_propagator_fund(p,q)			\
do{\
  int a; \
  for (a=0;a<NG;++a){ _g0_spinmatrix_g((p).c[a],(q).c[a]); }	\
} while(0)

#define _g1_propagator_fund(p,q)			\
do{\
  int a; \
  for (a=0;a<NG;++a){ _g1_spinmatrix_g((p).c[a],(q).c[a]); }	\
} while(0)

#define _g2_propagator_fund(p,q)			\
do{\
  int a; \
  for (a=0;a<NG;++a){ _g2_spinmatrix_g((p).c[a],(q).c[a]); }	\
} while(0)

#define _g3_propagator_fund(p,q)			\
do{\
  int a; \
  for (a=0;a<NG;++a){ _g3_spinmatrix_g((p).c[a],(q).c[a]); }	\
} while(0)

#define _g5_propagator_fund(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _g5_spinmatrix_g((p).c[a],(q).c[a]); }  \
} while(0)

#define _g5g0_propagator_fund(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _g5g0_spinmatrix_g((p).c[a],(q).c[a]); }  \
} while(0)

#define _g5g3_propagator_fund(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _g5g3_spinmatrix_g((p).c[a],(q).c[a]); }  \
} while(0)

#define _g5g1_propagator_fund(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _g5g1_spinmatrix_g((p).c[a],(q).c[a]); }  \
} while(0)

#define _g5g2_propagator_fund(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _g5g2_spinmatrix_g(p.c[a],q.c[a]); }  \
} while(0)


//P = Q Gamma
#define _propagator_fund_g0(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _spinmatrix_g_g0((p).c[a],(q).c[a]); }  \
} while(0)

#define _propagator_fund_g2(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _spinmatrix_g_g2((p).c[a],(q).c[a]); }  \
} while(0)

#define _propagator_fund_g5(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _spinmatrix_g_g5((p).c[a],(q).c[a]); }  \
} while(0)

#define _propagator_fund_g5g0(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _spinmatrix_g_g5g0((p).c[a],(q).c[a]); }  \
} while(0)

#define _g5g0g1_propagator_fund(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _g5g0g1_spinmatrix_g(p.c[a],q.c[a]); }  \
} while(0)


#define _g5g0g2_propagator_fund(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _g5g0g2_spinmatrix_g((p).c[a],(q).c[a]); }  \
} while(0)


#define _propagator_fund_g5g0g2(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _spinmatrix_g_g5g0g2((p).c[a],(q).c[a]); }  \
} while(0)

#define _g5g0g3_propagator_fund(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _g5g0g3_spinmatrix_g(p.c[a],q.c[a]); }  \
} while(0)


#define _propagator_fund_g5g3(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _spinmatrix_g_g5g3((p).c[a],(q).c[a]); }  \
} while(0)

#define _propagator_fund_g5g1(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _spinmatrix_g_g5g1((p).c[a],(q).c[a]); }  \
} while(0)


#define _g0g1_propagator_fund(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _g0g1_spinmatrix_g(p.c[a],q.c[a]); }  \
} while(0)


#define _g0g2_propagator_fund(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _g0g2_spinmatrix_g((p).c[a],(q).c[a]); }  \
} while(0)

#define _propagator_fund_g0g2(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _spinmatrix_g_g0g2((p).c[a],(q).c[a]); }  \
} while(0)

#define _g0g3_propagator_fund(p,q)			\
do{						\
  int a; \
  for (a=0;a<NG;++a){ _g0g3_spinmatrix_g(p.c[a],q.c[a]); }  \
} while(0)


#endif
