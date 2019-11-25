#include "suN.h"
#include "suN_repr_func.h"
#include "utils.h"
#include "stdio.h"
#include "stdlib.h"
#include "random.h"
#include "./commutator.h"

#ifdef GAUGE_SPN
 #define ALGDIM NG*(NG+1)/2
#elif defined(GAUGE_SUN)
 #define ALGDIM NG*NG-1
#else
 #error "Group not valid"
#endif

#ifdef GAUGE_SUN
double group_deviation(suNg sunmat){
    suNg tmp1,identity;
    double should_be_zero;
    _suNg_unit(identity);
    _suNg_times_suNg_dagger(tmp1,sunmat,sunmat);
    _suNg_sub_assign(tmp1,identity);
    _suNg_sqnorm(should_be_zero,tmp1);
    return should_be_zero;
}
#elif defined(GAUGE_SPN)
double group_deviation(suNg spnmat){
    suNg Omega,tmp1,tmp2,spntrans;
    double should_be_zero;
    _symplectic(Omega);
    _suNg_times_suNg(tmp1,Omega,spnmat);
    _suNg_transpose(spntrans,spnmat);
    _suNg_times_suNg(tmp2,spntrans,tmp1);

    _suNg_sub_assign(tmp2,Omega);

    _suNg_sqnorm(should_be_zero,tmp2);
    return should_be_zero;
}
#endif

void print_suNg(suNg *m){
    int i,j;
    double sumsq = 0;
    double traceRe = 0;

    printf("\n");
#ifdef GAUGE_SPN
    for(i=0;i<NG/2;++i){
#else
    for(i=0;i<NG;++i){
#endif
        for(j=0;j<NG;++j){
            printf("%+-.6f,%+-.6f ", m->c[i*NG+j].re, m->c[i*NG+j].im); 
            //printf("%+-.3f ", m->c[i*NF+j].im); 
            sumsq += m->c[i*NG+j].im*m->c[i*NG+j].im;
            sumsq += m->c[i*NG+j].re*m->c[i*NG+j].re;
            if(i==j) traceRe +=m->c[i*NG+j].re;
        }
        printf("\n\n\n");
# ifdef GAUGE_SPN
    }
#else
    }
#endif
    printf("\n");
    printf("SUM SQUARES: %f\n", sumsq);
    printf("TRACE RE: %f\n", traceRe);
}

int main(){

    suNg_algebra_vector A,B; 

    gauss(A.c,ALGDIM);
    gauss(B.c,ALGDIM);

    printf("A algebra vector:\n");
    for(int i = 0; i<ALGDIM; ++i) printf("%f ",A.c[i]); 
    printf("\n");
    printf("B algebra vector:\n");
    for(int i = 0; i<ALGDIM; ++i) printf("%f ",B.c[i]); 
    printf("\n");

    //A+A
    suNg_algebra_vector ApA;
    _algebra_vector_zero_g(ApA);
    _algebra_vector_add_assign_g(ApA,A);
    _algebra_vector_add_assign_g(ApA,A);
    //A+B
    suNg_algebra_vector ApB;
    _algebra_vector_zero_g(ApB);
    _algebra_vector_add_assign_g(ApB,A);
    _algebra_vector_add_assign_g(ApB,B);

    printf("-Test that exp(A) != 0 ...");
    {
        suNg expA;
        _suNg_unit(expA);
        ExpX(0.0, &A, &expA);
        double sq = 0.0;
        _suNg_sqnorm(sq,expA);
        if(sq < 0.0001){
            print_suNg(&expA);
            return 1;
        }
    }
    printf("OK\n");


    printf("-Test that exp(A) != 1 for A != 0 ...");
    {
        suNg expA;
        _suNg_unit(expA);
        ExpX(0.5, &A, &expA);
        suNg tmp;
        _suNg_unit(tmp);
        _suNg_sub_assign(tmp,expA);
        double diffsq = 0.0;
        _suNg_sqnorm(diffsq,tmp);
        if(diffsq < 0.0001){
            print_suNg(&expA);
            return 1;
        }
    }
    printf("OK\n");

    printf("-Test if exp((l+k)*a) == epx(la)exp(ka)...");
    {
        double diffsq;
        const double l = -0.1;
        const double k = 0.21;
        suNg explA;
        _suNg_unit(explA);
        ExpX(l, &A, &explA);
        suNg expkA;
        _suNg_unit(expkA);
        ExpX(k, &A, &expkA);
        suNg expkplA;
        _suNg_unit(expkplA);
        ExpX(k+l, &A, &expkplA);
 
        suNg diff;
        _suNg_times_suNg(diff,explA,expkA);
        _suNg_sub_assign(diff,expkplA);
        _suNg_sqnorm(diffsq,diff);
        if(diffsq > 1.0e-13)
        {
            printf("\nDiffsq: %f\n",diffsq);
            print_suNg(&explA);
            print_suNg(&expkA);
            for(int i = 0; i<ALGDIM; ++i) printf("%f ",A.c[i]); 
            printf("\n");
            suNg tmp;
            _suNg_times_suNg(tmp,explA,expkA);
            print_suNg(&tmp);
            print_suNg(&expkplA);
            return 1;
        }
    }
    printf("OK\n");

    printf("Test if exp(epxA^\\dagger B expA) = expA^\\dagger expB epxA...");
    {
        const double coeff1 = 0.35;
        const double coeff2 = -0.5;
        suNg Bmat;
        _algebra_represent(Bmat,B);
        suNg expA;
        _suNg_unit(expA);
        ExpX(coeff1, &A, &expA);
        suNg expAdaggerB;
        _suNg_dagger_times_suNg(expAdaggerB,expA,Bmat);
        suNg Brotmat;
        _suNg_times_suNg(Brotmat,expAdaggerB,expA);//epxA^\dagger B expA
        suNg_algebra_vector Brotated;
#ifdef GAUGE_SPN
        suNgfull Brotmat_full;
        _suNg_expand(Brotmat_full,Brotmat);
#else   
        suNg Brotmat_full = Brotmat;
#endif
        _algebra_project(Brotated,Brotmat_full);
        suNg expBrot;
        _suNg_unit(expBrot);
        ExpX(coeff2,&Brotated,&expBrot);
        suNg expB;
        _suNg_unit(expB);
        ExpX(coeff2, &B,&expB);
        suNg expAdaggerexpB;
        _suNg_dagger_times_suNg(expAdaggerexpB,expA,expB);
        suNg expBrot2;
        _suNg_times_suNg(expBrot2,expAdaggerexpB,expA);
        suNg diff = expBrot;
        _suNg_sub_assign(diff,expBrot2);
        double diffsq;
        _suNg_sqnorm(diffsq,diff);
        if(diffsq > 1.0e-13)
        {
            printf("\nDiffsq: %f\n",diffsq);
            print_suNg(&expBrot);
            print_suNg(&expBrot2);
            return 1;
        }
    }
    printf("OK\n");

    //printf("-test if Campbell-Baker-Hausdorff formula holds...");
    //{
    //    // to implement
    //}



    return 0;

}


