//#ifndef SUN_REPR_FUNC_FUND_H
//#define SUN_REPR_FUNC_FUND_H

//#define _FUND_NORM2 (+5.000000000000000e-01)
//#define _REPR_NORM2 (+1.000000000000000e+00)



#define _fund_algebra_project_FMAT(h,m) \
        (h).c[0] = -(-(m).c[0].im+(m).c[10].im); \
        (h).c[1] = -(m).c[15].im+(m).c[5].im; \
        (h).c[2] = +7.071067811865476e-01*(-(m).c[11].re+(m).c[14].re-(m).c[1].re+(m).c[4].re); \
        (h).c[3] = +7.071067811865476e-01*(-(m).c[11].im-(m).c[14].im+(m).c[1].im+(m).c[4].im); \
        (h).c[4] = +(m).c[2].im+(m).c[8].im; \
        (h).c[5] = -(m).c[2].re+(m).c[8].re; \
        (h).c[6] = +(m).c[13].im+(m).c[7].im; \
        (h).c[7] = -(-(m).c[13].re+(m).c[7].re); \
        (h).c[8] = +7.071067811865476e-01*(+(m).c[12].re-(m).c[3].re-(m).c[6].re+(m).c[9].re); \
        (h).c[9] = +7.071067811865476e-01*(+(m).c[12].im+(m).c[3].im+(m).c[6].im+(m).c[9].im);

