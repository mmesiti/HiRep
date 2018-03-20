#ifndef SUN_H
#define SUN_H
#include "./matrix.h"

string infinitesimal_evolution(const char* vname, const char* hname, const char* uname, const char* dtname);
string fundamental_algebra_represent(const char* mname, const char* hname);
string fundamental_algebra_project(const char* hname, const char* mname);


namespace group
{
    int N;
    int DIM;
    enum group_t{
        TYPESUN,
        TYPESON,
        TYPESPN
    };

    smatrix* T;
    string name;
    FLOATING Tnorm;

    void init(int n);
    void init(int n, group_t type, smatrix*& TOUT);
};

void group::init(int n){
#ifdef _GAUGE_SON_
	group::init(n,group::TYPESON,group::T);
#elif _GAUGE_SPN_
	group::init(n,group::TYPESPN,group::T);
#else 
	group::init(n,group::TYPESUN,group::T);
#endif
}
void group::init(int n, group::group_t type, smatrix*& TOUT)
{
    switch(type){
        case TYPESON:
            std::cerr << " Initializing group SO(" << n << ")..... ";
            break;
        case TYPESPN:
            std::cerr << " Initializing group USP(" << n << ")..... ";
            break;
        case TYPESUN:
            std::cerr << " Initializing group SU(" << n << ")..... ";
            break;
        default:
            break;
    }

    int A;
    int a, b;

    N = n;

    switch(type){
        case TYPESON:
            DIM = N*(N-1)/2;
            TOUT = new smatrix[DIM];

            A = 0;
            for(a = 0; a < N; a++) for(b = a+1; b < N; b++){
                TOUT[A].size = N;
                TOUT[A].set(a,b, complex(.0,1.));
                TOUT[A].set(b,a, complex(.0,-1.));
                A++;
            }
            Tnorm = 2.0;
            break;
        case TYPESPN:

            if( N%2 ==1 ) {
                cout << "\nMatrix size N must be even for gauge group SPN\n" ;
                exit(1);
            }

            DIM = (N*N+N)/2;
            TOUT = new smatrix[DIM];

            A = 0;
            for(a = 0; a < N/2; a++) {
                TOUT[A].size = N;
                TOUT[A].set(a,a,complex(sqrt(.5),.0));
                TOUT[A].set(a+N/2,a+N/2, -complex(sqrt(.5),.0));
                A++;
            }
            for(a = 0; a < N/2; a++) for(b = 0; b < N/2; b++)
                if(a > b)
                {
                    TOUT[A].size = N;
                    TOUT[A].set(a,b, complex(sqrt(.25),.0));
                    TOUT[A].set(b,a, complex(sqrt(.25),.0));
                    TOUT[A].set(a+N/2,b+N/2, -complex(sqrt(.25),.0));
                    TOUT[A].set(b+N/2,a+N/2, -complex(sqrt(.25),.0));
                    A++;
                }
                else if(a < b)
                {
                    TOUT[A].size = N;
                    TOUT[A].set(a,b, complex(.0,sqrt(.25)));
                    TOUT[A].set(b,a, complex(.0,-sqrt(.25)));
                    TOUT[A].set(a+N/2,b+N/2, complex(.0,sqrt(.25)));
                    TOUT[A].set(b+N/2,a+N/2, complex(.0,-sqrt(.25)));
                    A++;
                }

            for(a = 0; a < N/2; a++) {
                TOUT[A].size = N;
                TOUT[A].set(a,a+N/2, complex(sqrt(.5),.0));
                TOUT[A].set(a+N/2,a, complex(sqrt(.5),.0));
                A++;
                TOUT[A].size = N;
                TOUT[A].set(a,a+N/2, complex(.0,sqrt(.5)));
                TOUT[A].set(a+N/2,a, -complex(.0,sqrt(.5)));
                A++;
            }
            for(a = 0; a < N/2; a++) for(b = 0; b < N/2; b++)
                if(a > b)
                {
                    TOUT[A].size = N;
                    TOUT[A].set(a,b+N/2, complex(sqrt(.25),.0));
                    TOUT[A].set(b,a+N/2, complex(sqrt(.25),.0));
                    TOUT[A].set(a+N/2,b, complex(sqrt(.25),.0));
                    TOUT[A].set(b+N/2,a, complex(sqrt(.25),.0));
                    A++;
                }
                else if(a < b)
                {
                    TOUT[A].size = N;
                    TOUT[A].set(a,b+N/2, complex(.0,sqrt(.25)));
                    TOUT[A].set(b,a+N/2, complex(.0,sqrt(.25)));
                    TOUT[A].set(a+N/2,b, -complex(.0,sqrt(.25)));
                    TOUT[A].set(b+N/2,a, -complex(.0,sqrt(.25)));
                    A++;
                }
            Tnorm = 1.0;
            break;
        case TYPESUN:

            DIM = N*N-1;
            TOUT = new smatrix[DIM];

            A = 0;
            for(a = 0; a < N; a++) for(b = 0; b < N; b++)
                if(a > b)
                {
                    TOUT[A].size = N;
                    TOUT[A].set(a,b, complex(1.,.0));
                    TOUT[A].set(b,a, complex(1.,.0));
                    A++;
                }
                else if(a < b)
                {
                    TOUT[A].size = N;
                    TOUT[A].set(a,b, complex(.0,1.));
                    TOUT[A].set(b,a, complex(.0,-1.));
                    A++;
                }
                else if(a == b && a != 0)
                {
                    TOUT[A].size = N;
                    for(int k = 0; k < a; k++)
                        TOUT[A].set(k,k, complex(sqrt(2./(a*(a+1.))),.0));
                    TOUT[A].set(a,a, complex(-a*sqrt(2./(a*(a+1.))),.0));
                    A++;
                }
            Tnorm = 2.0;
            break;
        default:
            break;
    }

    //my changes below
    for (A=0;A<DIM;A++)
        TOUT[A].scale(sqrt(.5)/sqrt(Tnorm));
    Tnorm=0.5;
    
#ifndef NDEBUG 
    cerr << "OK\n";
#endif
}


string infinitesimal_evolution(const char* vname, const char* hname, const char* uname, const char* dtname)
{
    string RET;

    rvector H(group::DIM,hname);
    cmatrix U(group::N,uname);
    pmatrix M(group::N);
    pmatrix V(group::N);
    rvariable dt(dtname);

    dt.scale(complex(0.0,1.0));
    H.scale(dt);

    for(int A = 0; A < group::DIM; A++)
    {
        pmatrix T(group::T[A]);
        T.scale(H.get(A));
        M.add(T);
    }

    V.mult(M, U);

    RET = V.assignment("+=", vname);

    return RET;
}


string ExpX(const char* dtname,  const char* hname, const char* uname)
{
    ostringstream RET;

    rvector H(group::DIM,hname);
    pmatrix M(group::N);
    rvariable dt(dtname);

    dt.scale(complex(0.0,1.0));
    H.scale(dt);

    for(int A = 0; A < group::DIM; A++)
    {
        pmatrix T(group::T[A]);
        T.scale(H.get(A));
        M.add(T);
    }

    RET << 
        "\tdouble y[3];\n" << 
        "\tdouble s[" << group::N*(group::N-1)/2 << "][4];\n";

#ifdef _GAUGE_SON_
    RET <<
        "\tsuNgc ut, *u;\n\n"
        "\tfor (int i=0; i<NG*NG; ++i) { ut.c[i].re=r->c[i]; ut.c[i].im=0.; }\n"
        "\tfor (int i=0; i<NG*NG; ++i) { ut.c[i].re=r->c[i]; ut.c[i].im=0.; }\n"
        "\tu=&ut;\n\n";
#endif

#ifdef _GAUGE_SPN_
    RET <<
        "\tsuNgfull ut, *u;\n\n"
        "\tfor (int i=0; i<NG*NG/2; ++i) { ut.c[i].re=r->c[i].re; ut.c[i].im=r->c[i].im; }\n"
        "\tfor (int i=0; i<NG/2; ++i) {\n"
        "\t\tfor (int j=0; j<NG/2; ++j) {\n"
        "\t\t\tint ind = NG*i+j;\n"
        "\t\t\tut.c[NG*NG/2+NG/2+ind].re=r->c[ind].re; ut.c[NG*NG/2+NG/2+ind].im=-r->c[ind].im;\n"
        "\t\t\tut.c[NG*NG/2+ind].re=-r->c[ind+NG/2].re; ut.c[NG*NG/2+ind].im=r->c[ind+NG/2].im;\n"
        "\t\t}\n"
        "\t}\n"
        "\tu=&ut;\n\n";
#endif
	
    int k = 0;
    for(int j = 1; j < group::N; j++)
        for(int i = 0; i < j; i++)
        {
            polynomial tmp;
            pconstant ntmp(1.0/group::N);
            tmp = M.get(j,j);
            tmp.minus();
            tmp += M.get(i,i);
            tmp *= ntmp;
            RET <<
                "\ty[0] = " << M.get(i,j).str_imag() << ";\n" <<
                "\ty[1] = " << M.get(i,j).str_real() << ";\n" <<
                "\ty[2] = " << tmp.str_imag() << ";\n" <<
                "\tYtoU(s[" << k << "],y);\n";
            for(int p = 0; p < group::N; p++)
                RET << "\tsu2_rotate(s[" << k << "],&("
                    << uname << mindex(i,p,group::N) << "),&("
                    << uname << mindex(j,p,group::N) << "));\n";
            k++;
        }

    k = group::N*(group::N-1)/2 - 1;
    for(int j = group::N-1; j >= 1; j--)
        for(int i = j-1; i >= 0; i--)
        {
            for(int p = 0; p < group::N; p++)
                RET << "\tsu2_rotate(s[" << k << "],&("
                    << uname << mindex(i,p,group::N) << "),&("
                    << uname << mindex(j,p,group::N) << "));\n";
            k--;
        }
#ifdef _GAUGE_SON_
    RET<<"\n\tfor (int i=0; i<NG*NG; ++i) { r->c[i]=ut.c[i].re; }\n";
#endif
#ifdef _GAUGE_SPN_
    RET <<
        "\n\tfor (int i=0; i<NG*NG/2; ++i) { r->c[i].re=ut.c[i].re; r->c[i].im=ut.c[i].im; }\n";
#endif

    return RET.str();
}


string fundamental_algebra_represent(const char* mname, const char* hname)
{
	string RET;
	rvector H(group::DIM,hname);
	pmatrix M(group::N);
	pconstant I(complex(0.0,1.0));
	
	for(int A = 0; A < group::DIM; A++)
	{
		pmatrix iT(group::T[A]);
		iT.scale(H.get(A));
		iT.scale(I);
		M.add(iT);
	}
	
#ifndef _GAUGE_SPN_	
	RET = M.assignment("=", mname);
#else
 	RET = M.symplectic_compressed_assignment("=", mname);
#endif

	return RET;
}


string fundamental_algebra_project(const char* hname, const char* mname)
{
    string RET;
    pvector H(group::DIM);
    pmatrix *M;
    //	pmatrix adjM(group::N);
    pconstant I(complex(0.0,1.0));

#ifdef _GAUGE_SON_
    M = new rmatrix(group::N,mname);
#elif _GAUGE_SPN_
    M = new spmatrix(group::N,mname);
#else
    M = new cmatrix(group::N,mname);
#endif

    //	adjM = *M;
    //	adjM.adjoint();
    //	M->sub(adjM);

    for(int A = 0; A < group::DIM; A++)
    {
        pmatrix iT(group::T[A]);
        iT.scale(I);
        polynomial iTM;
        herm(iTM,iT,*M);
        iTM.real();
        iTM.scale(1.0/group::Tnorm);
        H.set(A, iTM);
    }

    RET = H.assignment("=", hname);


    delete M;

    return RET;
}
#endif
