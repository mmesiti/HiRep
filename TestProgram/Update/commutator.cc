#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace std;

#include "sun.h"

#ifdef _REPR_FUNDAMENTAL_
#include "fundamental.h"
#elif _REPR_ADJOINT_
#include "adjoint.h"
#elif _REPR_ANTISYMMETRIC_
#include "antisymmetric.h"
#elif _REPR_SYMMETRIC_
#include "symmetric.h"
#endif

#include "representation.h"

string commutator_over_i()
{

    ostringstream RET;
    RET << "void commutator(suNg_algebra_vector *c, suNg_algebra_vector*a, suNg_algebra_vector *b){\n";
    smatrix ifCAB[group::DIM];
    for(int C=0;C < group::DIM; C++) ifCAB[C].size = group::DIM;
    for(int A = 0; A < group::DIM; A++) for(int B=A+1; B < group::DIM; B++){
        smatrix COMM(group::N),BA(group::N);
        COMM.mult(group::T[A],group::T[B]);
        BA.mult(group::T[B],group::T[A]);
        BA.scale(-1.0);
        COMM.add(BA);
        for(int C=0;C < group::DIM; C++){
            smatrix tmp;
            complex coeff;
            tmp.mult(COMM,group::T[C]);
            trace(coeff,tmp);
            coeff *= 1.0/group::Tnorm;
            ifCAB[C].set(A,B,coeff);
            ifCAB[C].set(B,A,-coeff);
        }
    }
    for(int C=0;C < group::DIM; C++){
        RET << "\tc->c["<<C<<"] = "  ;
        for(int i=0;i<ifCAB[C].length;++i){
            int ai = ifCAB[C][i].index.row;
            int bi = ifCAB[C][i].index.col;
            double v = ifCAB[C][i].value.imag();
            if(v<0) RET << " -";
            else RET << " +"; 
            RET << fixed <<  setw(14) <<  setprecision(14) << fabs(v)<<"*"<<"a->c["
                << setw(2) << ai<<"]*b->c["
                << setw(2) << bi<<"]";
            if(i%2==0 && i != ifCAB[C].length-1) RET << "\n\t\t";
        }
        RET << ";\n";
    }
    RET << "}\n";

    return RET.str();
}

int main(int argc, char* argv[]) {
	if(argc != 2) {
		cerr << "Usage: "<<argv[0]<<" N\n";
		cerr << "N => number of colors\n";
		return 1;
	}
	
	int N = atoi(argv[1]);
	if(N < 2) {
		cerr << argv[0] <<": ERROR: number of color must be >1 !\n";
		return 2;
	}
	
	group::init(N);
	representation::init();



    string commutator_over_i_string = commutator_over_i();
    ofstream output("commutator.h");
    output << "/* This file has been generated automatically by " << argv[0] 
        << "*/"<<endl;
    output << "#ifndef COMMUTATOR_H" <<endl;
    output << "#define COMMUTATOR_H" <<endl;
    output << "#include \"suN.h\"" << endl;
    output << commutator_over_i_string ;
    output << "#endif" << endl;
    output.close();

	return 0;
}
