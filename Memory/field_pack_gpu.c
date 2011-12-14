/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
 *
 * File field_pack_gpu.c
 *
 * Functions for fields allocation on GPUs
 *
 *******************************************************************************/

#include <stdlib.h>
#include "suN_types.h"
#include "error.h"
#include "memory.h"
#include "global.h"
#include "spinor_field.h"
#include "geometry.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif
#include "gpu.h"


void spinor_field_togpuformat(spinor_field *out, spinor_field *in) {
    _DECLARE_INT_ITERATOR(ix);
    suNf_spinor *r=0;

    //check input and output type are the same
    error(out->type!=in->type,1,"spinor_field_togpuformat " __FILE__, "Spinors don't match!");
    
#ifdef UPDATE_EO
    if (in->type==&glattice) {
        // we call recursively this function twice
        // on the even and odd sublattices
        in->type=out->type=&glat_even;
        spinor_field_togpuformat(out, in);
    	in->type=out->type=&glat_odd;
        spinor_field_togpuformat(out, in);
        in->type=out->type=&glattice;
        return;
    }
#endif //UPDATE_EO
    
    _PIECE_FOR(in->type,ix) {
        const int start = in->type->master_start[_PIECE_INDEX(ix)];
        const int N = in->type->master_end[_PIECE_INDEX(ix)]-in->type->master_start[_PIECE_INDEX(ix)]+1;
        complex *cout=(complex*)(_FIELD_AT(out,start));
        _SITE_FOR(in->type,ix) {
        
        	r=_FIELD_AT(in,ix);
            
            for (int j=0; j<sizeof(*r)/sizeof(complex); ++j) {
            	cout[j*N]=((complex*)(r))[j];
            }
            ++cout;
       
    	}
    }
}

void spinor_field_tocpuformat(spinor_field *out, spinor_field *in) {
    _DECLARE_INT_ITERATOR(ix);
    suNf_spinor *r=0;

    //check input and output type are the same
    error(out->type!=in->type,1,"spinor_field_tocpuformat " __FILE__, "Spinors don't match!");
    
#ifdef UPDATE_EO
    if (in->type==&glattice) {
        // we call recursively this function twice
        // on the even and odd sublattices
        in->type=out->type=&glat_even;
        spinor_field_togpuformat(out, in);
    	in->type=out->type=&glat_odd;
        spinor_field_togpuformat(out, in);
        in->type=out->type=&glattice;
        return;
    }
#endif //UPDATE_EO    
    
    _PIECE_FOR(in->type,ix) {
        int start = in->type->master_start[_PIECE_INDEX(ix)];
        int N = in->type->master_end[_PIECE_INDEX(ix)]-in->type->master_start[_PIECE_INDEX(ix)]+1; 
        complex *cin=(complex*)(_FIELD_AT(in,start));
        _SITE_FOR(in->type,ix) {
            
        	r=_FIELD_AT(out,ix);
        	
            for (int j=0; j<sizeof(*r)/sizeof(complex); ++j) {
                ((complex*)(r))[j]=cin[j*N];
            }
            ++cin;
            
    	}
    }
}

void spinor_field_togpuformat_flt(spinor_field_flt *out, spinor_field_flt *in) {
  _DECLARE_INT_ITERATOR(ix);
  suNf_spinor_flt *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"spinor_field_togpuformat_flt " __FILE__, "Spinors don't match!");
  
#ifdef UPDATE_EO
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    spinor_field_togpuformat_flt(out, in);
    in->type=out->type=&glat_odd;
    spinor_field_togpuformat_flt(out, in);
    in->type=out->type=&glattice;
    return;
  }
#endif //UPDATE_EO
  
  _PIECE_FOR(in->type,ix) {
    const int start = in->type->master_start[_PIECE_INDEX(ix)];
    const int N = in->type->master_end[_PIECE_INDEX(ix)]-in->type->master_start[_PIECE_INDEX(ix)]+1;
    complex_flt *cout=(complex_flt*)(_FIELD_AT(out,start));
    _SITE_FOR(in->type,ix) {
      
      r=_FIELD_AT(in,ix);
      
      for (int j=0; j<sizeof(*r)/sizeof(complex_flt); ++j) {
        cout[j*N]=((complex_flt*)(r))[j];
      }
      ++cout;
      
    }
  }
}

void spinor_field_tocpuformat_flt(spinor_field_flt *out, spinor_field_flt *in) {
  _DECLARE_INT_ITERATOR(ix);
  suNf_spinor_flt *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"spinor_field_tocpuformat_flt " __FILE__, "Spinors don't match!");
  
#ifdef UPDATE_EO
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    spinor_field_tocpuformat_flt(out, in);
    in->type=out->type=&glat_odd;
    spinor_field_tocpuformat_flt(out, in);
    in->type=out->type=&glattice;
    return;
  }
#endif //UPDATE_EO    
  
  _PIECE_FOR(in->type,ix) {
    int start = in->type->master_start[_PIECE_INDEX(ix)];
    int N = in->type->master_end[_PIECE_INDEX(ix)]-in->type->master_start[_PIECE_INDEX(ix)]+1; 
    complex_flt *cin=(complex_flt*)(_FIELD_AT(in,start));
    _SITE_FOR(in->type,ix) {
      
      r=_FIELD_AT(out,ix);
      
      for (int j=0; j<sizeof(*r)/sizeof(complex_flt); ++j) {
        ((complex_flt*)(r))[j]=cin[j*N];
      }
      ++cin;
      
    }
  }
}


void gfield_togpuformat(suNg_field *out, suNg_field *in) {
  _DECLARE_INT_ITERATOR(ix);
  suNg *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"gield_togpuformat " __FILE__, "Gauge field types don't match!");
  
#ifdef UPDATE_EO
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    gfield_togpuformat(out, in);
    in->type=out->type=&glat_odd;
    gfield_togpuformat(out, in);
    in->type=out->type=&glattice;
    return;
  }
#endif //UPDATE_EO
  
  _PIECE_FOR(in->type,ix) {
    const int start = in->type->master_start[_PIECE_INDEX(ix)];
    const int N = in->type->master_end[_PIECE_INDEX(ix)]-in->type->master_start[_PIECE_INDEX(ix)]+1;
    double *cout=(double*)(_4FIELD_AT(out,start,0));
    _SITE_FOR(in->type,ix) {
      
      r=_4FIELD_AT(in,ix,0);
      
      for (int j=0; j<4*sizeof(*r)/sizeof(double); ++j) {
        cout[j*N]=((double*)(r))[j];
      }
      ++cout;
    }
  }
}

void gfield_tocpuformat(suNg_field *out, suNg_field *in) {
  _DECLARE_INT_ITERATOR(ix);
  suNg *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"gield_tocpuformat " __FILE__, "Gauge field types don't match!");
  
#ifdef UPDATE_EO
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    gfield_togpuformat(out, in);
    in->type=out->type=&glat_odd;
    gfield_togpuformat(out, in);
    in->type=out->type=&glattice;
    return;
  }
#endif //UPDATE_EO
  
  _PIECE_FOR(in->type,ix) {
    const int start = in->type->master_start[_PIECE_INDEX(ix)];
    const int N = in->type->master_end[_PIECE_INDEX(ix)]-in->type->master_start[_PIECE_INDEX(ix)]+1;
    double *cin=(double*)(_4FIELD_AT(out,start,0));
    _SITE_FOR(in->type,ix) {
      
      r=_4FIELD_AT(in,ix,0);
      
      for (int j=0; j<4*sizeof(*r)/sizeof(double); ++j) {
        ((double*)(r))[j]=cin[j*N];
      }
      ++cin;
      
    }
  }
}

void gfield_togpuformat_f(suNf_field *out, suNf_field *in) {
  gfield_togpuformat((suNg_field *)out, (suNg_field *)in);
}

void gfield_tocpuformat_f(suNf_field *out, suNf_field *in) {
  gfield_tocpuformat((suNg_field *)out, (suNg_field *)in);
}

void gfield_togpuformat_flt(suNg_field_flt *out, suNg_field_flt *in) {
  _DECLARE_INT_ITERATOR(ix);
  suNg_flt *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"gield_togpuformat " __FILE__, "Gauge field types don't match!");
  
#ifdef UPDATE_EO
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    gfield_togpuformat_flt(out, in);
    in->type=out->type=&glat_odd;
    gfield_togpuformat_flt(out, in);
    in->type=out->type=&glattice;
    return;
  }
#endif //UPDATE_EO
  
  _PIECE_FOR(in->type,ix) {
    const int start = in->type->master_start[_PIECE_INDEX(ix)];
    const int N = in->type->master_end[_PIECE_INDEX(ix)]-in->type->master_start[_PIECE_INDEX(ix)]+1;
    float *cout=(float*)(_4FIELD_AT(out,start,0));
    _SITE_FOR(in->type,ix) {
      
      r=_4FIELD_AT(in,ix,0);
      
      for (int j=0; j<4*sizeof(*r)/sizeof(float); ++j) {
        cout[j*N]=((float*)(r))[j];
      }
      ++cout;
      
    }
  }
}

void gfield_tocpuformat_flt(suNg_field_flt *out, suNg_field_flt *in) {
  _DECLARE_INT_ITERATOR(ix);
  suNg_flt *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"gield_tocpuformat " __FILE__, "Gauge field types don't match!");
  
#ifdef UPDATE_EO
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    gfield_togpuformat_flt(out, in);
    in->type=out->type=&glat_odd;
    gfield_togpuformat_flt(out, in);
    in->type=out->type=&glattice;
    return;
  }
#endif //UPDATE_EO
  
  _PIECE_FOR(in->type,ix) {
    const int start = in->type->master_start[_PIECE_INDEX(ix)];
    const int N = in->type->master_end[_PIECE_INDEX(ix)]-in->type->master_start[_PIECE_INDEX(ix)]+1;
    float *cin=(float*)(_4FIELD_AT(out,start,0));
    _SITE_FOR(in->type,ix) {
      
      r=_4FIELD_AT(in,ix,0);
      
      for (int j=0; j<4*sizeof(*r)/sizeof(double); ++j) {
        ((float*)(r))[j]=cin[j*N];
      }
      ++cin;
      
    }
  }
}

void gfield_togpuformat_f_flt(suNf_field_flt *out, suNf_field_flt *in) {
  gfield_togpuformat_flt((suNg_field_flt *)out, (suNg_field_flt *)in);
}

void gfield_tocpuformat_f_flt(suNf_field_flt *out, suNf_field_flt *in) {
  gfield_tocpuformat_flt((suNg_field_flt *)out, (suNg_field_flt *)in);
}


