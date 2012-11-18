/*******************************************************************************
*
* Computation of the observable E(t) evolved with the Wilson flow
*
*******************************************************************************/

#define MAIN_PROGRAM

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
#include "wilsonflow.h"

#include "cinfo.c"



#ifdef WITH_MPI
#error The Wilson Flow does not work with MPI
#endif

#if defined(BASIC_SF) || defined(ROTATED_SF)
#error The Wilson Flow does not work with SF
#endif





typedef struct _input_WF {
  double tmax;
  int nmeas;
  int nint;

  /* for the reading function */
  input_record_t read[4];

} input_WF;

#define init_input_WF(varname) \
{ \
  .read={\
    {"WF max integration time", "WF:tmax = %lf", DOUBLE_T, &((varname).tmax)},\
    {"WF number of measures", "WF:nmeas = %d", DOUBLE_T, &((varname).nmeas)},\
    {"WF number of integration steps between measures", "WF:nint = %d", INT_T, &((varname).nint)},\
    {NULL, NULL, 0, NULL}\
  }\
}


input_WF WF_var = init_input_WF(WF_var);

char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "wilsonflow.out";
enum { UNKNOWN_CNFG=0, DYNAMICAL_CNFG, QUENCHED_CNFG };


typedef struct {
  char string[256];
  int t, x, y, z;
  int nc, nf;
  double b, m;
  int n;
  int type;
} filename_t;


int parse_cnfg_filename(char* filename, filename_t* fn) {
  int hm;
  char *tmp = NULL;
  char *basename;

  basename = filename;
  while ((tmp = strchr(basename, '/')) != NULL) {
    basename = tmp+1;
  }            

/*#ifdef REPR_FUNDAMENTAL*/
/*#define repr_name "FUN"*/
/*#elif defined REPR_SYMMETRIC*/
/*#define repr_name "SYM"*/
/*#elif defined REPR_ANTISYMMETRIC*/
/*#define repr_name "ASY"*/
/*#elif defined REPR_ADJOINT*/
/*#define repr_name "ADJ"*/
/*#endif*/
  hm=sscanf(basename,"%*[^_]_%dx%dx%dx%d%*[Nn]c%dr%*[FSA]%*[UYSD]%*[NMYJ]%*[Nn]f%db%lfm%lfn%d",
      &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->nf),&(fn->b),&(fn->m),&(fn->n));
  if(hm==9) {
    fn->m=-fn->m; /* invert sign of mass */
    fn->type=DYNAMICAL_CNFG;
    return DYNAMICAL_CNFG;
  }
/*#undef repr_name*/

  double kappa;
  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%d%*[Nn]f%db%lfk%lfn%d",
      &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->nf),&(fn->b),&kappa,&(fn->n));
  if(hm==9) {
    fn->m = .5/kappa-4.;
    fn->type=DYNAMICAL_CNFG;
    return DYNAMICAL_CNFG;
  }

  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%db%lfn%d",
      &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->b),&(fn->n));
  if(hm==7) {
    fn->type=QUENCHED_CNFG;
    return QUENCHED_CNFG;
  }

  hm=sscanf(basename,"%*[^_]_%dx%dx%dx%d%*[Nn]c%db%lfn%d",
      &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->b),&(fn->n));
  if(hm==7) {
    fn->type=QUENCHED_CNFG;
    return QUENCHED_CNFG;
  }

  fn->type=UNKNOWN_CNFG;
  return UNKNOWN_CNFG;
}


void read_cmdline(int argc, char* argv[]) {
  int i, ai=0, ao=0, ac=0, al=0, am=0;
  FILE *list=NULL;

  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) ai=i+1;
    else if (strcmp(argv[i],"-o")==0) ao=i+1;
    else if (strcmp(argv[i],"-c")==0) ac=i+1;
    else if (strcmp(argv[i],"-l")==0) al=i+1;
    else if (strcmp(argv[i],"-m")==0) am=i;
  }

  if (am != 0) {
    print_compiling_info();
    exit(0);
  }

  if (ao!=0) strcpy(output_filename,argv[ao]);
  if (ai!=0) strcpy(input_filename,argv[ai]);

  error((ac==0 && al==0) || (ac!=0 && al!=0),1,"parse_cmdline [WF_measure.c]",
      "Syntax: mk_wilsonloops { -c <config file> | -l <list file> } [-i <input file>] [-o <output file>] [-m]");

  if(ac != 0) {
    strcpy(cnfg_filename,argv[ac]);
    strcpy(list_filename,"");
  } else if(al != 0) {
    strcpy(list_filename,argv[al]);
    error((list=fopen(list_filename,"r"))==NULL,1,"parse_cmdline [WF_measure.c]" ,
	"Failed to open list file\n");
    error(fscanf(list,"%s",cnfg_filename)==0,1,"parse_cmdline [WF_measure.c]" ,
	"Empty list file\n");
    fclose(list);
  }


}


int main(int argc,char *argv[]) {
  int i;
  char tmp[256];
  FILE* list;
  filename_t fpars;

  /* setup process id and communications */
  read_cmdline(argc, argv);
  setup_process(&argc,&argv);

  /* logger setup */
  /* disable logger for MPI processes != 0 */
  logger_setlevel(0,70);
  if (PID!=0) { logger_disable(); }
  if (PID==0) { 
    sprintf(tmp,">%s",output_filename); logger_stdout(tmp);
    sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
  }

  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS); 
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  lprintf("MAIN",0,"input file [%s]\n",input_filename); 
  lprintf("MAIN",0,"output file [%s]\n",output_filename); 
  if (strcmp(list_filename,"")!=0) lprintf("MAIN",0,"list file [%s]\n",list_filename); 
  else lprintf("MAIN",0,"cnfg file [%s]\n",cnfg_filename); 


  /* read & broadcast parameters */
  parse_cnfg_filename(cnfg_filename,&fpars);

  read_input(glb_var.read,input_filename);
  read_input(WF_var.read,input_filename);
  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;
 
  error(fpars.type==UNKNOWN_CNFG,1,"WF_measure.c","Bad name for a configuration file");
  error(fpars.nc!=NG,1,"WF_measure.c","Bad NG");

  lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_var.rlxd_level,glb_var.rlxd_seed);
  rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed+PID);
  srand(glb_var.rlxd_seed+PID);

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  init_BCs(NULL);
  
  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);

  
  lprintf("MAIN",0,"WF tmax: %e\n",WF_var.tmax);
  lprintf("MAIN",0,"WF number of measures: %d\n",WF_var.nmeas);
  lprintf("MAIN",0,"WF time lapse between measures: %e\n",WF_var.tmax/WF_var.nmeas);
  lprintf("MAIN",0,"WF number of integration intervals per measure: %d\n",WF_var.nint);
  lprintf("MAIN",0,"WF number of integration intervals: %d\n",WF_var.nint*WF_var.nmeas);
  lprintf("MAIN",0,"WF integration step: %e\n",WF_var.tmax/(WF_var.nmeas*WF_var.nint));

  
  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_mesons.c]" ,
	"Failed to open list file\n");
  }

  WF_initialize();
  
  i=0;
  while(1) {

    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    i++;

    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    /* NESSUN CHECK SULLA CONSISTENZA CON I PARAMETRI DEFINITI !!! */
    read_gauge_field(cnfg_filename);

    lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());

    full_plaquette();

    int k, n;
    double epsilon=WF_var.tmax/(WF_var.nmeas*WF_var.nint);
    double t=0.;

    double E=WF_E(u_gauge);
    double Esym=WF_Esym(u_gauge);
    lprintf("WILSONFLOW",0,"WF (ncnfg,t,E,t2*E,Esym,t2*Esym) = %d %e %e %e %e %e\n",i,t,E,t*t*E,Esym,t*t*Esym);
    for(n=0;n<WF_var.nmeas;n++) {
      for(k=0;k<WF_var.nint;k++) {
        WilsonFlow3(u_gauge,epsilon);
        t+=epsilon;
      }
      E=WF_E(u_gauge);
      Esym=WF_Esym(u_gauge);
      lprintf("WILSONFLOW",0,"WF (ncnfg,t,E,t2*E,Esym,t2*Esym) = %d %e %e %e %e %e\n",i,t,E,t*t*E,Esym,t*t*Esym);      
    }
    
    if(list==NULL) break;
  }

  if(list!=NULL) fclose(list);

  WF_free();
  
  free_BCs();
 
  free_gfield(u_gauge);

  finalize_process();
  
  return 0;
}

