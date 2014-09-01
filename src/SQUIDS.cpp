 /******************************************************************************
 *    This program is free software: you can redistribute it and/or modify     *
 *   it under the terms of the GNU General Public License as published by      *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   This program is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                             *   
 *   Authors:                                                                  *
 *      Carlos Arguelles (University of Wisconsin Madison)                     * 
 *         carguelles@icecube.wisc.edu                                         *
 *      Jordi Salvado (University of Wisconsin Madison)                        *
 *         jsalvado@icecube.wisc.edu                                           *
 ******************************************************************************/


#include "SQUIDS.h"

//#define CalNeuOscSUN_DEBUG

// Macros
#define SQR(x)      ((x)*(x))                        // x^2
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  // |x|^2
#define POW10(x)    (exp(M_LN10*(x)))                // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( (b) > 0.0 ? (fabs(a)) : (-fabs(a)) )
#define KRONECKER(i,j)  ( (i)==(j) ? 1 : 0 )


SQUIDS::SQUIDS(void){
  is_init=false;
}

SQUIDS::SQUIDS(int n,int ns,int nrh,int nsc){
  ini(n,ns,nrh,nsc);
}

void SQUIDS::ini(int n,int nsu,int nrh,int nsc){
  CoherentInt=false;
  NonCoherentInt=false;
  OtherInt=false;
  ScalarsInt=false;
  neu_and_aneu = false;
  AnyNumerics=(CoherentInt||NonCoherentInt||OtherInt||ScalarsInt);

  /*
    Setting the number of energy bins, number of components for the density matrix and
    number of scalar functions
  */

  nx=n;
  nsun=nsu;
  nrhos=nrh;
  nscalars=nsc;
  size_rho=SQR(nsun);
  size_state=(SQR(nsun)*nrhos+nscalars);
  /*
    Setting time units and initial time
  */
  tunit=1;
  t=0;

  //Allocate memeroy for the system
  Nsystem = nx*size_state;
  numeqn=Nsystem;
  system=new double[Nsystem];

  /*
    Initializing the SU algebra object, needed to compute algebraic operations like commutators,
    anticommutators, rotations and evolutions
  */

  SU.init(nsun);


  //Allocate memory
  x=new double[nx];
  delx=new double[nx];

  state=new SU_state[nx];
  dstate=new SU_state[nx];

  for(int ei = 0; ei < nx; ei++){
    state[ei].rho=new SU_vector[nrhos];
    dstate[ei].rho=new SU_vector[nrhos];
  }

  for(int ei = 0; ei < nx; ei++){
    for(int i=0;i<nrhos;i++){
      state[ei].rho[i].InitSU_vector(nsun,&(system[ei*size_state+i*size_rho]));
    }
    if(nscalars>0){
      state[ei].scalar=&(system[ei*size_state+nrhos*size_rho]);
    }
  }


  //default errors
  rel_error=1e-20;
  abs_error=1e-20;

  h        = tunit*1e-16;
  h_min    = tunit*1e-40;
  h_max    = tunit*5.0;

  // setting up GSL ODE solver
  s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, numeqn);
  c = gsl_odeiv2_control_y_new(abs_error,rel_error);
  e = gsl_odeiv2_evolve_alloc(numeqn);


  sys.function = &RHS;
  sys.jacobian = NULL;
  sys.dimension = (size_t)numeqn;
  sys.params = this;

  //  gsl_odeiv2_step_set_driver(s,d);


  is_init=true;
};

void SQUIDS::set_deriv_system_pointer(double *p){
  deriv_system=p;
  for(int ei = 0; ei < nx; ei++){
    for(int i=0;i<nrhos;i++){
      dstate[ei].rho[i].InitSU_vector(nsun,&(p[ei*size_state+i*size_rho]));
    }
    dstate[ei].scalar=&(p[ei*size_state+nrhos*size_rho]);
  }
}

void SQUIDS::set_system_pointer(double *p){
  system=p;
  for(int ei = 0; ei < nx; ei++){
    for(int i=0;i<nrhos;i++){
      state[ei].rho[i].InitSU_vector(nsun,&(p[ei*size_state+i*size_rho]));
    }
    state[ei].scalar=&(p[ei*size_state+nrhos*size_rho]);
  }
}



SQUIDS::~SQUIDS(void){
  if(is_init)
    free();
}


void SQUIDS::free(void){
  if(is_init){
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);

    delete x;
    delete delx;



    delete system;

  }else{
    cerr << "free call error, is not initialized" << endl;
    exit(1);
  }
}

/*
  This functions set's the grid of points in where we have a
  density matrix and a set of scalars
 */

int SQUIDS::Set_xrange(double xi, double xf, string type){
  if (xi == xf){
    x[0] = xi;
    return 0;
  }

  if(type=="linear" || type=="Linear" || type=="lin" || type=="Lin"){
    for(int e1 = 0; e1 < nx; e1++){
      x[e1]=xi+(xf-xi)*(double)e1/(double)(nx-1);
    }
  }else if(type=="log" || type=="Log"){
    double xmin_log,xmax_log;
    if (xi < 1.0e-10 ){
      cerr << "SQUIDS::Set_xrange : Xmin too small for log scale" << endl;
      exit(1);
    }else{
        xmin_log = log(xi);
    }
    xmax_log = log(xf);

    double step_log = (xmax_log - xmin_log)/double(nx);
    int i=0;
    for(int e1 = 0; e1 < nx; e1++){
      double X=xmin_log+(xmax_log-xmin_log)*(double)e1/(double)(nx-1);
      x[e1]=exp(X);
    }
  }else{
    cerr << "SQUIDS::Set_xrange : Not well deffined X range" << endl;
    exit(1);
  }

  for(int e1 = 1; e1 < nx; e1++){
    delx[e1] = x[e1] - x[e1-1];
  }
  return 0;
}


double SQUIDS::GetExpectationValue(SU_vector & op,  int nrh, int i){
  SU_vector h0=H0(x[i]);
  return state[i].rho[nrh]*op.SUEvolve(h0,t-t_ini);
}

double SQUIDS::GetExpectationValueD(SU_vector & op, int nrh,  double xi){
  SU_vector h0=H0(xi);
  int xid;
  for(unsigned int i = 0; i < nx; i++){
    if ( xi >= x[i] && xi <= x[i+1]){
      xid = i;
      break;
    }else{
      if(i==nx-1){
	cerr << "SQUIDS::GetExpectationValueD : x value not in the array.";
      }
    }
  }

  return (state[xid].rho[nrh] + 
	   (state[xid+1].rho[nrh]-
	    state[xid].rho[nrh])*((xi-x[xid])/(x[xid+1]-x[xid])))*op.SUEvolve(h0,t-t_ini);
}

int SQUIDS::Get_i(double xi){
  double xl, xr;
  int nr=nx-1;
  int nl=0;

  xl=x[nl];
  xr=x[nr];

  if(xi>xr || xi<xl){
    cout << xr << "  " << xl << endl;
    cerr << " Error IDS::Get_i :  value " << xi  <<" out of bounds" << endl;
    exit(1);
  }

  while((nr-nl)>1){
    if(((nr-nl)%2)!=0){
      if(nr<nx-1)nr++;
      else if(nl>0)nl--;
    }
    if(xi<(xl+(xr-xl)/2)){
      nr=nl+(nr-nl)/2;
      xr=x[nr];
    }else{
      nl=nl+(nr-nl)/2;
      xl=x[nl];
    }
  }
  return nl;
}

void SQUIDS::Set(string name,const gsl_odeiv2_step_type * opt){
  if(name=="GSL_Step" || name=="GSL_step" || name=="step"){ 
    gsl_odeiv2_step_free(s);
    s = gsl_odeiv2_step_alloc(opt, numeqn);  
  }
}

void SQUIDS::Set(string name,bool opt){
  if(name=="CoherentInteractions" || name=="coherentinteractoins" || name=="CohInt" || name=="H1"){
    CoherentInt=opt;
    AnyNumerics=(CoherentInt||NonCoherentInt||OtherInt||ScalarsInt);
  }else if(name=="NonCoherentInteractions" || name=="noncoherentinteractoins" || name=="NonCohInt"){
    NonCoherentInt=opt;
    AnyNumerics=(CoherentInt||NonCoherentInt||OtherInt||ScalarsInt);
  }else if(name=="OtherInt" || name=="otherint" || name=="OInt"){
    OtherInt=opt;
    AnyNumerics=(CoherentInt||NonCoherentInt||OtherInt||ScalarsInt);
  }else if(name=="ScalarInteractions" || name=="scalarinteractions" || name=="ScalInt"){
    ScalarsInt=opt;
    AnyNumerics=(CoherentInt||NonCoherentInt||OtherInt||ScalarsInt);
  }else if(name=="AntiNeutrinos" || name=="antinu" || name=="AntiNu" || name=="anu"){
    neu_and_aneu=opt;
  }else{
    cerr << "DMKS::Set : Option (" << name <<  ") not found. " << endl;
    exit(1);
  }
}


void SQUIDS::Set(string name,double opt){
  if(!params.Set(name,opt)){
    if(name=="h_min"){
      h_min=opt;
    }else if(name=="h_max"){
      h_max=opt;
    }else if(name=="h"){
      h=opt;
    }else if(name=="rel_error"){
      rel_error=opt;
      gsl_odeiv2_control_free(c);
      c = gsl_odeiv2_control_y_new(abs_error,rel_error);
    }else if(name=="abs_error"){
      abs_error=opt;
      gsl_odeiv2_control_free(c);
      c = gsl_odeiv2_control_y_new(abs_error,rel_error);
    }else  if(name=="t" || name=="T"){
      t=opt;
      t_ini=opt;
    }else if(name=="Units" || name=="units" || name=="tunits"|| name=="Tunits"){
      tunit=opt;
    }else{
      cerr << "DMKS::Set : Option (" << name <<  ") not found. " << endl;
      exit(1);
    }
  }
}


void SQUIDS::Set(string name,int opt){
  if(name=="nx"){
    if(opt!=nx){
      this->free();
      ini(opt,nsun,nrhos,nscalars);
    }
  }else if(name=="nsun"){
    if(opt!=nsun){
      this->free();
      ini(nx,opt,nrhos,nscalars);
    }
  }else if(name=="nrhos"){
    if(opt!=nrhos){
      this->free();
      ini(nx,nsun,opt,nscalars);
    }
  }else if(name=="nscalars"){
    if(opt!=nscalars){
      this->free();
      ini(nx,nsun,nrhos,opt);
    }
  }else{
    cerr << "DMKS::Set : Option (" << name <<  ") not found. " << endl;
    exit(1);
  }
}




int SQUIDS::Derive(double at){
  t=at;

  EvolveProjectors(at);
  for(int ei = 0; ei < nx; ei++){
    for(int i = 0; i < nrhos; i++){
      index_rho=i;
      // Density matrix
      // Coherent interaction
      if(CoherentInt){
  	//dstate[ei].rho[i] = state[ei].rho[i]*(-1);
	//dstate[ei].rho[i] = SU.iCommutator(HI(ei,t),state[ei].rho[i]);
  dstate[ei].rho[i] = SU.iCommutator(state[ei].rho[i],HI(ei,t));
  //double rho2 = state[ei].rho[i]*state[ei].rho[i];
  //std::cout << t/params.km << " " << rho2 << std::endl;
  //if( rho2 > 100.0 ) exit(1);
  /*
  if( rho2 > 1.001 and rho2 < 2 ) {
    std::cout << dstate[ei].rho[i] << std::endl;
    std::cout << HI(ei,t) << std::endl;
  }
  */
	//cout << ei << " " << i << " drho  " << dstate[ei].rho[i] << endl;
      }
      // Non coherent interaction
      if(NonCoherentInt)
  	dstate[ei].rho[i] -= SU.ACommutator(GammaRho(ei,t),state[ei].rho[i]);
      // Other possible interaction, for example involving the Scalars or non linear terms in rho.
      if(OtherInt)
  	dstate[ei].rho[i] += InteractionsRho(ei,t);
      //Scalars
    }
    if(ScalarsInt){
      for(int is=0;is<nscalars;is++){
  	index_scalar=is;
  	dstate[ei].scalar[is] = -state[ei].scalar[is]*GammaScalar(ei,t);
  	dstate[ei].scalar[is] += InteractionsScalar(ei,t);
      }
    }
  }

  return GSL_SUCCESS;
}



int SQUIDS::EvolveSUN(double ti, double tf){
  t_ini=ti;
  t_end=tf;
  if(AnyNumerics){
    int gsl_status = GSL_SUCCESS;
    t=ti;

    // defining ODE extra variables

    double x;                       // ODE independent variable
    double x_aux;                   // ODE2 independent variable
    double x_ini = ti;  // initial position
    double x_end = tf;  // final position
    // step sizes


    //    gsl_odeiv2_driver_set_hmax(d,h_max);


#ifdef CalNeuOscSUN_DEBUG
    printf("GSL paramers :\n");

    printf("x_ini : %lf \n", x_ini/tunit);
    printf("x_end : %lf \n", x_end/tunit);
    printf("h : %g \n", h/tunit);
    printf("h_min : %g \n", h_min/tunit);
#endif


    // initial position
    x = x_ini;
    //x_aux = x_ini;

    // ODE system error control
    c = gsl_odeiv2_control_y_new(abs_error,rel_error);

#ifdef CalNeuOscSUN_DEBUG
    int count = 0;
    int count_step = 100000;
#endif

#ifdef CalNeuOscSUN_DEBUG
    printf("Start calculation.\n");
#endif

    double * gsl_sys = system;

    while (x < x_end){
      t=x;
      double x_inter = x + (x_end-x_ini)/1000.0;

      gsl_status = gsl_odeiv2_evolve_apply(e,c,s,&sys,&x,x_inter,&h,system);


#ifdef CalNeuOscSUN_DEBUG
      if(count%count_step == 0){
	printf("x_current : %lf \n", x/tunit);
	printf("h : %g \n", h/tunit);
      }
#endif

      if( gsl_status != GSL_SUCCESS ){
	fprintf(stderr,"CalNeuOscGSL: Error in GSL ODE solver.\n");
	break;
      }

      if(h < h_min){
        h = h_min;
      }

      if(h > h_max){
        h = h_max;
      }

#ifdef CalNeuOscSUN_DEBUG
      count++;
#endif
    }

    x = x_end;

#ifdef CalNeuOscSUN_DEBUG
    printf("End calculation. x_final :  %lf \n",x/tunit);
#endif

  }else{
    EvolveProjectors(tf);
    t=tf;
  }
  return GSL_SUCCESS;
}


int RHS(double t ,const double *state_dbl_in,double *state_dbl_out,void *par){
  SQUIDS *dms=(SQUIDS*)par;
  dms->set_deriv_system_pointer(state_dbl_out);
  dms->Derive(t);
  return  0;
}

