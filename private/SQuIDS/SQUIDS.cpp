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


#include "SQuIDS/SQUIDS.h"

//#define CalNeuOscSUN_DEBUG

SQUIDS::SQUIDS(void){
  is_init=false;
}

SQUIDS::SQUIDS(int n,int ns,int nrh,int nsc, double ti = 0.0){
  ini(n,ns,nrh,nsc,ti);
}

void SQUIDS::ini(int n,int nsu,int nrh,int nsc, double ti=0.0){
  CoherentInt=false;
  NonCoherentInt=false;
  OtherInt=false;
  ScalarsInt=false;
  neu_and_aneu = false;
  AnyNumerics=(CoherentInt||NonCoherentInt||OtherInt||ScalarsInt);

  adaptive_step=true;

  /*
    Setting the number of energy bins, number of components for the density matrix and
    number of scalar functions
  */

  nx=n;
  nsun=nsu;
  nrhos=nrh;
  nscalars=nsc;
  size_rho=nsun*nsun;
  size_state=(size_rho*nrhos+nscalars);
  /*
    Setting time units and initial time
  */
  tunit=1;
  t_ini=ti;
  t=ti;
  nsteps=1000;

  //Allocate memeroy for the system
  Nsystem = nx*size_state;
  numeqn=Nsystem;
  system.reset(new double[Nsystem]);

  /*
    Initializing the SU algebra object, needed to compute algebraic operations like commutators,
    anticommutators, rotations and evolutions
  */

  //Allocate memory
  x.reset(new double[nx]);
  delx.reset(new double[nx]);

  state.reset(new SU_state[nx]);
  dstate.reset(new SU_state[nx]);

  for(int ei = 0; ei < nx; ei++){
    state[ei].rho.reset(new SU_vector[nrhos]);
    dstate[ei].rho.reset(new SU_vector[nrhos]);
  }

  for(int ei = 0; ei < nx; ei++){
    for(int i=0;i<nrhos;i++){
      //state[ei].rho[i].InitSU_vector(nsun,&(system[ei*size_state+i*size_rho]));
      state[ei].rho[i]=SU_vector(nsun,&(system[ei*size_state+i*size_rho]));
      dstate[ei].rho[i]=SU_vector(nsun,nullptr);
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
  h_max    = tunit*1e10;

  // setting up GSL ODE solver

  step = (gsl_odeiv2_step_type*)gsl_odeiv2_step_rkf45;

  sys.function = &RHS;
  sys.jacobian = NULL;
  sys.dimension = (size_t)numeqn;
  sys.params = this;

  is_init=true;
};

void SQUIDS::set_deriv_system_pointer(double *p){
  deriv_system=p;
  for(int ei = 0; ei < nx; ei++){
    for(int i=0;i<nrhos;i++){
      dstate[ei].rho[i].SetBackingStore(&(p[ei*size_state+i*size_rho]));
    }
    dstate[ei].scalar=&(p[ei*size_state+nrhos*size_rho]);
  }
}

SQUIDS::~SQUIDS(){}

/*
  This functions sets the grid of points in where we have a
  density matrix and a set of scalars
 */

int SQUIDS::Set_xrange(double xi, double xf, std::string type){
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
      throw std::runtime_error("SQUIDS::Set_xrange : Xmin too small for log scale");
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
    throw std::runtime_error("SQUIDS::Set_xrange : Not well deffined X range");
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
	throw std::runtime_error("SQUIDS::GetExpectationValueD : x value not in the array.");
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

  if(xi>xr || xi<xl)
    throw std::runtime_error(" Error SQUIDS::Get_i :  value  out of bounds");

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

void SQUIDS::Set_GSL_step(const gsl_odeiv2_step_type * opt){
  step = (gsl_odeiv2_step_type *) opt;  
}

void SQUIDS::Set_AdaptiveStep(bool opt){
  adaptive_step=opt;
}
void SQUIDS::Set_CoherentInteractions(bool opt){
  CoherentInt=opt;
  AnyNumerics=(CoherentInt||NonCoherentInt||OtherInt||ScalarsInt);
}
void SQUIDS::Set_NonCoherentInteractions(bool opt){
  NonCoherentInt=opt;
  AnyNumerics=(CoherentInt||NonCoherentInt||OtherInt||ScalarsInt);
}
void SQUIDS::Set_OtherInteractions(bool opt){
  OtherInt=opt;
  AnyNumerics=(CoherentInt||NonCoherentInt||OtherInt||ScalarsInt);
}
void SQUIDS::Set_ScalarInteractions(bool opt){
  ScalarsInt=opt;
  AnyNumerics=(CoherentInt||NonCoherentInt||OtherInt||ScalarsInt);
}

void SQUIDS::Set_h_min(double opt){
  h_min=opt;
  if(h<h_min){
    if(h_max<50.0*h_min){
      h=(h_min+h_max)/2.0;
    }else{
      h=h_min*10.0;
    }
  }
}
void SQUIDS::Set_h_max(double opt){
  h_max=opt;
  if(h>h_max){
    if(h_max<50.0*h_min){
      h=(h_min+h_max)/2.0;
    }else{
      h=h_min*10.0;
    }
  }
}

void SQUIDS::Set_h(double opt){
  h=opt;
}

void SQUIDS::Set_rel_error(double opt){
  rel_error=opt;
}

void SQUIDS::Set_abs_error(double opt){
  abs_error=opt;
}

void SQUIDS::Set_t(double opt){
  t=opt;
  t_ini=opt;
}

void SQUIDS::Set_units(double opt){
  tunit=opt;
}

void SQUIDS::Set_nx(int opt){
  if(opt!=nx)
    ini(opt,nsun,nrhos,nscalars);
}

void SQUIDS::Set_nsun(int opt){
  if(opt!=nsun)
    ini(nx,opt,nrhos,nscalars);
}

void SQUIDS::Set_NumSteps(int opt){
  nsteps=opt;
}

void SQUIDS::Set_nrhos(int opt){
  if(opt!=nrhos)
    ini(nx,nsun,opt,nscalars);
}
void SQUIDS::Set_nscalars(int opt){
  ini(nx,nsun,nrhos,opt);
}

int SQUIDS::Derive(double at){
  t=at;
  PreDerive(at);
  for(int ei = 0; ei < nx; ei++){
    for(int i = 0; i < nrhos; i++){
      index_rho=i;
      // Density matrix
      // Coherent interaction
      if(CoherentInt)
        dstate[ei].rho[i] = iCommutator(state[ei].rho[i],HI(ei,t));
      
      // Non coherent interaction
      if(NonCoherentInt)
        dstate[ei].rho[i] -= ACommutator(GammaRho(ei,t),state[ei].rho[i]);
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

int SQUIDS::EvolveSUN(double dt){
  if(AnyNumerics){
    int gsl_status = GSL_SUCCESS;


#ifdef CalNeuOscSUN_DEBUG
    printf("GSL paramers :\n");
    
    printf("x_ini : %lf \n", t/tunit);
    printf("x_end : %lf \n", (t+dt)/tunit);
    printf("h : %g \n", h/tunit);
    printf("h_min : %g \n", h_min/tunit);
#endif

    // initial time
    
    // ODE system error control
    gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys,step,h,abs_error,rel_error);
    gsl_odeiv2_driver_set_hmin(d,h_min);
    gsl_odeiv2_driver_set_hmax(d,h_max);
    gsl_odeiv2_driver_set_nmax(d,0);

#ifdef CalNeuOscSUN_DEBUG
    printf("Start calculation.\n");
#endif
    
    double* gsl_sys = system.get();
    
    if(adaptive_step){
      gsl_status = gsl_odeiv2_driver_apply(d, &t, t+dt, gsl_sys);
    }else{
      gsl_status = gsl_odeiv2_driver_apply_fixed_step(d, &t, dt/nsteps , nsteps , gsl_sys);
    }
    
    gsl_odeiv2_driver_free(d);
    if( gsl_status != GSL_SUCCESS ){
      throw std::runtime_error("SQUIDS::EvolveSUN: Error in GSL ODE solver.");
    }

#ifdef CalNeuOscSUN_DEBUG
    printf("End calculation. x_final :  %lf \n",x/tunit);
#endif
  }else{
    t+=dt;
    PreDerive(t);
  }

  return GSL_SUCCESS;
}

int RHS(double t ,const double *state_dbl_in,double *state_dbl_out,void *par){
  SQUIDS *dms=(SQUIDS*)par;
  dms->set_deriv_system_pointer(state_dbl_out);
  dms->Derive(t);
  return 0;
}
