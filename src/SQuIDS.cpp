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
 *      Christopher Weaver (University of Wisconsin Madison)                   *
 *         chris.weaver@icecube.wisc.edu                                       *
 ******************************************************************************/

#include "SQuIDS.h"
#include <cmath>
#include <limits>
#include <algorithm>

///\brief Auxiliary function used for the GSL interface
int RHS(double ,const double*,double*,void*);

SQuIDS::SQuIDS():
CoherentRhoTerms(false),
NonCoherentRhoTerms(false),
OtherRhoTerms(false),
GammaScalarTerms(false),
OtherScalarTerms(false),
AnyNumerics(false),
is_init(false),
adaptive_step(true),
nsteps(1000),
step(gsl_odeiv2_step_rkf45),
h(std::numeric_limits<double>::epsilon()),
h_min(std::numeric_limits<double>::min()),
h_max(std::numeric_limits<double>::max()),
abs_error(1e-20),
rel_error(1e-20)
{
  sys.function = &RHS;
  sys.jacobian = NULL;
  sys.dimension = 0;
  sys.params = this;
}

SQuIDS::SQuIDS(unsigned int n, unsigned int ns, unsigned int nrh, unsigned int nsc, double ti):
SQuIDS(){
  ini(n,ns,nrh,nsc,ti);
}

SQuIDS::SQuIDS(SQuIDS&& other):
CoherentRhoTerms(other.CoherentRhoTerms),
NonCoherentRhoTerms(other.NonCoherentRhoTerms),
OtherRhoTerms(other.OtherRhoTerms),
GammaScalarTerms(other.GammaScalarTerms),
OtherScalarTerms(other.OtherScalarTerms),
AnyNumerics(other.AnyNumerics),
is_init(other.is_init),
adaptive_step(other.adaptive_step),
x(std::move(other.x)),
t(other.t),
t_ini(other.t_ini),
nsteps(other.nsteps),
size_rho(other.size_rho),
size_state(other.size_state),
system(std::move(other.system)),
step(other.step),
sys(other.sys),
h(other.h),
h_min(other.h_min),
h_max(other.h_max),
abs_error(other.abs_error),
rel_error(other.rel_error),
dstate(std::move(other.dstate)),
nx(other.nx),
nsun(other.nsun),
nrhos(other.nrhos),
nscalars(other.nscalars),
params(std::move(other.params)),
state(std::move(other.state))
{
  sys.params=this;
  other.is_init=false; //other is no longer usable, since we stole its contents
}

void SQuIDS::ini(unsigned int n, unsigned int nsu, unsigned int nrh, unsigned int nsc, double ti){
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
  t_ini=ti;
  t=ti;

  //Allocate memeroy for the system
  unsigned int numeqn=nx*size_state;
  system.reset(new double[numeqn]);
  sys.dimension = static_cast<size_t>(numeqn);

  /*
    Initializing the SU algebra object, needed to compute algebraic operations like commutators,
    anticommutators, rotations and evolutions
  */

  //Allocate memory
  x.resize(nx);

  state.reset(new SU_state[nx]);
  dstate.reset(new SU_state[nx]);

  for(unsigned int ei = 0; ei < nx; ei++){
    state[ei].rho.reset(new SU_vector[nrhos]);
    dstate[ei].rho.reset(new SU_vector[nrhos]);
  }

  for(unsigned int ei = 0; ei < nx; ei++){
    for(unsigned int i=0;i<nrhos;i++){
      state[ei].rho[i]=SU_vector(nsun,&(system[ei*size_state+i*size_rho]));
      dstate[ei].rho[i]=SU_vector(nsun,nullptr);
    }
    if(nscalars>0){
      state[ei].scalar=&(system[ei*size_state+nrhos*size_rho]);
    }
  }

  is_init=true;
};

void SQuIDS::set_deriv_system_pointer(double *p){
  deriv_system=p;
  for(unsigned int ei = 0; ei < nx; ei++){
    for(unsigned int i=0;i<nrhos;i++){
      dstate[ei].rho[i].SetBackingStore(&(p[ei*size_state+i*size_rho]));
    }
    dstate[ei].scalar=&(p[ei*size_state+nrhos*size_rho]);
  }
}

SQuIDS::~SQuIDS(){}

SQuIDS& SQuIDS::operator=(SQuIDS&& other){
  if(&other==this)
    return(*this);
  
  CoherentRhoTerms=other.CoherentRhoTerms;
  NonCoherentRhoTerms=other.NonCoherentRhoTerms;
  OtherRhoTerms=other.OtherRhoTerms;
  GammaScalarTerms=other.GammaScalarTerms;
  OtherScalarTerms=other.OtherScalarTerms;
  AnyNumerics=other.AnyNumerics;
  is_init=other.is_init;
  adaptive_step=other.adaptive_step;
  x=std::move(other.x);
  t=other.t;
  t_ini=other.t_ini;
  nsteps=other.nsteps;
  size_rho=other.size_rho;
  size_state=other.size_state;
  system=std::move(other.system);
  step=other.step;
  sys=other.sys;
  h=other.h;
  h_min=other.h_min;
  h_max=other.h_max;
  abs_error=other.abs_error;
  rel_error=other.rel_error;
  dstate=std::move(other.dstate);
  nx=other.nx;
  nsun=other.nsun;
  nrhos=other.nrhos;
  nscalars=other.nscalars;
  params=std::move(other.params);
  state=std::move(other.state);
  sys.params=this;
  other.is_init=false; //other is no longer usable, since we stole its contents
  
  return(*this);
}

/*
  This functions sets the grid of points in where we have a
  density matrix and a set of scalars
 */

void SQuIDS::Set_xrange(double xi, double xf, std::string type){
  if (xi == xf){
    x[0] = xi;
    return;
  }

  if(type=="linear" || type=="Linear" || type=="lin" || type=="Lin"){
    for(unsigned int e1 = 0; e1 < nx; e1++){
      x[e1]=xi+(xf-xi)*static_cast<double>(e1)/static_cast<double>(nx-1);
    }
  }else if(type=="log" || type=="Log"){
    double xmin_log,xmax_log;
    if (xi < 1.0e-10 ){
      throw std::runtime_error("SQUIDS::Set_xrange : Xmin too small for log scale");
    }else{
        xmin_log = log(xi);
    }
    xmax_log = log(xf);

    for(unsigned int e1 = 0; e1 < nx; e1++){
      double X=xmin_log+(xmax_log-xmin_log)*static_cast<double>(e1)/static_cast<double>(nx-1);
      x[e1]=exp(X);
    }
  }else{
    throw std::runtime_error("SQUIDS::Set_xrange : Not well deffined X range");
  }
}

double SQuIDS::GetExpectationValue(SU_vector op, unsigned int nrh, unsigned int i) const{
  SU_vector h0=H0(x[i],nrh);
  return state[i].rho[nrh]*op.Evolve(h0,t-t_ini);
}

double SQuIDS::GetExpectationValueD(SU_vector op, unsigned int nrh, double xi) const{
  SU_vector h0=H0(xi,nrh);
  int xid=-1;
  for(unsigned int i = 0; i < nx; i++){
    if ( xi >= x[i] && xi <= x[i+1]){
      xid = i;
      break;
    }
  }
  if(xid==-1){
    throw std::runtime_error("SQUIDS::GetExpectationValueD : x value not in the array.");
  }

  return (state[xid].rho[nrh] +
	   (state[xid+1].rho[nrh]-
	    state[xid].rho[nrh])*((xi-x[xid])/(x[xid+1]-x[xid])))*op.Evolve(h0,t-t_ini);
}

void SQuIDS::Set_xrange(const std::vector<double>& xs){
  if(xs.size()!=nx)
    throw std::runtime_error("SQUIDS::Set_xrange : wrong number of x values");
  if(!std::is_sorted(xs.begin(),xs.end()))
    throw std::runtime_error("SQUIDS::Set_xrange : x values must be sorted");
  x=xs;
}

unsigned int SQuIDS::Get_i(double xi) const{
  double xl, xr;
  unsigned int nr=nx-1;
  unsigned int nl=0;

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

void SQuIDS::Set_GSL_step(gsl_odeiv2_step_type const* opt){
  step = opt;
}

void SQuIDS::Set_AdaptiveStep(bool opt){
  adaptive_step=opt;
}
void SQuIDS::Set_CoherentRhoTerms(bool opt){
  CoherentRhoTerms=opt;
  AnyNumerics=(CoherentRhoTerms||NonCoherentRhoTerms||OtherRhoTerms||GammaScalarTerms||OtherScalarTerms);
}
void SQuIDS::Set_NonCoherentRhoTerms(bool opt){
  NonCoherentRhoTerms=opt;
  AnyNumerics=(CoherentRhoTerms||NonCoherentRhoTerms||OtherRhoTerms||GammaScalarTerms||OtherScalarTerms);
}
void SQuIDS::Set_OtherRhoTerms(bool opt){
  OtherRhoTerms=opt;
  AnyNumerics=(CoherentRhoTerms||NonCoherentRhoTerms||OtherRhoTerms||GammaScalarTerms||OtherScalarTerms);
}
void SQuIDS::Set_GammaScalarTerms(bool opt){
  GammaScalarTerms=opt;
  AnyNumerics=(CoherentRhoTerms||NonCoherentRhoTerms||OtherRhoTerms||GammaScalarTerms||OtherScalarTerms);
}
void SQuIDS::Set_OtherScalarTerms(bool opt){
  OtherScalarTerms=opt;
  AnyNumerics=(CoherentRhoTerms||NonCoherentRhoTerms||OtherRhoTerms||GammaScalarTerms||OtherScalarTerms);
}

void SQuIDS::Set_h_min(double opt){
  h_min=opt;
  if(h<h_min){
    if(h_max<50.0*h_min){
      h=(h_min+h_max)/2.0;
    }else{
      h=h_min*10.0;
    }
  }
}
void SQuIDS::Set_h_max(double opt){
  h_max=opt;
  if(h>h_max){
    if(h_max<50.0*h_min){
      h=(h_min+h_max)/2.0;
    }else{
      h=h_min*10.0;
    }
  }
}

void SQuIDS::Set_h(double opt){
  h=opt;
}

void SQuIDS::Set_rel_error(double opt){
  rel_error=opt;
}

void SQuIDS::Set_abs_error(double opt){
  abs_error=opt;
}

void SQuIDS::Set_NumSteps(unsigned int opt){
  nsteps=opt;
}

int SQuIDS::Derive(double at){
  t=at;
  PreDerive(at);
  for(unsigned int ei = 0; ei < nx; ei++){
    // Density matrix
    for(unsigned int i = 0; i < nrhos; i++){
      // Coherent interaction
      if(CoherentRhoTerms)
        dstate[ei].rho[i] = iCommutator(state[ei].rho[i],HI(ei,i,t));
      else
        dstate[ei].rho[i].SetAllComponents(0.);
      
      // Non coherent interaction
      if(NonCoherentRhoTerms)
        dstate[ei].rho[i] -= ACommutator(GammaRho(ei,i,t),state[ei].rho[i]);
      // Other possible interaction, for example involving the Scalars or non linear terms in rho.
      if(OtherRhoTerms)
        dstate[ei].rho[i] += InteractionsRho(ei,i,t);
    }
    //Scalars
    for(unsigned int is=0;is<nscalars;is++){
      dstate[ei].scalar[is]=0.;
      if(GammaScalarTerms)
        dstate[ei].scalar[is] += -state[ei].scalar[is]*GammaScalar(ei,is,t);
      if(OtherScalarTerms)
        dstate[ei].scalar[is] += InteractionsScalar(ei,is,t);
    }
  }

  return GSL_SUCCESS;
}

int SQuIDS::Evolve(double dt){
  if(AnyNumerics){
    int gsl_status = GSL_SUCCESS;

    // initial time
    
    // ODE system error control
    gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys,step,h,abs_error,rel_error);
    gsl_odeiv2_driver_set_hmin(d,h_min);
    gsl_odeiv2_driver_set_hmax(d,h_max);
    gsl_odeiv2_driver_set_nmax(d,0);
    
    double* gsl_sys = system.get();
    
    if(adaptive_step){
      gsl_status = gsl_odeiv2_driver_apply(d, &t, t+dt, gsl_sys);
    }else{
      gsl_status = gsl_odeiv2_driver_apply_fixed_step(d, &t, dt/nsteps , nsteps , gsl_sys);
    }
    
    gsl_odeiv2_driver_free(d);
    if( gsl_status != GSL_SUCCESS ){
      throw std::runtime_error("SQUIDS::Evolve: Error in GSL ODE solver.");
    }
  }else{
    t+=dt;
    PreDerive(t);
  }

  return GSL_SUCCESS;
}

int RHS(double t ,const double *state_dbl_in,double *state_dbl_out,void *par){
  SQuIDS *dms=static_cast<SQuIDS*>(par);
  dms->set_deriv_system_pointer(state_dbl_out);
  dms->Derive(t);
  return 0;
}
