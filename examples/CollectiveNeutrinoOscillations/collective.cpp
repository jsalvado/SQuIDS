#include "collective.h"


void collective::init(double m, double wmin, double wmax, int Nbins){
  mu=m;
  ini(Nbins,2,1,0,0.0);
  w_min=wmin;
  w_max=wmax;

  Set_xrange(wmin,wmax,"lin");

  //  params.SetMixingAngle(0,1,params.pi/4.0);  
  Set_CoherentInteractions(true);  
  ex=SU_vector::Component(nsun,1);
  ey=SU_vector::Component(nsun,2);
  ez=SU_vector::Component(nsun,3);

  // evol_b0_proj.reset(new SU_vector[nx*nsun]);
  // evol_b1_proj.reset(new SU_vector[nx*nsun]);
  // b0_proj.reset(new SU_vector[nsun]);
  // b1_proj.reset(new SU_vector[nsun]);

  // for(int i = 0; i < nsun; i++){
  //   b0_proj[i]=SU_vector::Projector(nsun,i);
  //   b1_proj[i]=SU_vector::Projector(nsun,i);
  //   b1_proj[i].RotateToB1(params);

  //   for(int ei = 0; ei < nx; ei++){
  //     evol_b0_proj[i*nx + ei]=SU_vector::Projector(nsun,i);
  //     evol_b1_proj[i*nx + ei]=SU_vector::Projector(nsun,i);
  //     evol_b1_proj[i*nx + ei].RotateToB1(params);
  //   }
  // }

  B=ez;

  // set initial conditions for the density mattrix.
  for(int ei = 0; ei < nx; ei++){
    double th=0.1;
    double w=Get_x(ei);
    //if(fabs(w)>0.1){
      if(w>0){
	state[ei].rho[0] = (ez*sin(th) + ey*cos(th))*(1/(2*w*w)*Fermi(1/(2.0*w)));
      }else{
	state[ei].rho[0] = (ez*sin(th) + ey*cos(th))*(-1/(2*w*w)*Fermi(-1/(2.0*w)));
      }
    }
  //}
    
}


void collective::PreDerive(double t){
  P=state[0].rho[0];
  for(int ei = 1; ei < nx; ei++){
    P+=state[ei].rho[0];
  }
}


SU_vector collective::H0(double x){
  return B*x;
}


SU_vector collective::HI(int ix,double t){
  return P*(mu*(w_max-w_min));
  //std::cout << P << std::endl;
  //return P;
}


double collective::Fermi(double EoverT){
  return 1.0/(exp(EoverT)+1.0);
}
