#include "collective.h"


void collective::init(double m, double wmin, double wmax, int Nbins){
  mu=m;
  ini(Nbins,2,1,0);
  Set_xrange(wmin,wmax,"lin");

  //  params.Set_th12(params.pi/4.0);  
  Set_CoherentInteractions(true);  

  evol_b0_proj.reset(new SU_vector[nx*nsun]);
  evol_b1_proj.reset(new SU_vector[nx*nsun]);
  b0_proj.reset(new SU_vector[nsun]);
  b1_proj.reset(new SU_vector[nsun]);

  for(int i = 0; i < nsun; i++){
    b0_proj[i]=SU_vector::Projector(nsun,i);
    b1_proj[i]=SU_vector::Projector(nsun,i);
    b1_proj[i].RotateToB1(&params);

    for(int ei = 0; ei < nx; ei++){
      evol_b0_proj[i*nx + ei]=SU_vector::Projector(nsun,i);
      evol_b1_proj[i*nx + ei]=SU_vector::Projector(nsun,i);
      evol_b1_proj[i*nx + ei].RotateToB1(&params);
    }
  }

  B=b0_proj[1];

  // set initial conditions for the density mattrix.
  state[0].rho[0] = b1_proj[0];
    
}


void collective::PreDerive(double t){
  for(int i = 0; i < nsun; i++){
    for(int ei = 0; ei < nx; ei++){
      evol_b0_proj[i*nx + ei] = b0_proj[i].SUEvolve(B,t-t_ini);
      evol_b1_proj[i*nx + ei] = b1_proj[i].SUEvolve(B,t-t_ini);
    }
  }
}



SU_vector collective::H0(double x){
  return B;
}


SU_vector collective::HI(int ix,double t){
  P=state[0].rho[0];
  for(int ei = 0; ei < nx; ei++){
    P+=state[0].rho[0];
  }
  P=P*(mu/(double)nx);
  return P;
}

