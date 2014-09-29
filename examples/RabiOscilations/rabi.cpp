#include "rabi.h"


void rabi::init(double D_E, double wi, double Am){
  Delta_E=D_E;
  w=wi;
  A=Am;
  Set("th12",params.pi/4.0);  
  ini(1,2,1,0);

  Set("H1",true);  

  evol_b0_proj.reset(new SU_vector[nx*nsun]);
  evol_b1_proj.reset(new SU_vector[nx*nsun]);
  b0_proj.reset(new SU_vector[nsun]);
  b1_proj.reset(new SU_vector[nsun]);

  for(int i = 0; i < nsun; i++){
    b0_proj[i].InitSU_vector("Proj",i,nsun);
    b1_proj[i].InitSU_vector("Proj",i,nsun);
    b1_proj[i].RotateToB1(&params);

    for(int ei = 0; ei < nx; ei++){
      evol_b0_proj[i*nx + ei].InitSU_vector("Proj",i,nsun);
      evol_b1_proj[i*nx + ei].InitSU_vector("Proj",i,nsun);
      evol_b1_proj[i*nx + ei].RotateToB1(&params);
    }
  }


  suH0.InitSU_vector(nsun);
  d.InitSU_vector(nsun);
  d0.InitSU_vector(nsun);

  d0=(b1_proj[0]-b1_proj[1]);
  suH0 = b0_proj[1]*Delta_E;

  // set initial conditions for the density mattrix.
  state[0].rho[0] = b0_proj[0];
    
}


void rabi::PreDerive(double t){
  SU_vector h0(nsun);
  for(int i = 0; i < nsun; i++){
    for(int ei = 0; ei < nx; ei++){
      //h0=H0(x[ei]);
      evol_b0_proj[i*nx + ei] = b0_proj[i].SUEvolve(suH0,t-t_ini);
      evol_b1_proj[i*nx + ei] = b1_proj[i].SUEvolve(suH0,t-t_ini);
      //evol_b0_proj[i*nx + ei] = b0_proj[i].SUEvolve(h0,t-t_ini);      
    }
  }
}



SU_vector rabi::H0(double x){
  return suH0;
}


SU_vector rabi::HI(int ix,double t){
  d=(evol_b1_proj[0]-evol_b1_proj[1]);
  //  cout << d << endl;
  return (d*A)*cos(w*t);
}


void rabi::set_evol(void){
  for(int i = 0; i < nx; i++){
    SU_vector h0=H0(x[i]);
    for(int nrh=0;nrh<nrhos;nrh++){
      state[i].rho[nrh]=state[i].rho[nrh].SUEvolve(h0,(t_ini-t_end));
    }
  }
}
