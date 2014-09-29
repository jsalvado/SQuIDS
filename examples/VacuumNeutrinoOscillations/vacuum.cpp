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

 /******************************************************************************
 / Very simple example that computes the evolution for vacuum oscilatios of    /
 / neutrinos. The numerical integration is not needed for this since is        /
 / already solved analyticaly using the SUN_vector evolution.                  /
 ******************************************************************************/



#include "vacuum.h"


void vacuum::init(int n,int ns, double Ein, double Efin){
  //initialize SQUID with one density matrix and zero scalar functions
  ini(n,ns,1,0);
  //initialize the SU_vector that contines the delta_m
  DM2.InitSU_vector(nsun);
  //set the energy range in log scale
  Set_xrange(Ein, Efin,"log");


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


  for(int i = 1; i < nsun; i++){
    DM2 += (b0_proj[i])*gsl_matrix_get(params.dmsq,i,0);
  }


  //set initial conditions for the density mattrix.  
  for(int ei = 0; ei < nx; ei++){
    state[ei].rho[0]=b1_proj[0];
    //    state[ei].rho[0].RotateToB1(&params);
    
    for(int i=0;i<nscalars;i++){
      state[ei].scalar[i]=0;
    }
  }
}








//funtion that sets the parameters for the Vacuum oscilation 
//is necesari just to update the mass matrix, and mixings every time.
void vacuum::SetVacuum(string name,double opt){
  Set(name,opt);
  DM2 = (b0_proj[1])*gsl_matrix_get(params.dmsq,1,0);
  for(int i = 2; i < nsun; i++){
    DM2 += (b0_proj[i])*gsl_matrix_get(params.dmsq,i,0);
  }
  for(int i = 0; i < nsun; i++){
    b1_proj[i]=b0_proj[i];
    b1_proj[i].RotateToB1(&params);
    
    for(int ei = 0; ei < nx; ei++){
      evol_b1_proj[i*nx + ei]=b0_proj[i];
      evol_b1_proj[i*nx + ei].RotateToB1(&params);
    }
  }
  
  for(int ei = 0; ei < nx; ei++){
    state[ei].rho[0]=b1_proj[0];
    state[ei].rho[0].RotateToB1(&params);
  }
}  


//Function that returns the H0 operator
SU_vector vacuum::H0(double x){
  return DM2*(0.5/x);
}


//function that returns the flux for the neutrino with flavo "i" and energy "e"
double vacuum::Get_flux(int i,double e){
  return GetExpectationValueD(b1_proj[i],0,e);
}


//function that returns the flux for the neutrino with flavo "i" and energy "e"
double vacuum::Get_flux(int i,int e){
  return GetExpectationValue(b1_proj[i],0,e);
}


void vacuum::set_evol(void){
  for(int i = 0; i < nx; i++){
    SU_vector h0=H0(x[i]);
    for(int nrh=0;nrh<nrhos;nrh++){
      state[i].rho[nrh]=state[i].rho[nrh].SUEvolve(h0,(t_ini-t_end));
    }
  }
}


