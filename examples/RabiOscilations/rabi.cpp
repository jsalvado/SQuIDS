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
 *      Christopher Weaver (University of Wisconsin Madison)                   * 
 *         chris.weaver@icecube.wisc.edu                                       *
 *      Jordi Salvado (University of Wisconsin Madison)                        *
 *         jsalvado@icecube.wisc.edu                                           *
 ******************************************************************************/


#include "rabi.h"


void rabi::init(double D_E, double wi, double Am){
  Delta_E=D_E;
  w=wi;
  A=Am;
  params.SetMixingAngle(0,1,params.pi/4);
  ini(1,2,1,0,0);

  Set_CoherentRhoTerms(true);

  evol_b0_proj.reset(new SU_vector[nx*nsun]);
  evol_b1_proj.reset(new SU_vector[nx*nsun]);
  b0_proj.reset(new SU_vector[nsun]);
  b1_proj.reset(new SU_vector[nsun]);

  for(int i = 0; i < nsun; i++){
    b0_proj[i]=SU_vector::Projector(nsun,i);
    b1_proj[i]=SU_vector::Projector(nsun,i);
    b1_proj[i].RotateToB1(params);

    for(int ei = 0; ei < nx; ei++){
      evol_b0_proj[i*nx + ei]=SU_vector::Projector(nsun,i);
      evol_b1_proj[i*nx + ei]=SU_vector::Projector(nsun,i);
      evol_b1_proj[i*nx + ei].RotateToB1(params);
    }
  }

  suH0=SU_vector(nsun);
  d0=(b1_proj[0]-b1_proj[1]);
  suH0 = b0_proj[1]*Delta_E;

  // set initial conditions for the density matrix.
  state[0].rho[0] = b0_proj[0];
}

void rabi::PreDerive(double t){
  for(int i = 0; i < nsun; i++){
    for(int ei = 0; ei < nx; ei++){
      evol_b0_proj[i*nx + ei] = b0_proj[i].Evolve(suH0,t-t_ini);
      evol_b1_proj[i*nx + ei] = b1_proj[i].Evolve(suH0,t-t_ini);
    }
  }
}

SU_vector rabi::H0(double x, unsigned int irho) const{
  return suH0;
}

SU_vector rabi::HI(unsigned int ix, unsigned int irho, double t) const{
  return (A*cos(w*t))*(evol_b1_proj[0]-evol_b1_proj[1]);
}
