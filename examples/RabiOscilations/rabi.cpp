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
#include <cmath>

void rabi::init(double D_E, double wi, double Am){
  params.SetEnergyDifference(1,D_E);
  w=wi;
  A=Am;
  params.SetMixingAngle(0,1,params.pi/4);
  ini(1,2,1,0,0);

  Set_CoherentRhoTerms(true);

  b0_proj.reset(new SU_vector[nsun]);
  b1_proj.reset(new SU_vector[nsun]);

  for(int i = 0; i < nsun; i++){
    b0_proj[i]=SU_vector::Projector(nsun,i);
    b1_proj[i]=SU_vector::Projector(nsun,i);
    b1_proj[i].RotateToB1(params);
  }

  suH0=SU_vector(nsun);
  suH0 = b0_proj[1]*params.GetEnergyDifference(1);
  
  d0=(b1_proj[0]-b1_proj[1]);
  d.reset(new SU_vector[nx]);
  for(unsigned int i=0; i<nx; i++)
    d[i]=d0;

  // set initial conditions for the density matrix.
  state[0].rho[0] = b0_proj[0];
}

void rabi::PreDerive(double t){
  for(int ei = 0; ei < nx; ei++)
    d[ei] = d0.Evolve(suH0,t-Get_t_initial());
}

SU_vector rabi::H0(double x, unsigned int irho) const{
  return suH0;
}

SU_vector rabi::HI(unsigned int ix, unsigned int irho, double t) const{
  return (A*cos(w*t))*d[ix];
}
