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

 /******************************************************************************
  Very simple example that computes the evolution for vacuum oscilatios of 
  neutrinos. The numerical integration is not needed for this since is already
  solved using the SUN_vector evolution.
 ******************************************************************************/



#ifndef __VACUUM_H
#define __VACUUM_H

#include <iostream>
#include <float.h>
#include <math.h>
#include <complex>
#include <vector>

#include <SQUIDS.h>

class vacuum: public SQUIDS {
 private:
  SU_vector DM2;
 public:
  std::unique_ptr<SU_vector[]> b0_proj;
  std::unique_ptr<SU_vector[]> b1_proj;
  std::unique_ptr<SU_vector[]> evol_b0_proj;
  std::unique_ptr<SU_vector[]> evol_b1_proj;


  vacuum(void){};
  //constructor and initialization
  //int -> number of energy bins
  //int -> number of flavors
  //double -> initial energy
  //double -> final energy
  vacuum(int n1,int n2,double d1, double d2){init(n1,n2,d1,d2);};
  void init(int,int,double, double);

  //Function to set the parameters.
  void SetVacuum(Const par);  
  //H0 operator
  SU_vector H0(double);
  void EvolveProjectors(double t);
  void set_evol(void);

  
  //get the final flux, the initial one is flat in energy
  double Get_flux(int,double);
  double Get_flux(int,int);

};

#endif
