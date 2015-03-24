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
  Very simple example that computes the evolution for vacuum oscillations of 
  neutrinos. The numerical integration is not needed for this since is already
  solved using the SUN_vector evolution.
 ******************************************************************************/

#ifndef VACUUM_H
#define VACUUM_H

#include <iostream>
#include <float.h>
#include <math.h>
#include <complex>
#include <vector>

#include <SQuIDS.h>

class vacuum: public squids::SQuIDS {
 private:
  squids::SU_vector DM2;
 public:
  std::unique_ptr<squids::SU_vector[]> b0_proj;
  std::unique_ptr<squids::SU_vector[]> b1_proj;

  vacuum(){};
  //constructor and initialization
  //nbins -> number of energy bins
  //nflavor -> number of flavors
  //Ein -> initial energy
  //Efin -> final energy
  vacuum(unsigned int nbins, unsigned int nflavor, double Eini, double Efin){init(nbins,nflavor,Eini,Efin);};
  void init(unsigned int nbins, unsigned int nflavor, double Eini, double Efin);

  //Function to set the parameters.
  void SetVacuum(squids::Const par);
  //H0 operator
  squids::SU_vector H0(double E, unsigned int irho) const;
  
  //get the final flux, the initial one is flat in energy
  double Get_flux(int,double);
};

#endif
