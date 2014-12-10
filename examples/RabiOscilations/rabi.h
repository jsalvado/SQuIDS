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

 /*****************************************************************************
 / Simple example computting two rabi models with and whithout dettuning       /
 / reproducing some plots form Phys. Rev. B 84, 075107 J. I. Fuks, N. Helbig,  /
 / I. V. Tokatly, and A. Rubio arXiv:1101.2880                                 /
 ******************************************************************************/


#ifndef __RABI_H
#define __RABI_H

#include <iostream>
#include <float.h>
#include <math.h>
#include <complex>
#include <vector>

#include <SQUIDS.h>

class rabi: public SQUIDS {
 private:
  //Hamiltonian no external field
  SU_vector suH0;

  //Energy difference
  double Delta_E;
  //Laser Frequency
  double w;
  //Laser Amplitude
  double A;

 public:
  //Dipole operator
  SU_vector d0;
  //Projectors in the different basis
  std::unique_ptr<SU_vector[]> b0_proj;
  std::unique_ptr<SU_vector[]> b1_proj;
  std::unique_ptr<SU_vector[]> evol_b0_proj;
  std::unique_ptr<SU_vector[]> evol_b1_proj;

  //Function evaluated before any derivative
  void PreDerive(double t);
  //Constructors and initializer
  // 
  rabi(){};
  rabi(double D_E, double wi, double Am){init(D_E,wi,Am);};
  void init(double D_E, double wi, double Am);
  
  //Time independent hamiltonian
  SU_vector H0(double) const;
  //Time dependent hamiltonian(Laser)
  SU_vector HI(int ix,double t) const;

};

#endif
