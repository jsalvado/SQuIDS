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


/*******************************************************************************
 * Example of the collective phenomena for neturinos,                          *
 * based on Georg G. Raffelt  Phys.Rev. D83 (2011) 105022.                     *
 * arXiv:1103.2891                                                             *
 ******************************************************************************/

#ifndef __COLLECTIVE_H
#define __COLLECTIVE_H




#include <iostream>
#include <float.h>
#include <math.h>
#include <complex>
#include <vector>

#include <SQUIDS.h>

class collective: public SQUIDS {
 private:
  //Hamiltonian vectors for the H0 part and the self interacting term
  SU_vector B;
  SU_vector P;

  //self interaction strengh sqrt(2) G_F n_\nu
  double mu;
  double mu_i,mu_f;
  
  //range of the spectrum w=T/E
  double w_min, w_max;
  
  //The object is designed to produce the evolution varying mu linearly this is the value of 
  //the time periond for this variation.
  double period;
  
  //Mixing angle
  double theta;
  
  //sets on and off the progress bar
  bool bar;
  
  //shows the progress bar
  void progressbar( int percent, double mu);

  //function that is evaluated before computting the derivatives, it bassically computes the vector P 
  void PreDerive(double t);

  //Hamilitonian of the system, in this case we don't use any time idependent separation(non in the interaction picutre formalism)
  SU_vector HI(int ix,double t);

  //Fermi distribution
  double Fermi(double EoverT);

 public:
  //basis
  SU_vector ex,ey,ez;

  collective(void){};

  //Constructor and initializer
  // mu -> value of the self interacting strengh, is proportional to the neutrino density
  // th -> Mixing angle
  // wmin and wmax -> deffines the range of values for w
  // Nbins number of bins in the spectra
  collective(double mu,double th, double wmin, double wmax, int Nbins){init(mu,th, wmin, wmax, Nbins);};
  void init(double mu,double th, double wmin, double wmax, int Nbins);

  //Function that sets the value of mu
  void Set_mu(double m){mu=m;}

  //Function that solve the evolution changing mu from mu_i to mu_f in the time period given by period
  // bar shows the progress bar if it's true.
  void Adiabatic_mu(double mu_i, double mu_f, double period, bool bar);

};

#endif
