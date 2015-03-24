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
 / Very simple example that computes the evolution for vacuum oscilatios of    /
 / neutrinos. The numerical integration is not needed for this since is        /
 / already solved analyticaly using the SUN_vector evolution.                  /
 ******************************************************************************/


#include "vacuum.h"
#include <fstream>

int main(){
  squids::Const units;
  
  //Number of energy bins
  unsigned int Nenergy=1000;
  //Number of flavors
  unsigned int Nflavor=3;
  //Energy Range
  double Emin=10*units.MeV, Emax=10*units.GeV;
  //declaration of the object
  vacuum V0(Nenergy,Nflavor,Emin,Emax);

  V0.Evolve(1000*units.km);

  std::ofstream file("oscillations.dat");

  const int nu_e=0, nu_mu=1, nu_tau=2;
  for(double lE=log(Emin); lE<log(Emax); lE+=0.0001){
    double E=exp(lE);
    file << E/units.GeV << "  " << V0.Get_flux(nu_e,E) << "  " <<
      V0.Get_flux(nu_mu,E) << "  " << V0.Get_flux(nu_tau,E) << std::endl;
  }

  std::cout << std::endl <<  "Done! " << std::endl <<  "Do you want to run the gnuplot script? yes/no" << std::endl;
  std::string plt;
  std::cin >> plt;

  if(plt=="yes" || plt=="y"){
    return system("./plot.plt");
  }
  
  return 0;
}
