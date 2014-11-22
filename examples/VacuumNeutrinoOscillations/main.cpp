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
#include <fstream>

#define Kilometer 0.197

int main(){
  //declaration of the object
  vacuum V0;
  //Number of energy bins
  int Nenergy=1000;
  //Name of the output file
  std::string plt;
  //Initialization of the object
  V0.init(Nenergy,3,0.0005,10);

  V0.EvolveSUN(0,1000*Kilometer);

  std::ofstream file("oscillations.dat");

  for(double lE=log(0.0005); lE<log(10); lE+=0.0001){
    double E=exp(lE);
    file << E << "  " << V0.Get_flux(0,E) << "  " <<
      V0.Get_flux(1,E) << "  " << V0.Get_flux(2,E) << std::endl;
  }

  std::cout << std::endl <<  "Done! " << std::endl <<  "Do you want to run the gnuplot script? yes/no" << std::endl;
  std::cin >> plt;

  if(plt=="yes" || plt=="y"){
    return system("./plot.plt");
  }
  
  return 0;
}
