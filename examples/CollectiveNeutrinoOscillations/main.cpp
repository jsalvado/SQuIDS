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

#include "collective.h"
#include <fstream>



int main(){
  //Parameters
  double mu=10.0;
  double mu2=10;
  double wmin=-2;
  double wmax=2;
  double th=0.01;
  double th2=0.01;
  int Nbins=200;

  
  collective ColNus(mu,th,wmin,wmax,Nbins);
  collective ColNus2(mu2,th2,wmin,wmax,Nbins);
  
  std::ofstream file;


  //Evolution from mu=10 to mu=0 in a time period of 100
  ColNus.Adiabatic_mu(10,0,100,true);

  //write the ouput in the file
  //col 1 value of w
  //col 2 expectation value of ez for the evolved system
  //col 3 expectation value of ez for the non evolved system
  SU_vector o=ColNus.ez;
  file.open("collective.dat");
  for(int w=0;w<Nbins;w++){
    file << std::scientific << ColNus.Get_x(w) << "\t" 
	 << ColNus.GetExpectationValue(o,0,w)<<"  " <<  
      ColNus2.GetExpectationValue(o,0,w) << std::endl;
  }
  
  file.close();

  //runs the gnuplot script if yes
  std::string plt;
  std::cout << std::endl <<  "Done! " << std::endl <<  "Do you want to run the gnuplot script? yes/no" << std::endl;
  std::cin >> plt;
  
  if(plt=="yes" || plt=="y"){
    return system("./plot.plt");
  }
  
  return 0;


}
