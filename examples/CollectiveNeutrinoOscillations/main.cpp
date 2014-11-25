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
        Chris
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

void progressbar( int percent){
  static int last=-1;
  if(percent==last)
    return;
  last=percent;
  std::string bar;
  
  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }
  
  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%  " << std::flush;
}


int main(){
  // Declaration of the objects
  double mu=10.0;
  double mu2=0;
  double wmin=-4;
  double wmax=4;
  int Nbins=1000;
  collective ColNus(mu,wmin,wmax,Nbins);
  collective ColNus2(mu2,wmin,wmax,Nbins);
  std::cout << "Computing Collective" << std::endl;
  std::ofstream file("collective.dat");

  // Evolve and save the evolution
  // for(double t=0;t<tf;t+=dt){
  //  progressbar(100*t/tf);
  double dt=10.23;
  ColNus.EvolveSUN(dt);
  //ColNus2.EvolveSUN(dt);
  SU_vector o=ColNus.ez;
  for(int w=0;w<Nbins;w++){
    file << std::scientific << ColNus.Get_x(w) << "\t" 
  	 << ColNus.GetExpectationValue(o,0,w)/ 
     ColNus2.GetExpectationValue(o,0,w) << std::endl;//"  " 
    // 	 << ColNus.GetExpectationValue(ColNus.b1_proj[1],0,w) << std::endl;
  }
  
  file.close();


}
