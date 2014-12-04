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
  rabi R0,Rd;
  // de-tuning
  double del;

  // delta time for the prints
  double dt=0.01;
  // Final time
  double tf=1200;


  // Tuned Rabi system
  R0.init(10,10,0.1);

  // Setting the errors
  R0.Set_rel_error(1e-5);
  R0.Set_abs_error(1e-5);

  std::cout << "Rabi system with frequency of 10 initialized." << std::endl;
  std::cout << "give the value for the detuning: " << std::endl;
  //std::cin >> del;
  del=1;

  // un-tuned Rabi system
  Rd.init(10,10+del,0.1);
  // Setting the errors
  Rd.Set_rel_error(1e-5);
  Rd.Set_abs_error(1e-5);

  std::cout << "Computing rabi" << std::endl;
  std::ofstream file("rabi.dat");

  // Evolve and save the evolution
  for(double t=0;t<tf;t+=dt){
    progressbar(100*t/tf);
    R0.EvolveSUN(dt);
    file << t << "\t" << R0.GetExpectationValue(R0.d0,0,0) << "  " 
	 << R0.GetExpectationValue(R0.evol_b0_proj[0],0,0) << "  " 
	 << R0.GetExpectationValue(R0.evol_b0_proj[1],0,0) << std::endl;
  }
  file.close();
  file.open("rabi_detuned.dat");
  std::cout << std::endl << "Computing detuned rabi" << std::endl;
  for(double t=0;t<tf;t+=dt){
    progressbar(100*t/tf);
    Rd.EvolveSUN(dt);
    file << t << "\t" << Rd.GetExpectationValue(Rd.d0,0,0) << "  " 
	 << Rd.GetExpectationValue(Rd.b0_proj[0],0,0) << "  " 
	 << Rd.GetExpectationValue(Rd.b0_proj[1],0,0) << std::endl;
  }
  file.close();
  
  //Ask whether to run the gnuplot script
  std::string plt;
  std::cout << std::endl <<  "Done! " << std::endl <<
   "  Do you want to run the gnuplot script? yes/no" << std::endl;
  //std::cin >> plt;
  if(plt=="yes" || plt=="y")
    return system("./plot.plt");
  return 0;
}
