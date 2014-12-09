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
 *      Christopher Weaver (University of Wisconsin Madison)                   *
 *         chris.weaver@icecube.wisc.edu                                       *
 ******************************************************************************/


/*******************************************************************************
 * Example of the collective phenomena for neturinos,                          *
 * based on Georg G. Raffelt  Phys.Rev. D83 (2011) 105022.                     *
 * arXiv:1103.2891                                                             *
 ******************************************************************************/

#include "collective.h"

void collective::progressbar( int percent, double mu){
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
  std::cout<< percent << "%   mu: " << mu << std::flush;
}

void collective::init(double m,double th, double wmin, double wmax, int Nbins){
  mu=m;
  ini(Nbins,2,1,0,0.0);
  w_min=wmin;
  w_max=wmax;
  theta=th;

  Set_xrange(wmin,wmax,"lin");

  Set_CoherentInteractions(true);  

  ex=SU_vector::Generator(nsun,1);
  ey=SU_vector::Generator(nsun,2);
  ez=SU_vector::Generator(nsun,3);

  B=ez;

  // set initial conditions for the density matrix.
  double Norm=0;

  for(int ei = 0; ei < nx; ei++){
    double w=Get_x(ei);
    if(w>0) Norm+=(1.0/(w*w)*Fermi(1/(2*w)));
  }
  Norm=Norm*(w_max-w_min)/((double)nx);

  for(int ei = 0; ei < nx; ei++){
    double w=Get_x(ei);
      if(w>0){
        state[ei].rho[0] = (ey*sin(theta) + ez*cos(theta))*(1.0/(Norm*w*w)*Fermi(1/(2*w)));
      }else{
        state[ei].rho[0] = (ey*sin(theta) + ez*cos(theta))*(-0.7/(Norm*w*w)*Fermi(-1/(2*w)));
      }
    }

  //setting errors and step function for the GSL
  Set_rel_error(1e-8);
  Set_abs_error(1e-8);
  Set_h(1e-10);
  Set_GSL_step(gsl_odeiv2_step_rk4);
}


void collective::PreDerive(double t){
  P=state[0].rho[0];
  for(int ei = 1; ei < nx; ei++){
    P+=state[ei].rho[0];
  }
}

SU_vector collective::HI(int ix,double t){
  mu = mu_f+(mu_i-mu_f)*(1.0-t/period);
  if(bar)
    progressbar(100*t/period, mu);
  return Get_x(ix)*B+P*(mu*(w_max-w_min)/(double)nx);
}

void collective::Adiabatic_mu(double mui, double muf, double per, bool b){
  period=per;
  bar=b;
  mu_i=mui;
  mu_f=muf;
  EvolveSUN(period);
  if(bar)
    std::cout << std::endl;  
}

double collective::Fermi(double EoverT){
  return 1.0/(exp(EoverT)+1.0);
}
