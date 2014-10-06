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



#ifndef __SQUIDS_H
#define __SQUIDS_H

//#define CalNeuOscSUN_DEBUG

#include "const.h"
#include "SUNalg.h"

#include <iostream>
#include <float.h>
#include <math.h>
#include <complex>
#include <vector>
#include <limits>
#include <memory>


#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_odeiv2.h>


struct SU_state
{
  std::unique_ptr<SU_vector[]> rho;
  double* scalar; //not owned
};



//density matrix kinetic equation solver
class SQUIDS {
 protected:
  bool CoherentInt,NonCoherentInt,OtherInt,ScalarsInt,AnyNumerics;
  bool is_init;
  bool adaptive_step;
 
  std::unique_ptr<double[]> x;
  std::unique_ptr<double[]> delx;
  double t;
  double t_ini;
  double t_end;

  double tunit;
  int nx;
  int nsun;
  int nrhos;
  int nscalars;

  int nsteps;

  int index_rho;
  int index_scalar;

  int size_rho;
  int size_state;

  bool neu_and_aneu;

  std::unique_ptr<double[]> system;
  double* deriv_system; //not owned

  int Nsystem;

  Const params;

  std::unique_ptr<SU_state[]> state;


  int numeqn;
  //SU_alg SU;


  // setting up GSL ODE solver
  gsl_odeiv2_step_type* step; //not owned

  gsl_odeiv2_system sys;

  double h;
  double h_min;
  double h_max;
  double abs_error;
  double rel_error;

    
  //***************************************************************
  //sets the deriv system pointer
  void set_deriv_system_pointer(double*);
  //interface function called by GSL
  friend int RHS(double ,const double*,double*,void*);

 public:
  std::unique_ptr<SU_state[]> dstate;
  //****************
  //Constructors
  //****************
  SQUIDS();
  //***************************************************************
  //int -> Number of components in the array "x"
  //int -> Dimension of SU(n)
  //int -> Number of density matrix in every "x" site
  //int -> Number of scalars in every "x" site
  SQUIDS(int,int ,int,int);

  //***************************************************************
  virtual ~SQUIDS();

  //***************************************************************
  //Initializer, the same argumens as the constructor
  void ini(int,int ,int, int );

  //***************************************************************
  //set the range of values for the array "x"
  //double -> x_min
  //double -> x_max
  //string -> "log or lin" type of scale
  int Set_xrange(double , double, string);

  //***************************************************************
  //returns the closes position in the array x for the value given
  int Get_i(double);

  //***************************************************************
  //returns de value in the position "i"
  double Get_x(int i){return x[i];};

  //***************************************************************
  //virtual functions defined in the dervied class
  //H0 not time dependend evolution operator.
  virtual SU_vector H0(double){SU_vector V(nsun); return V;}
  //H1 time deppenden evolution operator.
  virtual SU_vector HI(int ix,double t){SU_vector V(nsun); return V;}
  //Atenuation and/or decoerence operator
  virtual SU_vector GammaRho(int ix, double t){SU_vector V(nsun); return V;}
  //Other function containing other possible operations, like non linear terms in rho
  //or terms involving the scalar functions
  virtual SU_vector InteractionsRho(int ix,double t){SU_vector V(nsun); return V;}
  //Atenuation for the scalar functions
  virtual double GammaScalar(int ix, double t){return 0.0;}
  //Other possible interaction terms for the scalar fucntions.
  virtual double InteractionsScalar(int ix,double t){return 0.0;}

  virtual void PreDerive(double){};

  //***************************************************************
  //Computes all the derivatives.
  int Derive(double);

  //***************************************************************
  //GSL evolution of the state.
  //double -> initial time.
  //double -> final time.
  int EvolveSUN(double, double);

  //***************************************************************
  //functions to set parameters in the object.
  // string -> parameter
  // other -> value
  //the possible parameters are:
  // ******* bool *************
  // CoherentInteractions
  // NonCoherentInteractions 
  // OtherInt
  // ScalarInteractions
  // AntiNeutrinos
  // ******* double *********
  // t -> time
  // Units -> time unit
  // The parameters in the struc "const" are also available from this functions
  // ******* int ************
  // nx
  // nsun
  // nrhos
  // nscalars
  // The parameters in the struc "const" are also available from this functions
  void Set(string,bool);
  void Set(string,double);
  void Set(string,int);
  void Set(string,const gsl_odeiv2_step_type*);
  

  //***************************************************************
  //returns the expectation value for a given operator
  //SU_vector -> operator 
  //int -> index of rho
  //int -> index in the array "x"
  double GetExpectationValue(SU_vector&,int,int);

  //***************************************************************
  //returns the expectation value for a given operator, the same as the other,
  //but using the exact value for "x" for the H0 evolution (very usefull fur fast
  //oscilations in H0.
  //SU_vector -> operator 
  //int -> index of rho
  //double -> value of x
  double GetExpectationValueD(SU_vector&,int,double);


};

//Auxiliar function used for the GSL
int RHS(double ,const double*,double*,void*);

#endif
