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

/**
 \mainpage
 
 Introduction
 ------------
 
 Simple Quantum Integro-Differential Solver (SQuIDS) is a C++ code designed to 
 solve semi-analytically the evolution of a set of density matrices and scalar 
 functions. SQuIDS provides a base class from which users can derive new classes 
 to include new non-trivial terms from the right hand sides of density matrix 
 equations. The code was designed in the context of solving neutrino oscillation 
 problems, but can be applied in any problem that involves solving the quantum 
 evolution of a collection of particles with Hilbert space of dimension up to 6.
 
 Dependencies
 ------------
 
 SQuIDS is written using features of C++11 for efficiency, so it requires a
 reasonably new C++ compiler. GCC 4.7 or newer or Clang 3.3 or newer is 
 recommended, but the code is written to be portable and should work with 
 other compilers as well. 
 
 The GNU Scientific Library (GSL),
 version 1.14 or newer is needed by SQuIDS to perform differential equation 
 integration. The latest version can be found at 
 https://www.gnu.org/software/gsl/ .
 
 Finally, gnuplot (http://www.gnuplot.info) version 4.0 or newer is required 
 to generate plots from the example programs, but is not needed by the library 
 itself.
 
 Building
 --------
 
 First, run the included configure script: 
 
     ./config
 
 It may be useful to run it with the `--help` flag in order to find out all
 available options. 
 
 After configuring the build, simply compile with the generated makefile:
 
     make
 
 You can verify that the library is working correctly by running the unit 
 tests with `make test`, and you can install the built library with 
 `make install`.
 
 Usage
 -----
 
 To learn about using the library, see the examples in the `examples` 
 subdirectory. 
 
 In general, implementing a system using this library requires defining a
 new subclass of the SQUIDS base class which encapsulates the problem, and
 using the SU_vector class to represent quantum mechanical states and 
 operators.
 
 */

#ifndef SQUIDS_H
#define SQUIDS_H

#include "SUNalg.h"

#include <iosfwd>
#include <vector>
#include <memory>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

///\brief SQuIDS main class

//density matrix kinetic equation solver
class SQUIDS {
 protected:
  
  ///\brief Structure that contains the node state
  struct SU_state
  {
    ///\brief Vector of SU(N) vectors that represents the quantum part of the state
    std::unique_ptr<SU_vector[]> rho;
    ///\brief Vector of scalars that represents the classic part of the state
    double* scalar; //not owned
  };
  
  bool CoherentInt,NonCoherentInt,OtherInt,ScalarsInt,AnyNumerics;
  bool is_init;
  bool adaptive_step;
 
  std::unique_ptr<double[]> x;
  std::unique_ptr<double[]> delx;
  double t;
  double t_ini;

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

  std::unique_ptr<double[]> system;
  double* deriv_system; //not owned

  int Nsystem;

  Const params;

  std::unique_ptr<SU_state[]> state;

  int numeqn;

  // setting up GSL ODE solver
  gsl_odeiv2_step_type* step; //not owned

  gsl_odeiv2_system sys;

  double h;
  double h_min;
  double h_max;
  double abs_error;
  double rel_error;
    
  //***************************************************************
  ///\brief Sets the derivative system pointer for GSL use
  void set_deriv_system_pointer(double*);
  //interface function called by GSL
  friend int RHS(double ,const double*,double*,void*);

 public:
  ///\todo
  std::unique_ptr<SU_state[]> dstate;
  //****************
  //Constructors
  //****************
  SQUIDS();
  //***************************************************************
  ///\brief Constructs a SQUIDS object
  ///
  ///\param nx Number of components in the array "x"
  ///\param dim Dimension of SU(n)
  ///\param nrho Number of density matrix in every "x" site
  ///\param nscalar Number of scalars in every "x" site
  ///\param ti initial value for the evolution parameter t

  SQUIDS(int nx,int dim,int nrho,int nscalar, double ti=0.0);

  //***************************************************************
  virtual ~SQUIDS();

  //***************************************************************
  ///\brief Initializes a SQUIDS object
  ///
  ///\param nx Number of components in the array "x"
  ///\param dim Dimension of SU(n)
  ///\param nrho Number of density matrix in every "x" site
  ///\param nscalar Number of scalars in every "x" site
  ///\param ti initial value for the evolution parameter t
  void ini(int nx,int dim,int nrho, int nscalar, double ti=0.0);


  //***************************************************************
  ///\brief Set the range of values for the array "x"
  ///\param xini  x_min
  ///\param xend  x_max
  ///\param scale "log or lin" type of scale
  int Set_xrange(double xini, double xend, std::string scale);

  //***************************************************************
  ///\brief Returns the closes position in the array x for the value given
  ///\param x value of x to look for
  int Get_i(double x) const;

  //***************************************************************
  ///\brief Returns de value in the position "i"
  ///\param i node position
  double Get_x(int i) const{return x[i];};

  //***************************************************************
  //virtual functions defined in the dervied class
  ///\brief H0 time independent evolution operator
  virtual SU_vector H0(double x) const{ return SU_vector(nsun);}
  ///\brief H1 time dependent evolution operator
  virtual SU_vector HI(int ix,double t) const{ return SU_vector(nsun);}
  ///\brief Attenuation and/or decoherence operator
  ///\param ix Index in the x-array
  ///\param t time
  virtual SU_vector GammaRho(int ix, double t) const{ return SU_vector(nsun);}
  ///\brief Function containing other possible operations, like non linear terms in rho
  ///or terms involving the scalar functions
  ///\param ix Index in the x-array
  ///\param t time
  virtual SU_vector InteractionsRho(int ix,double t) const{ return SU_vector(nsun);}
  ///\brief Attenuation for the scalar functions
  ///\param ix Index in the x-array
  ///\param t time
  virtual double GammaScalar(int ix, double t) const{return 0.0;}
  ///\brief Other possible interaction terms for the scalar functions.
  ///\param ix Index in the x-array
  ///\param t time
  virtual double InteractionsScalar(int ix,double t) const{return 0.0;}
  ///\brief Function to be evaluated before the derivative
  ///\param t time
  ///
  /// This function enables the user to perform operations or updates before the derivative.
  virtual void PreDerive(double t){};

  //***************************************************************
  ///\brief Computes the right hand side of the kinetic equation.
  ///\param t time
  int Derive(double t);

  //***************************************************************
  ///\brief Numerical evolution of the state using GSL
  ///\param dt evolution time interval.
  int Evolve(double dt);

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
  ///\brief Sets the GSL stepper function (numerical method)
  ///\param opt GSL step function
  void Set_GSL_step(const gsl_odeiv2_step_type * opt);

  ///\brief Turns on and off adaptive runge-kutta stepping
  ///\param opt If true: uses adaptive stepping, else: it does not.
  void Set_AdaptiveStep(bool opt);
  ///\brief Activate coherent interaction
  void Set_CoherentInteractions(bool opt);
  ///\brief Activate noncoherent interaction
  void Set_NonCoherentInteractions(bool opt);
  ///\brief Activate other SU_vector interactions
  void Set_OtherInteractions(bool opt);
  ///\brief Activate other scalar interactions
  void Set_ScalarInteractions(bool opt);
  ///\brief Set the minimum runge-kutta step
  void Set_h_min(double opt);
  ///\brief Set the maximum runge-kutta step
  void Set_h_max(double opt);
  ///\brief Set the initial runge-kutta step
  void Set_h(double opt);
  ///\brief Set the numerical relative error
  void Set_rel_error(double opt);
  ///\brief Set the numerical absolute error
  void Set_abs_error(double opt);
  ///\brief Set the time scale
  void Set_units(double opt);
  ///\brief Set the number of steps when not using adaptive stepping
  void Set_NumSteps(int opt);
  
  //***************************************************************
  ///\brief Returns the expectation value for a given operator for a give state irho in a node ix.
  ///\param op operator
  ///\param irho index of rho
  ///\param ix index in the array "x"
  double GetExpectationValue(SU_vector op,int irho,int ix) const;

  //***************************************************************
  ///\brief Returns the expectation value for a given operator for the rho given by irho 
  /// and using linear interpolation in "x"
  ///\param op operator 
  ///\param irho index of rho
  ///\param x value of x
  double GetExpectationValueD(SU_vector op,int irho,double x) const;
  
  ///\brief Returns the current time of the system
  double Get_t() const{ return(t); }
};

#endif
