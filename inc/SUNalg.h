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


#ifndef __SUNALG_H
#define __SUNALG_H

#include <assert.h>
#include <vector>
#include <cstring>
#include <iostream>
#include <float.h>
#include <math.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>

#include "const.h"

//The code only works up to 6 dim Hilber space 
#define MAXSIZE 36

class SU_vector{
  friend class SU_alg;
private:
  int dim;
  int size;
  double *components;
  bool isinit;
  bool isinit_d;
public:
  //***************
  // Constructors 
  //***************
  
  //void and copy constructors
  SU_vector();
  SU_vector( const SU_vector& V);

  //reference constructor, the memory for the vector should be already 
  //reserved and point to the double pointer.
  SU_vector(int,double*);

  //Different constructors
  SU_vector(int);
  SU_vector(gsl_matrix_complex*);
  SU_vector(std::vector<double>);

  //standard constructors the options are:
  //string -> (Projector,Identity,PosProjector,NegProjector,Component)
  //another short names are also valid.
  //int -> option, for example for projector is the component
  //int -> dimmension of the algebra ex: for SU(4) will be 4
  SU_vector(string,int,int);
  // destructor
  ~SU_vector();
  // initializers, the arguments are the same as in the constructors.
  void InitSU_vector(int);
  void InitSU_vector(const SU_vector& V);
  void InitSU_vector(int,double*);
  void InitSU_vector(gsl_matrix_complex*);
  void InitSU_vector(std::vector<double>);
  void InitSU_vector(string,int,int);

  //*************
  // Functions
  //*************
  //multiplication for double.
  SU_vector Rescale(double);
  //Change of bases with a rotation on the components given by the integers,
  //double -> angle
  //double -> phase
  SU_vector Rotate(int,int,double, double);
  //Const, contines a standar set of angles, this funcion does the set of rotations
  //to change the bases.
  void RotateToB1(const Const*);
  void RotateToB0(const Const*);

  //equivalent in Matrix notation to the trace of the product for the two operators
  double SUTrace(SU_vector*,SU_vector*);
  double SUTrace(SU_vector&,const SU_vector&);

  //returns the dimension
  int Dim(void) const {return dim;};

  //Prints the SU_vecto
  void SUPrint (void) const;

  //Evolves the SU_vector with the evolution operator given as an argument.
  //The vector that gives the evolution should be diagonal (linear combinations of the
  //weight of the SU(N) algebra).
  SU_vector SUEvolve(SU_vector&,double);


  //**********
  //operators
  //**********
  bool operator ==(const SU_vector&);

  //Scalar product of two SU_vectors.
  //Is equivalent to the trace of the matrix multiplication.
  double operator*(const SU_vector&);
  SU_vector operator*(const double);
  SU_vector & operator =(const SU_vector&);
  SU_vector & operator+=(const SU_vector&);
  SU_vector & operator-=(const SU_vector&);
  SU_vector operator +(const SU_vector&);
  SU_vector operator -(const SU_vector&);


  // vector like functions
  double& operator[](int i) {assert(i < dim*dim); return components[i];};
  const double& operator[](int i) const {assert(i < dim*dim); return const_cast<SU_vector&>(*this)[i];};

  //ostream overload operator
  friend ostream& operator<<(ostream&, const SU_vector&);  


};

//**************************************************************
//SU(n) algebra object, it contains the functions to compute the 
//the algebra operations, i*commutator and anticommutator
//**************************************************************

class SU_alg{
private:
  int dim;
public:
  SU_alg(void){};
  SU_alg(int d){dim=d;};
  void init(int d){dim=d;}
  SU_vector iCommutator(const SU_vector&,const SU_vector&);
  SU_vector ACommutator(const SU_vector&,const SU_vector&);
};


//Left and right handed multiplication operators for SU_vectors.
SU_vector operator*(const double x, const SU_vector &other);
SU_vector operator*(const SU_vector &other, const double x);


#endif
