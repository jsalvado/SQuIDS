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


#ifndef SQUIDS_SUNALG_H
#define SQUIDS_SUNALG_H

#include <cassert>
#include <vector>
#include <iosfwd>
#include <cmath>
#include <memory>
#include <stdexcept>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>

#include "const.h"
#include "detail/ProxyFwd.h"

//The code only works up to 6 dim Hilbert space
#define MAXSIZE 36

///\brief A vector represented in the SU(n) basis
///
/// An object in SU(n) has n^2 components. By default SU_vector will
/// automatically allocate sufficient 'backing storage' to contain
/// these. It is, however, possible to specify that an SU_vector
/// should treat some externally provided buffer as its backing
/// storage. In this case the size and lifetime of that buffer are
/// the resposibility of the user, so users are encouraged to avoid
/// using this mode unless it is required by their application, as
/// its use is more difficult and requires much greater care.
/// The external storage mode is primarily useful because it allows
/// interfacing with other, low-level numerical codes efficiently.
///
/// SU_vector provides overloaded mathematical operators so that algebra can
/// be written in a natural way.
///
/// The SQuIDs library has a limited ability to optimize away temporary objects.
/// That is, all operations of the forms
///     v1 [Op1]= v2 [Op2] v3;
///     v1 [Op1]= s * v2;
/// where v1, v2, and v3 are pre-existing objects of type SU_vector (and s is
/// a scalar) are performed without allocating memory. Op1 may +, -, or nothing
/// (normal assignment), and Op2 may be +, -, time evolution, a commutator or
/// an anticommutator.
///
/// This optimization is inhibited when v1 aliases v2 or v3 (they are the same
/// objects or they otherwise refer to the same backing storage).
/// This has no influence on the correctness of writing complex expressions in
/// terms of subexpressions: These will still be correctly evaluated, but memory
/// will be allocated for the results of the subexpressions, making the
/// calculation slower than if this can be avoided. It is expected that users
/// will write expressions in the form they find most natural, and only if
/// performance optimization is required consider restructuring code to take
/// deliberate advantage of this optimization. In that case, the following
/// techniques may be useful:
///
/// 1. If a calculation involving subexpressions is performed in a loop, it is
/// advantageous to manually create sufficient temporaries for all
/// subexpressions outside of the loop and then split the complex expression
/// into a series of basic operations whose results are stored to the
/// temporaries. This ensures that allocation will be performed only once per
/// temporary before entering the loop, rather than once per temporary, per
/// loop iteration. For example:
///
///     //assuming size N arrays of SU_vector state, v1, v2, v3, and v4,
///     //and a floating-point t
///     for(unsigned int i=0; i<N; i++)
///         state[i] += v1[i].SUEvolve(v2[i],t) + v3[i].SUEvolve(v4[i],t);
///
///    In this code the addition on the right hand side can be performed without
/// allocation, but each of the evolution operations must allocate a temporary,
/// so 2*N allocations and deallocations must occur. This can be reduced to 2
/// allocations and deallocation by rewriting in this form:
///
///     SU_vector temp1, temp2;
///     for(unsigned int i=0; i<N; i++){
///       temp1 = v1[i].SUEvolve(v2[i],t);
///       temp2 = v3[i].SUEvolve(v4[i],t)
///       state[i] += temp1 + temp2;
///     }
/// 2. If a calculation has an SU_vector calulation as a subexpression but
/// otherwise operates on scalars it can be useful to rewite the expression so
/// that the vector calculation forms the top level if possible:
///
///     //assuming SU_vectors v1 and v2 and scalars s1 and s2
///     v1 = s1*(s2*v2);
///
///    can be better reassociated as:
///
///     v1 = (s1*s2)*v2;
///
class SU_vector{
private:
  unsigned int dim;
  unsigned int size; //dim^2
  double *components;
  bool isinit; //this vector has been initialized and owns its storage
  bool isinit_d; //this vector has been initialized and does not own its storage
  
  ///\brief Internal implementation of optimized assignment
  template<typename WrapperType, typename ProxyType>
  SU_vector& assignProxy(const ProxyType& proxy){
    if(components==proxy.suv1.components || components==proxy.suv2.components) //beware of aliasing
      return(WrapperType::apply(*this,(SU_vector)proxy)); //evaluate via a temporary
    //check whether sizes match
    if(this->size!=proxy.suv1.size){
      if(isinit_d) //can't resize
        throw std::runtime_error("Non-matching dimensions in assignment to SU_vector with external storage");
      if(!WrapperType::allowTargetResize)
        throw std::runtime_error("Non-matching dimensions in SU_vector increment/decrement");
      //can resize
      if(isinit)
        delete[] components;
      dim=proxy.suv1.dim;
      size=proxy.suv1.size;
      components=new double[size];
      isinit=true;
    }
    //evaluate in place
    proxy.compute(detail::vector_wrapper<WrapperType>{dim,components});
    return(*this);
  }
  
public:
  //***************
  // Constructors 
  //***************
  
  ///\brief Default constructor
  ///
  /// Constructs an empty SU_vector with no size.
  SU_vector();
  
  ///\brief Copy constructor
  ///
  /// The newly constructed SU_vector will allocate its own storage
  /// which it will manage automatically.
  SU_vector( const SU_vector& V);
  
  ///\brief Move constructor
  ///
  /// If V owned its own storage it will be taken by the newly
  /// constructed SU_vector, and V will be left empty, as if default
  /// constructed.
  SU_vector( SU_vector&& V);
  
  ///\brief Construct a vector with a particular dimension
  ///
  /// Creates a new SU_vector with the specified dimension and initializes it
  /// to zero.
  ///
  ///\param dim The dimension of the vector
  ///\pre dim <= 6
  explicit SU_vector(unsigned int dim);

  ///\brief Construct a vector with external storage
  ///
  /// The newly constructed SU_vector will treat the specified data buffer as
  /// its backing storage. The user is responsible for ensuring that this buffer
  /// is large enough (at least dim*dim), and has a lifetime at least as long
  /// as the constructed vector and any vectors which inherit the use of this
  /// buffer from it by move construction or assignment. The contents of the
  /// buffer will not be modified during construction.
  ///
  ///\param dim The dimension of the SU_vector
  ///\param storage The data buffer which the SU_vector will use
  ///\pre dim <= 6
  SU_vector(unsigned int dim, double* storage);

  ///\brief Construct an SU_vector from a GSL matrix
  ///
  ///\param m The matrix whose data should be copied
  ///\pre m is hermitian
  ///\pre m->size1 <= 6
  SU_vector(gsl_matrix_complex* m);
  
  ///\brief Construct an SU_vector from existing data
  ///
  /// The newly constructed SU_vector will allocate its own storage, but it
  /// will copy its component information from data.
  ///\param data The existing vector components to use
  ///\pre data has a square size (4, 9, 16, 25, or 36)
  SU_vector(std::vector<double> data);

  // destructor
  ~SU_vector();

  //*************
  // Functions
  //*************
  
  ///\brief Overwrite all components with a single value
  ///
  ///\param v The value with which to fill
  void SetAllComponents(double v);
  
  ///\brief Get a copy of the SU_vector's components
  std::vector<double> GetComponents() const;
  
  ///brief Scale components by a constant
  ///
  ///\param a The scaling factor to use
  SU_vector Rescale(double a);
  
  //Change of bases with a rotation on the components given by the integers,
  //double -> angle
  //double -> phase
  ///\todo Needs proper documentation
  SU_vector Rotate(unsigned int, unsigned int, double, double);
  //Const, contines a standar set of angles, this funcion does the set of rotations
  //to change the bases.
  ///\todo Needs proper documentation
  void RotateToB1(const Const&);
  ///\todo Needs proper documentation
  void RotateToB0(const Const&);

  //equivalent in Matrix notation to the trace of the product for the two operators
  ///\todo Needs proper documentation
  double SUTrace(SU_vector*,SU_vector*);
  ///\todo Needs proper documentation
  double SUTrace(SU_vector&,const SU_vector&);

  ///\brief Gets the dimension of the SU_vector
  unsigned int Dim() const {return dim;};

  ///\brief Compute the time evolution of the SU_vector
  ///
  /// Evolves the SU_vector with the evolution operator given by op.
  ///\param op The evolution operator
  ///\param time The time over which to do the evolution
  ///\pre op must be diagonal
  ///\returns An object convertible to an SU_vector
  detail::EvolutionProxy SUEvolve(const SU_vector& op, double time) const;

  //**********
  //operators
  //**********
  ///\brief Equality comparison
  bool operator ==(const SU_vector&);

  ///\brief Scalar product of two SU_vectors.
  ///
  /// Equivalent to the trace of the matrix multiplication.
  double operator*(const SU_vector&);
  
  ///\brief Multiplication by a scalar
  ///\returns An object convertible to an SU_vector
  detail::MultiplicationProxy operator*(const double) const;
  
  ///\brief Assignment
  ///
  /// Assignment will fail if this vector uses external storage and the
  /// dimensions of the two vectors differ.
  SU_vector & operator =(const SU_vector&);
  
  ///\brief Move assignement
  ///
  /// If this vector is empty or owns its own storage it will switch to
  /// using whatever storage was used by other, causing other to relinquish
  /// any ownership of that storage.
  /// If, however, this vector is non-empty and uses external storage, it
  /// will copy other's data rather than shifting its storage. In this case
  /// if the dimensions of the two vectors differ the assignment will fail.
  SU_vector & operator =(SU_vector&& other);
  
  ///\brief Incrementing assignment
  SU_vector & operator+=(const SU_vector&);
  
  ///\brief Decrementing assignment
  SU_vector & operator-=(const SU_vector&);
  
  ///\brief Addition
  ///\returns An object convertible to an SU_vector
  detail::AdditionProxy operator +(const SU_vector&) const;
  
  ///\brief Subtraction
  ///\returns An object convertible to an SU_vector
  detail::SubtractionProxy operator -(const SU_vector& other) const;
  
  ///\brief Optimized assignment from the result of an arithmetic expression
  template<typename ProxyType, REQUIRE_EVALUATION_PROXY>
  SU_vector& operator=(const ProxyType& proxy){
    return(assignProxy<detail::AssignWrapper>(proxy));
  }
  
  ///\brief Optimized incrementing assignment from the result of an
  /// arithmetic expression
  template<typename ProxyType, REQUIRE_EVALUATION_PROXY>
  SU_vector& operator+=(const ProxyType& proxy){
    return(assignProxy<detail::IncrementWrapper>(proxy));
  }
  
  ///\brief Optimized decrementing assignment from the result of an
  /// arithmetic expression
  template<typename ProxyType, REQUIRE_EVALUATION_PROXY>
  SU_vector& operator-=(const ProxyType& proxy){
    return(assignProxy<detail::DecrementWrapper>(proxy));
  }
  
  ///\brief Set the external storage used by this SU_vector
  void SetBackingStore(double* storage);

  ///\brief Array-like indexing
  ///\pre i < dimension^2
  double& operator[](unsigned int i) {assert(i < dim*dim); return components[i];};
  ///\brief Array-like indexing
  ///\pre i < dimension^2
  const double& operator[](unsigned int i) const {assert(i < dim*dim); return components[i];}

  ///\brief Construct a projection operator
  ///
  ///\param d The dimension of the operator
  ///\param i The dimension selected by the projection
  static SU_vector Projector(unsigned int d, unsigned int i);
  
  ///\brief Construct an identity operator
  ///
  ///\param d The dimension of the operator
  static SU_vector Identity(unsigned int d);
  
  ///\todo Needs proper documentation
  static SU_vector PosProjector(unsigned int d, unsigned int i);
  
  ///\todo Needs proper documentation
  static SU_vector NegProjector(unsigned int d, unsigned int i);
  
  ///\todo Needs proper documentation
  static SU_vector Component(unsigned int d, unsigned int i);
  
  template<typename Op>
  friend struct detail::EvaluationProxy;
  friend struct detail::EvolutionProxy;
  friend struct detail::AdditionProxy;
  friend struct detail::SubtractionProxy;
  friend struct detail::MultiplicationProxy;
  friend struct detail::iCommutatorProxy;
  friend struct detail::ACommutatorProxy;
  
  friend detail::iCommutatorProxy iCommutator(const SU_vector&,const SU_vector&);
  friend detail::ACommutatorProxy ACommutator(const SU_vector&,const SU_vector&);
  
  //ostream overload operator
  friend ostream& operator<<(ostream&, const SU_vector&);
};

//Commutator of two SU_vectors
detail::iCommutatorProxy iCommutator(const SU_vector&,const SU_vector&);

//Anticommutator of two SU_vectors
detail::ACommutatorProxy ACommutator(const SU_vector&,const SU_vector&);

//Multiplication of an SU_vectors by a scalar from the left.
detail::MultiplicationProxy operator*(double x, const SU_vector& v);

#include "detail/ProxyImpl.h"

#endif
