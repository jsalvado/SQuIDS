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


#ifndef SQUIDS_SUNALG_H
#define SQUIDS_SUNALG_H

#if __cplusplus < 201103L
#error C++11 compiler required. Update your compiler and use the flag -std=c++11
#endif

#include <cassert>
#include <vector>
#include <iosfwd>
#include <cmath>
#include <memory>
#include <stdexcept>

#include <gsl/gsl_complex.h>

#include "const.h"
#include "SU_inc/dimension.h"
#include "detail/ProxyFwd.h"

namespace squids{

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
/// allocations and deallocations by rewriting in this form:
///
///     SU_vector temp1, temp2;
///     for(unsigned int i=0; i<N; i++){
///       temp1 = v1[i].SUEvolve(v2[i],t);
///       temp2 = v3[i].SUEvolve(v4[i],t)
///       state[i] += temp1 + temp2;
///     }
///
/// 2. If a calculation has an SU_vector calulation as a subexpression, but
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
      return(WrapperType::apply(*this,static_cast<SU_vector>(proxy))); //evaluate via a temporary
    //check whether sizes match
    if(this->size!=proxy.suv1.size){
      if(isinit_d) //can't resize
        throw std::runtime_error("Non-matching dimensions in assignment to SU_vector with external storage");
      if(!WrapperType::allowTargetResize)
        throw std::runtime_error("Non-matching dimensions in SU_vector assignment");
      //can resize
      if(isinit)
        delete[] components;
      dim=proxy.suv1.dim;
      size=proxy.suv1.size;
      if(proxy.mayStealArg1()){ //if the operation is component-wise and suv1 is an rvalue
        components=proxy.suv1.components; //take suv1's backing storage
        isinit=proxy.suv1.isinit;
        isinit_d=proxy.suv1.isinit_d;
        if(isinit)
          const_cast<SU_vector&>(proxy.suv1).isinit=false; //complete the theft
      }
      else{
        components=new double[size];
        isinit=true;
      }
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
  
  ///\brief construct a vector from the result of a vector arithmetic expression
  ///
  /// The newly constructed vector will take the data buffer from one of the
  /// operands, if possible, to elide memory allocation.
  template<typename ProxyType>
  SU_vector(ProxyType&& proxy, REQUIRE_EVALUATION_PROXY_FPARAM):
  dim(proxy.suv1.dim),
  size(proxy.suv1.size),
  components(proxy.mayStealArg1() ? proxy.suv1.components : new double[size]),
  isinit(proxy.mayStealArg1() ? proxy.suv1.isinit : true),
  isinit_d(proxy.mayStealArg1() ? proxy.suv1.isinit_d : false)
  {
    if(components==proxy.suv1.components && proxy.suv1.isinit)
      const_cast<SU_vector&>(proxy.suv1).isinit=false; //complete the theft
    proxy.compute(detail::vector_wrapper<detail::AssignWrapper>{dim,components});
  }

  ///\brief Construct an SU_vector from a GSL matrix
  ///
  ///\param m The matrix whose data should be copied
  ///\pre m is hermitian
  ///\pre m->size1 <= 6
  SU_vector(const gsl_matrix_complex* m);
  
  ///\brief Construct an SU_vector from existing data
  ///
  /// The newly constructed SU_vector will allocate its own storage, but it
  /// will copy its component information from data.
  ///\param data The existing vector components to use
  ///\pre data has a square size (4, 9, 16, 25, or 36)
  SU_vector(const std::vector<double>& data);

  // destructor
  ~SU_vector(){
    if(isinit)
      delete[] components;
  }

  //*************
  // Functions
  //*************
  
  ///\brief Overwrite all components with a single value
  ///
  ///\param v The value with which to fill
  void SetAllComponents(double v);
  
  ///\brief Get a copy of the SU_vector's components
  std::vector<double> GetComponents() const;
  
  ///\brief Returns a rotated SU_vector with a rotation in the ij-subspace
  ///\param i subspace index
  ///\param j subspace index
  ///\param theta rotation angle in radians
  ///\param delta complex phase
  SU_vector Rotate(unsigned int i, unsigned int j, double theta, double delta) const;
  ///\brief Same as RotateToB0, but with reversed angles
  ///\param params Contains the rotation parameters
  ///
  /// If the SU_vector is in the B0 basis, it transforms it to the B1 representation.
  void RotateToB1(const Const& params);
  ///\brief Rotates a SU_vector given the angles in Const
  ///\param params Contains the rotation parameters
  ///
  /// If the SU_vector is in the B1 basis, it transforms it to the B0 representation.
  void RotateToB0(const Const& params);

  ///\brief Gets the dimension of the SU_vector
  unsigned int Dim() const {return dim;}

  ///\brief Gets the number of components in the vector
  unsigned int Size() const { return size; }

  ///\brief Applies unitary transformation of the form
  /// Exp(-Op)vExp(Op) where Op is represented by v.
  SU_vector UTransform(const SU_vector& v) const;

  ///\brief Returns the Eigen vector set that diagonalizes
  /// the given SU_vector
  ///\brief Construct a GSL matrix from a SU_vector
  std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> GetGSLMatrix() const;

  ///\brief Compute the time evolution of the SU_vector
  ///
  /// Evolves the SU_vector with the evolution operator given by op.
  ///\param op The evolution operator
  ///\param time The time over which to do the evolution
  ///\pre op must be diagonal
  ///\returns An object convertible to an SU_vector
  detail::EvolutionProxy Evolve(const SU_vector& op, double time) const;


  //**********
  //operators
  //**********
  ///\brief Equality comparison
  bool operator ==(const SU_vector&) const;

  ///\brief Scalar product of two SU_vectors.
  ///
  /// Equivalent to the trace of the matrix multiplication.
  double operator*(const SU_vector&) const;
  
  ///\brief Multiplication by a scalar
  ///\returns An object convertible to an SU_vector
  detail::MultiplicationProxy operator*(const double) const &;
  
  ///\brief Multiplication by a scalar
  ///\returns An object convertible to an SU_vector
  detail::MultiplicationProxy operator*(const double) &&;
  
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
  
  ///\brief Multiplying assignment
  SU_vector & operator*=(double);
  
  ///\brief Dividing assignment
  SU_vector & operator/=(double);
  
  ///\brief Addition
  ///\returns An object convertible to an SU_vector
  detail::AdditionProxy operator +(const SU_vector&) const &;
  
  ///\brief Addition
  ///\returns An object convertible to an SU_vector
  detail::AdditionProxy operator +(SU_vector&&) const &;
  
  ///\brief Addition
  ///\returns An object convertible to an SU_vector
  detail::AdditionProxy operator +(const SU_vector&) &&;
  
  ///\brief Addition
  ///\returns An object convertible to an SU_vector
  detail::AdditionProxy operator +(SU_vector&&) &&;
  
  ///\brief Subtraction
  ///\returns An object convertible to an SU_vector
  detail::SubtractionProxy operator -(const SU_vector& other) const &;
  
  ///\brief Subtraction
  ///\returns An object convertible to an SU_vector
  detail::SubtractionProxy operator -(const SU_vector& other) &&;
  
  ///\brief Negation
  ///\returns An object convertible to an SU_vector
  detail::NegationProxy operator -() const &;
  
  ///\brief Negation
  ///\returns An object convertible to an SU_vector
  detail::NegationProxy operator -() &&;
  
  ///\brief Optimized assignment from the result of an arithmetic expression
  template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
  SU_vector& operator=(const ProxyType& proxy){
    return(assignProxy<detail::AssignWrapper>(proxy));
  }
  
  ///\brief Optimized incrementing assignment from the result of an
  /// arithmetic expression
  template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
  SU_vector& operator+=(const ProxyType& proxy){
    return(assignProxy<detail::IncrementWrapper>(proxy));
  }
  
  ///\brief Optimized decrementing assignment from the result of an
  /// arithmetic expression
  template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
  SU_vector& operator-=(const ProxyType& proxy){
    return(assignProxy<detail::DecrementWrapper>(proxy));
  }
  
  ///\brief Set the external storage used by this SU_vector
  void SetBackingStore(double* storage){
    if(isinit)
      delete[] components;
    components = storage;
    isinit=false;
    isinit_d=true;
  }

  ///\brief Array-like indexing
  ///\pre i < dimension^2
  double& operator[](unsigned int i) {assert(i < size); return components[i];}
  ///\brief Array-like indexing
  ///\pre i < dimension^2
  const double& operator[](unsigned int i) const {assert(i < size); return components[i];}

  ///\brief Construct a projection operator
  ///
  ///\param d The dimension of the operator
  ///\param i The dimension selected by the projection
  static SU_vector Projector(unsigned int d, unsigned int i);
  
  ///\brief Construct an identity operator
  ///
  ///\param d The dimension of the operator
  static SU_vector Identity(unsigned int d);
  
  ///\brief Constructs a projector to the upper subspace of dimension i
  ///
  ///\param d The dimension of the operator
  ///\param i The subspace dimension
  ///
  /// \f$ NegProj = diag(1,...,1,0,...,0) \f$ where the last one is at the \f$i-1\f$ entry.
  static SU_vector PosProjector(unsigned int d, unsigned int i);
  
  ///\brief Constructs a projector to the lower subspace of dimension i
  ///
  ///\param d The dimension of the operator
  ///\param i The subspace dimension
  ///
  /// \f$ NegProj = diag(0,...,0,1,...,1) \f$ where the first one is at the \f$d-i\f$ entry.
  static SU_vector NegProjector(unsigned int d, unsigned int i);
  
  ///\brief Creates a SU_vector corresponding to the \f$i\f$SU_N basis generator.
  //
  // Tts represented by \f$ v = (0,...,1,...,0)\f$ where 1 is in the \f$i\f$ component.
  static SU_vector Generator(unsigned int d, unsigned int i);
  
  template<typename Op>
  friend struct detail::EvaluationProxy;
  friend struct detail::EvolutionProxy;
  friend struct detail::AdditionProxy;
  friend struct detail::SubtractionProxy;
  friend struct detail::NegationProxy;
  friend struct detail::MultiplicationProxy;
  friend struct detail::iCommutatorProxy;
  friend struct detail::ACommutatorProxy;
  
  friend detail::iCommutatorProxy iCommutator(const SU_vector&,const SU_vector&);
  friend detail::ACommutatorProxy ACommutator(const SU_vector&,const SU_vector&);
  friend double SUTrace(const SU_vector&,const SU_vector&);
  
  //overloaded output operator
  friend std::ostream& operator<<(std::ostream&, const SU_vector&);
};

///\brief Returns the trace of the product of the SU_vector matrix representations
///
/// Defines a scalar product of SU_vectors.
double SUTrace(const SU_vector&,const SU_vector&);

///\brief Commutator of two SU_vectors
detail::iCommutatorProxy iCommutator(const SU_vector&,const SU_vector&);

///\brief Anticommutator of two SU_vectors
detail::ACommutatorProxy ACommutator(const SU_vector&,const SU_vector&);

///\brief Multiplication of an SU_vector by a scalar from the left.
detail::MultiplicationProxy operator*(double x, const SU_vector& v);

///\brief Multiplication of an SU_vector by a scalar from the left.
detail::MultiplicationProxy operator*(double x, SU_vector&& v);
  
} //namespace squids

#include "detail/ProxyImpl.h"

#endif
