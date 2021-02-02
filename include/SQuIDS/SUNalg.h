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

//work around strange interaction between pgi and libc++
#if defined(__PGI)
  #define _LIBCPP_ASSERT
#endif

#include <cassert>
#include <cmath>
#include <functional>
#include <iosfwd>
#include <memory>
#include <stdexcept>
#include <vector>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "const.h"
#include "SU_inc/dimension.h"
#include "detail/ProxyFwd.h"
#include "detail/MatrixExp.h"
#if SQUIDS_USE_STORAGE_CACHE
  #include "detail/Cache.h"
#endif

namespace squids{

///\brief Returns the trace of the product of the SU_vector matrix representations
///
/// Defines a scalar product of SU_vectors.
template<unsigned int Flags=0>
double SUTrace(const SU_vector&,const SU_vector&);

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
/// Concerning alignment: Some SU_vector operations are more efficient when the
/// destination and operands' underlying storage has a particular alignment.
/// Specifically, for best performance even dimensioned SU_vectors should be 32
/// byte aligned, while odd dimensioned SU_vectors should have the second
/// component 32 byte aligned. SU_vector::make_aligned will automatically create
/// SU_vectors satisfying these criteria, as will the other static factory
/// functions.
class SU_vector{
private:
  unsigned int dim;
  unsigned int size; //dim^2
  double* components;
  unsigned char ptr_offset; //offset between components and the beginning of the allocation which must be released
  bool isinit; //this vector has been initialized and owns its storage
  bool isinit_d; //this vector has been initialized and does not own its storage

  ///\brief Internal implementation of optimized assignment
  template<typename WrapperType, typename ProxyType>
  SQUIDS_ALWAYS_INLINE SU_vector& assignProxy(const ProxyType& proxy){
    using traits=detail::operation_traits<ProxyType>;
    
    if(!traits::elementwise && !traits::no_alias_target &&
       (components==proxy.suv1.components ||
        (traits::vector_arity==2 && components==proxy.suv2.components))) //beware of aliasing
      return(WrapperType::apply(*this,static_cast<SU_vector>(proxy))); //evaluate via a temporary
    //check whether sizes match
    if(!traits::equal_target_size && this->size!=proxy.suv1.size){
      if(isinit_d) //can't resize
        throw std::runtime_error("Non-matching dimensions in assignment to SU_vector with external storage");
      if(!WrapperType::allowTargetResize)
        throw std::runtime_error("Non-matching dimensions in SU_vector assignment");
      //can resize
      if(isinit)
        deallocate_mem();
      dim=proxy.suv1.dim;
      size=proxy.suv1.size;
      if(proxy.mayStealArg1()){ //if the operation is component-wise and suv1 is an rvalue
        components=proxy.suv1.components; //take suv1's backing storage
        ptr_offset=proxy.suv1.ptr_offset;
        isinit=proxy.suv1.isinit;
        isinit_d=proxy.suv1.isinit_d;
        if(isinit)
          const_cast<SU_vector&>(proxy.suv1).isinit=false; //complete the theft
      }
      else if(proxy.mayStealArg2()){ //if the operation is component-wise and suv2 is an rvalue
        components=proxy.suv2.components; //take suv2's backing storage
        ptr_offset=proxy.suv2.ptr_offset;
        isinit=proxy.suv2.isinit;
        isinit_d=proxy.suv2.isinit_d;
        if(isinit)
          const_cast<SU_vector&>(proxy.suv2).isinit=false; //complete the theft
      }
      else{
        alloc_aligned(dim,size,components,ptr_offset);
        isinit=true;
      }
    }
    //evaluate in place
    proxy.compute(detail::vector_wrapper<WrapperType>{dim,components});
    return(*this);
  }
  
#if SQUIDS_USE_STORAGE_CACHE
  struct mem_cache_entry{
    double* storage;
    unsigned char offset;
    mem_cache_entry():storage(nullptr){}
    mem_cache_entry(double* p, uint8_t o):storage(p),offset(o){}
  };
  ///A cache of previously used backing storage blocks
  static
  #ifdef SQUIDS_THREAD_LOCAL
  SQUIDS_THREAD_LOCAL //one cache per thread if supported
  #endif
  detail::cache<mem_cache_entry,32> storage_cache[SQUIDS_MAX_HILBERT_DIM+1];
#endif
  
  ///A helper function which tries to put a memory block into the cache rather
  ///than deleting it.
  void deallocate_mem(){
#if SQUIDS_USE_STORAGE_CACHE
    bool cached=false;
    if(((intptr_t)(components+dim%2))%32 == 0) //only try to save aligned storage
      cached=storage_cache[dim].insert(mem_cache_entry{components,ptr_offset});
    if(!cached)
#endif
      delete[] (components-ptr_offset);
  }
  
  ///A helper function which fetches a memory block from the cache if possible,
  ///and otherwise allocates a new block with optimal alignment
  ///\param dim The dimension of the vector for which the storage is to be allocated
  ///\param size dim squared
  ///\param components The pointer to be set to the aligned, allocated memory
  ///\param ptr_offset The location where the necessary alignment offset should be stored
  static void alloc_aligned(unsigned int dim, unsigned int size,
                            double*& components, unsigned char& ptr_offset){
#if SQUIDS_USE_STORAGE_CACHE
    mem_cache_entry cache_result=storage_cache[dim].get();
    if(cache_result.storage){
      components=cache_result.storage;
      ptr_offset=cache_result.offset;
    }
    else{
#endif
      size_t before=size%2;
      size_t maxHeadroom=(32/sizeof(double))-1/*+before*/;
      components=new double[size+maxHeadroom];
      ptr_offset=(intptr_t)(components+before)%32; //bytes
      if(ptr_offset){
        ptr_offset=(32-ptr_offset)/sizeof(double); //convert to units of doubles
        assert(ptr_offset<=maxHeadroom);
        components+=ptr_offset;
      }
      assert((intptr_t)(components+before)%32 == 0);
#if SQUIDS_USE_STORAGE_CACHE
    }
#endif
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
  isinit(proxy.mayStealArg1() ? proxy.suv1.isinit : true),
  isinit_d(proxy.mayStealArg1() ? proxy.suv1.isinit_d : false)
  {
    if(proxy.mayStealArg1()){
      components=proxy.suv1.components;
      ptr_offset=proxy.suv1.ptr_offset;
    }
    else
      alloc_aligned(dim,size,components,ptr_offset);
    
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

  ///\brief Construct an SU_vector from a GSL matrix
  ///
  ///\param m The matrix whose data should be copied
  ///\pre m is hermitian
  ///\pre m->size1 <= 6
  SU_vector(std::unique_ptr<gsl_matrix_complex, void (*)(gsl_matrix_complex *)> &&m):SU_vector(m.get()){};

  ///\brief Construct an SU_vector from existing data
  ///
  /// The newly constructed SU_vector will allocate its own storage, but it
  /// will copy its component information from data.
  ///\param data The existing vector components to use
  ///\pre data has a square size (4, 9, 16, 25, or 36)
  SU_vector(const std::vector<double>& data);

  ///\brief Constructs a vector with a particular dimension, with optimal
  ///       alignment.
  static SU_vector make_aligned(unsigned int dim, bool zero_fill=true);
  
  // destructor
  ~SU_vector(){
    if(isinit)
      deallocate_mem();
  }

  //*************
  // Functions
  //*************

  ///\brief Overwrite all components with a single value
  ///
  ///\param v The value with which to fill
  void SetAllComponents(double v){
    for(unsigned int i=0; i<size; i++)
      components[i] = v;
  }

  ///\brief Get a copy of the SU_vector's components
  std::vector<double> GetComponents() const;

  ///\brief Returns a rotated SU_vector rotating by the matrix m.
  ///\param rotation_matrix gsl_complex_matrix which rotates to
  SU_vector Rotate(const gsl_matrix_complex* rotation_matrix) const;

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

  ///\brief It does the tranformation given by the concatenation of two rotations weighted with a diagonal matrics Yd
  ///\brief V^{\dagger}Yd W (SU_vector) W^{dagger}Yd V
  ///\param paramsV Mixings for the unitary tranformation V
  ///\param paramsYd Diagonal matrix with the weights (yukawas like)
  ///\param paramsW Mixings for the unitary tranformation W
  void WeightedRotation(const Const& paramV, const SU_vector& Yd, const Const& paramW);  
  void WeightedRotation(gsl_matrix_complex* V, const SU_vector& Yd, gsl_matrix_complex* W);

  void Transpose(void);

  ///\brief It does the transformation equivalent to transpose the SU_vector in the matrix representation
  void Tr(void);


  ///\brief Applies unitary transformation given by the complex matrix em
  /// em*v*em^\dagger where Op is represented by v.
  SU_vector UTransform(gsl_matrix_complex* em) const;
  ///\brief Applies unitary transformation given by the complex matrix em
  /// em^\dagger*v*em where Op is represented by v.
  SU_vector UDaggerTransform(gsl_matrix_complex* em) const;

  ///\brief It returns as a SU_vector the complex part of the corresponding complex matrix represantation
  SU_vector Imag(void) const;
  ///\brief It returns as a SU_vectorthe real part of the corresponding complex matrix represantation
  SU_vector Real(void) const;


  ///\brief Applies unitary transformation of the form
  /// Exp(-scale*Op)*(this)*Exp(scale*Op) where Op is represented by v.
  SU_vector UTransform(const SU_vector& v, gsl_complex scale = GSL_COMPLEX_ONE) const;



  ///\brief Returns the the eigen values and eigen vectors
  /// unitary transformation that diagonalizes the matrix
  /// represented by the SU_vector
  std::pair<std::unique_ptr<gsl_vector,void (*)(gsl_vector*)>,
  std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)>> GetEigenSystem(bool order = true) const;

  ///\brief Construct a GSL matrix from a SU_vector
  std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> GetGSLMatrix() const;
  
  ///\brief Construct a GSL matrix from a SU_vector
  ///\param m the matrix into which the result is stored
  ///\pre m must be a square matrix of size this->Dim()
  void GetGSLMatrix(gsl_matrix_complex* m) const;

  ///\brief Compute the time evolution of the SU_vector
  ///
  /// Evolves the SU_vector with the evolution operator given by op.
  ///\param op The evolution operator
  ///\param time The time over which to do the evolution
  ///\pre op must be diagonal
  ///\returns An object convertible to an SU_vector
  detail::EvolutionProxy Evolve(const SU_vector& op, double time) const{
    return(detail::EvolutionProxy{op,*this,time});
  }
  
  ///\brief Get the buffer size required by PrepareEvolve
  size_t GetEvolveBufferSize() const{
    return(dim*(dim-1));
  }
  
  ///\brief Precompute operator dependent elements of an evolution
  ///
  /// Much of the calculation needed by Evolve depends only on the evolution
  /// operator and the time. If several SU_vectors will be evolved over the same
  /// time by the same operator, this work can be shared by precomputing these
  /// parts of the calculation into a buffer, which can then be passed several
  /// times to Evolve.
  /// Note that PrepareEvolve should be called on the evolution operator, and
  /// then Evolve should be called on the state(s) being evolved.
  ///\param buffer The buffer where the intermediate results will be stored.
  ///              Must be at least as large as the result of GetEvolveBufferSize
  ///\param t The time over which the evolution will be performed.
  void PrepareEvolve(double* buffer, double t) const{
    auto& suv1=*this;
    size_t offset=GetEvolveBufferSize()/2;
    double* CX=buffer;
    double* SX=buffer+offset;
    double term;
#include "SU_inc/PreEvolutionSelect.txt"
  }

  ///\brief Precompute operator dependent elements of an evolution
  ///
  /// Much of the Calculation needed by Evolve depends only on the evaolution
  /// operator and the time. If several SU_vectors will be evolved over the same
  /// time by the same operator, this work can be shared by precomputing these
  /// parts of the calculation into a buffer, which can then be passed several
  /// times to Evolve.
  /// Note that PrepareEvolve should be called on the evolution operator, and
  /// then Evolve should be called on the state(s) being evolved.
  ///\param buffer The buffer where the intermediate results will be stored.
  ///              Must be at least as large as the result of GetEvolveBufferSize
  ///\param t The time over which the evolution will be performed.
  ///\param scale Averaging frequency scale. If oscillation frequency is larger
  ///              than scale the oscillation will be averaged out.
  ///\param avr Vector of booleans signaling which scales were averaged out.
  void PrepareEvolve(double* buffer, double t, double scale, std::vector<bool>& avr) const{
    auto& suv1=*this;
    size_t offset=GetEvolveBufferSize()/2;
    double* CX=buffer;
    double* SX=buffer+offset;
    double term;
#include "SU_inc/PreEvolutionSelectAvg.txt"
  }

  void PrepareEvolve(double* buffer, double t_start, double t_end) const{
    auto& suv1=*this;
    size_t offset=GetEvolveBufferSize()/2;
    double* CX=buffer;
    double* SX=buffer+offset;
    double alpha;
    double range = t_end - t_start;
#include "SU_inc/PreEvolutionSelectAvgRange.txt"
  }
  
  ///\brief Compute low-pass filter for pre-computed sine and cosine evaluations
  ///
  /// This function applies a low-pass filter to the pre-evolution buffer computed by
  /// PrepareEvolve to filter out high frequencies. This is distinct from the already
  /// existing averaging mechanism because this function filters the frequency, not the
  /// number of rotations (i.e. frequency times time). This allows this function to be
  /// used both during evaluation of probabilities from the evolved state as well as
  /// during the integration of the interaction state equation. A linear ramp can be
  /// used to soften the cut-off of the filter.
  ///\param buffer  The buffer with evaluated sine/cosine values from PrepareEvolve.
  ///\param cutoff  Cut-off frequency of the filter. Sine and cosine evaluations 
  ///               with input frequencies higher than this will be set to zero
  ///\param scale   Distance in frequency between cut-off and pass-through. A value of 
  ///               greater than zero allows for a gradual transition to cut-off, rather
  ///               than applying a hard step function.
  void LowPassFilter(double* buffer, double cutoff, double scale) const{
    auto& suv1=*this;
    size_t offset=GetEvolveBufferSize()/2;
    double* CX=buffer;
    double* SX=buffer+offset;
    double term;
    int i;
#include "SU_inc/LowPassFilterSelect.txt"    
  }

  ///\brief Compute averaging filter for pre-computed sine and cosine evaluations.
  ///
  /// This function applies an averaging filter to the pre-evolution buffer computed by
  /// PrepareEvolve to filter out high frequencies in the same manner as the averaging
  /// mechanism, meaning that it applies a cut-off to the number of rotations rather
  /// than frequency. The difference to the already existing mechanism is that this
  /// filter can be applied as a post-processing step (like the low-pass filter), which
  /// means that it can also be applied after averaging over a distance. This filter is
  /// NOT appropriate for the RHS of the state density evolution! Use the low-pass
  /// filter instead for that purpose.
  ///\param buffer  The buffer with evaluated sine/cosine values from PrepareEvolve.
  ///\param cutoff  Cut-off frequency of the filter. Sine and cosine evaluations 
  ///               with input (frequencies * time) higher than this will be set to zero
  ///\param scale   Distance in (frequency * time) between cut-off and pass-through. A
  ///               value of greater than zero allows for a gradual transition to
  ///               cut-off, rather than applying a hard step function.
  void AvgRampFilter(double* buffer, double cutoff, double scale) const{
    auto& suv1=*this;
    size_t offset=GetEvolveBufferSize()/2;
    double* CX=buffer;
    double* SX=buffer+offset;
    double term;
    int i;
#include "SU_inc/AvgRampFilterSelect.txt"
  }
  
  
  ///\brief Compute the time evolution of the SU_vector
  ///
  ///\param buffer A buffer which has been filled by a previous call to
  ///              PrepareEvolve on the evaolution operator.
  ///\returns An object convertible to an SU_vector
  detail::FastEvolutionProxy Evolve(const double* buffer) const{
    return(detail::FastEvolutionProxy{*this,buffer});
  }


  //**********
  //operators
  //**********
  ///\brief Equality comparison
  bool operator==(const SU_vector&) const;

  ///\brief Scalar product of two SU_vectors.
  ///
  /// Equivalent to the trace of the matrix multiplication.
  SQUIDS_ALWAYS_INLINE double operator*(const SU_vector& other) const{
    if(size!=other.size)
      throw std::runtime_error("Non-matching dimensions in SU_vector inner product");
    return SUTrace<>(*this,other);
  }

  ///\brief Multiplication by a scalar
  ///\returns An object convertible to an SU_vector
  detail::MultiplicationProxy operator*(const double x) const &{
    return(detail::MultiplicationProxy{*this,x});
  }

  ///\brief Multiplication by a scalar
  ///\returns An object convertible to an SU_vector
  detail::MultiplicationProxy operator*(const double x) &&{
    return(detail::MultiplicationProxy{*this,x,detail::Arg1Movable});
  }

  ///\brief Assignment
  ///
  /// Assignment will fail if this vector uses external storage and the
  /// dimensions of the two vectors differ.
  SU_vector & operator=(const SU_vector&);

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
  detail::AdditionProxy operator+(const SU_vector& other) const &{
    if(size!=other.size)
      throw std::runtime_error("Non-matching dimensions in SU_vector addition");
    return(detail::AdditionProxy{*this,other});
  }

  ///\brief Addition
  ///\returns An object convertible to an SU_vector
  detail::AdditionProxy operator+(SU_vector&& other) const &{
    if(size!=other.size)
      throw std::runtime_error("Non-matching dimensions in SU_vector addition");
    //exploit commutativity and put the movable object first
    return(detail::AdditionProxy{other,*this,detail::Arg1Movable});
  }

  ///\brief Addition
  ///\returns An object convertible to an SU_vector
  detail::AdditionProxy operator+(const SU_vector& other) &&{
    if(size!=other.size)
      throw std::runtime_error("Non-matching dimensions in SU_vector addition");
    return(detail::AdditionProxy{*this,other,detail::Arg1Movable});
  }

  ///\brief Addition
  ///\returns An object convertible to an SU_vector
  detail::AdditionProxy operator+(SU_vector&& other) &&{
    if(size!=other.size)
      throw std::runtime_error("Non-matching dimensions in SU_vector addition");
    return(detail::AdditionProxy{*this,other,detail::Arg1Movable|detail::Arg2Movable});
  }

  ///\brief Subtraction
  ///\returns An object convertible to an SU_vector
  detail::SubtractionProxy operator-(const SU_vector& other) const &{
    if(size!=other.size)
      throw std::runtime_error("Non-matching dimensions in SU_vector subtraction");
    return(detail::SubtractionProxy{*this,other});
  }

  ///\brief Subtraction
  ///\returns An object convertible to an SU_vector
  detail::SubtractionProxy operator -(const SU_vector& other) &&{
    if(size!=other.size)
      throw std::runtime_error("Non-matching dimensions in SU_vector subtraction");
    return(detail::SubtractionProxy{*this,other,detail::Arg1Movable});
  }

  ///\brief Negation
  ///\returns An object convertible to an SU_vector
  detail::NegationProxy operator-() const &{
    return(detail::NegationProxy{*this});
  }

  ///\brief Negation
  ///\returns An object convertible to an SU_vector
  detail::NegationProxy operator-() &&{
    return(detail::NegationProxy{*this,detail::Arg1Movable});
  }

  ///\brief Optimized assignment from the result of an arithmetic expression
  template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
  SQUIDS_ALWAYS_INLINE SU_vector& operator=(const ProxyType& proxy){
    return(assignProxy<detail::AssignWrapper>(proxy));
  }

  ///\brief Optimized incrementing assignment from the result of an
  /// arithmetic expression
  template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
  SQUIDS_ALWAYS_INLINE SU_vector& operator+=(const ProxyType& proxy){
    return(assignProxy<detail::IncrementWrapper>(proxy));
  }

  ///\brief Optimized decrementing assignment from the result of an
  /// arithmetic expression
  template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
  SQUIDS_ALWAYS_INLINE SU_vector& operator-=(const ProxyType& proxy){
    return(assignProxy<detail::DecrementWrapper>(proxy));
  }

  ///\brief Set the external storage used by this SU_vector
  void SetBackingStore(double* storage){
    if(isinit)
      deallocate_mem();
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

  ///\brief Construct a projection operator.
  /// The resulting vector will have optimal alignement.
  ///
  ///\param d The dimension of the operator
  ///\param i The dimension selected by the projection
  static SU_vector Projector(unsigned int d, unsigned int i);

  ///\brief Construct an identity operator.
  /// The resulting vector will have optimal alignement.
  ///
  ///\param d The dimension of the operator
  static SU_vector Identity(unsigned int d);

  ///\brief Constructs a projector to the upper subspace of dimension i.
  /// The resulting vector will have optimal alignement.
  ///
  ///\param d The dimension of the operator
  ///\param i The subspace dimension
  ///
  /// \f$ NegProj = diag(1,...,1,0,...,0) \f$ where the last one is at the \f$i-1\f$ entry.
  static SU_vector PosProjector(unsigned int d, unsigned int i);

  ///\brief Constructs a projector to the lower subspace of dimension i.
  /// The resulting vector will have optimal alignement.
  ///
  ///\param d The dimension of the operator
  ///\param i The subspace dimension
  ///
  /// \f$ NegProj = diag(0,...,0,1,...,1) \f$ where the first one is at the \f$d-i\f$ entry.
  static SU_vector NegProjector(unsigned int d, unsigned int i);

  ///\brief Creates a SU_vector corresponding to the \f$i\f$SU_N basis generator.
  /// The resulting vector will have optimal alignement.
  //
  // Tts represented by \f$ v = (0,...,1,...,0)\f$ where 1 is in the \f$i\f$ component.
  static SU_vector Generator(unsigned int d, unsigned int i);

  template<typename Op>
  friend struct detail::EvaluationProxy;
  friend struct detail::EvolutionProxy;
  friend struct detail::FastEvolutionProxy;
  friend struct detail::AdditionProxy;
  friend struct detail::SubtractionProxy;
  friend struct detail::NegationProxy;
  friend struct detail::MultiplicationProxy;
  friend struct detail::iCommutatorProxy;
  friend struct detail::ACommutatorProxy;
  template<typename Op>
  friend struct detail::BinaryElementwiseOpProxy;

  friend struct detail::SU_vector_operator_access;

  //overloaded output operator
  friend std::ostream& operator<<(std::ostream&, const SU_vector&);
  
#if SQUIDS_USE_STORAGE_CACHE
  ///\brief Testing function which discards all cached memory blocks
  ///In normal use there is no reason to call this function, as it prevents the
  ///quick reuse of memory which has previously been allocated
  static void clear_mem_cache(){
    for(unsigned int dim=1; dim<=SQUIDS_MAX_HILBERT_DIM; dim++){
      mem_cache_entry cache_result;
      while(true){
        cache_result=storage_cache[dim].get();
        if(cache_result.storage)
          delete[] (cache_result.storage-cache_result.offset);
        else
          break;
      }
    }
  }
#endif
};

namespace detail{
//a helper class which gives externally defined operators access to the internals of SU_vector
struct SU_vector_operator_access{
  struct view{
    unsigned int& dim;
    unsigned int& size;
    double*& components;
    view(SU_vector& v):
    dim(v.dim),size(v.size),components(v.components){}
  };
  struct const_view{
    const unsigned int& dim;
    const unsigned int& size;
    const double* components;
    const_view(const SU_vector& v):
    dim(v.dim),size(v.size),components(v.components){}
  };
  static view make_view(SU_vector& v){ return(view(v)); }
  static const_view make_view(const SU_vector& v){ return(const_view(v)); }
};
} //namespace detail

///\brief Commutator of two SU_vectors
template<typename=void>
detail::iCommutatorProxy iCommutator(const SU_vector& suv1,const SU_vector& suv2){
  if(suv1.Dim()!=suv2.Dim())
    throw std::runtime_error("Commutator error, non-matching dimensions ");
  return(detail::iCommutatorProxy{suv1,suv2});
}

///\brief Anticommutator of two SU_vectors
template<typename=void>
detail::ACommutatorProxy ACommutator(const SU_vector& suv1,const SU_vector& suv2){
  if(suv1.Dim()!=suv2.Dim())
    throw std::runtime_error("Anti Commutator error: non-matching dimensions ");
  return(detail::ACommutatorProxy{suv1,suv2});
}

///\brief Multiplication of an SU_vector by a scalar from the left.
template<typename=void>
detail::MultiplicationProxy operator*(double x, const SU_vector& v){
  return(detail::MultiplicationProxy{v,x});
}

///\brief Multiplication of an SU_vector by a scalar from the left.
template<typename=void>
detail::MultiplicationProxy operator*(double x, SU_vector&& v){
  return(detail::MultiplicationProxy{v,x,detail::Arg1Movable});
}
  
///\brief Apply an arbitrary binary operation element-wise to two SU_vectors.
///\param op the operation to apply
///\param suv1 the 'left' or first operand
///\param suv2 the 'right' or second operand
///\return (an object convertible to) an SU_vector whose components are equal to
///        the result of applying op to the corresponding elements of suv1 and
///        suv2.
template<typename Op>
detail::BinaryElementwiseOpProxy<Op> ElementwiseOperation(Op op, const SU_vector& suv1, const SU_vector& suv2){
  if(suv1.Dim()!=suv2.Dim())
    throw std::runtime_error("Non-matching dimensions in SU_vector element-wise operation");
  return(detail::BinaryElementwiseOpProxy<Op>{op,suv1,suv2});
}
///\brief Apply an arbitrary binary operation element-wise to two SU_vectors.
///This overload can reuse the storage of the first operand.
///\param op the operation to apply
///\param suv1 the 'left' or first operand
///\param suv2 the 'right' or second operand
///\return (an object convertible to) an SU_vector whose components are equal to
///        the result of applying op to the corresponding elements of suv1 and
///        suv2.
template<typename Op>
detail::BinaryElementwiseOpProxy<Op> ElementwiseOperation(Op op, SU_vector&& suv1, const SU_vector& suv2){
  if(suv1.Dim()!=suv2.Dim())
    throw std::runtime_error("Non-matching dimensions in SU_vector element-wise operation");
  return(detail::BinaryElementwiseOpProxy<Op>{op,suv1,suv2,detail::Arg1Movable});
}
///\brief Apply an arbitrary binary operation element-wise to two SU_vectors.
///This overload can reuse the storage of the second operand.
///\param op the operation to apply
///\param suv1 the 'left' or first operand
///\param suv2 the 'right' or second operand
///\return (an object convertible to) an SU_vector whose components are equal to
///        the result of applying op to the corresponding elements of suv1 and
///        suv2.
template<typename Op>
detail::BinaryElementwiseOpProxy<Op> ElementwiseOperation(Op op, const SU_vector& suv1, SU_vector&& suv2){
  if(suv1.Dim()!=suv2.Dim())
    throw std::runtime_error("Non-matching dimensions in SU_vector element-wise operation");
  return(detail::BinaryElementwiseOpProxy<Op>{op,suv1,suv2,detail::Arg2Movable});
}
///\brief Apply an arbitrary binary operation element-wise to two SU_vectors.
///This overload can reuse the storage of the first or second operand.
///\param op the operation to apply
///\param suv1 the 'left' or first operand
///\param suv2 the 'right' or second operand
///\return (an object convertible to) an SU_vector whose components are equal to
///        the result of applying op to the corresponding elements of suv1 and
///        suv2.
template<typename Op>
detail::BinaryElementwiseOpProxy<Op> ElementwiseOperation(Op op, SU_vector&& suv1, SU_vector&& suv2){
  if(suv1.Dim()!=suv2.Dim())
    throw std::runtime_error("Non-matching dimensions in SU_vector element-wise operation");
  return(detail::BinaryElementwiseOpProxy<Op>{op,suv1,suv2,detail::Arg1Movable|detail::Arg2Movable});
}

//TODO: It should be possible to do suitable forwarding without four distinct copies of the function!
///\brief Element-wise product of two SU_vectors
template<typename=void>
detail::BinaryElementwiseOpProxy<std::multiplies<double>> ElementwiseProduct(const SU_vector& suv1, const SU_vector& suv2){
  return(ElementwiseOperation(std::multiplies<double>(),suv1,suv2));
}
///\brief Element-wise product of two SU_vectors
template<typename=void>
detail::BinaryElementwiseOpProxy<std::multiplies<double>> ElementwiseProduct(SU_vector&& suv1, const SU_vector& suv2){
  return(ElementwiseOperation(std::multiplies<double>(),std::move(suv1),suv2));
}
///\brief Element-wise product of two SU_vectors
template<typename=void>
detail::BinaryElementwiseOpProxy<std::multiplies<double>> ElementwiseProduct(const SU_vector& suv1, SU_vector&& suv2){
  return(ElementwiseOperation(std::multiplies<double>(),suv1,std::move(suv2)));
}
///\brief Element-wise product of two SU_vectors
template<typename=void>
detail::BinaryElementwiseOpProxy<std::multiplies<double>> ElementwiseProduct(SU_vector&& suv1, SU_vector&& suv2){
  return(ElementwiseOperation(std::multiplies<double>(),std::move(suv1),std::move(suv2)));
}

///\brief Gets the exponential of a GSL complex matrix
void gsl_complex_matrix_exponential(gsl_matrix_complex *eA, const gsl_matrix_complex *A, unsigned int dimx);
} //namespace squids

#include "detail/ProxyImpl.h"

#endif
