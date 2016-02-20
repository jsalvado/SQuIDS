#ifndef SQUIDS_DETAIL_PROXYFWD_H
#define SQUIDS_DETAIL_PROXYFWD_H

#define REQUIRE_EVALUATION_PROXY_CORE typename std::enable_if<std::is_base_of<detail::EvaluationProxy<ProxyType>,ProxyType>::value>
#define REQUIRE_EVALUATION_PROXY_TPARAM typename= REQUIRE_EVALUATION_PROXY_CORE
#define REQUIRE_EVALUATION_PROXY_FPARAM REQUIRE_EVALUATION_PROXY_CORE ::type* =nullptr

//Try to figue out how to spell hints to the optimizer
#ifdef __clang__
  #if __has_builtin(__builtin_assume)
    #define SQUIDS_COMPILER_ASSUME(axiom) __builtin_assume(axiom)
  #else
    #if __has_builtin(__builtin_unreachable)
      #define SQUIDS_COMPILER_ASSUME(axiom) if(not (axiom)) __builtin_unreachable()
    #else
      #define SQUIDS_COMPILER_ASSUME(axiom) do{}while(0) //not available
    #endif
  #endif
  #if __has_builtin(__builtin_assume_aligned)
    #define SQUIDS_POINTER_IS_ALIGNED(ptr,alignment) \
      ptr=(decltype(ptr))__builtin_assume_aligned(ptr,alignment)
  #else
    #define SQUIDS_POINTER_IS_ALIGNED(ptr,alignment) do{}while(0) //not available
  #endif
#endif
#if defined(__GNUC__) && !defined(__clang__)
  //All GCC versions otherwise able to compile the code should support this
  #define SQUIDS_COMPILER_ASSUME(axiom) if(not (axiom)) __builtin_unreachable()
  #define SQUIDS_POINTER_IS_ALIGNED(ptr,alignment) \
    ptr=(decltype(ptr))__builtin_assume_aligned(ptr,alignment)
#endif
#ifdef __INTEL_COMPILER
  //Documentation seems to be lacking on when this became supported or how to tell
  #define SQUIDS_COMPILER_ASSUME(axiom) __assume(axiom)
  #define SQUIDS_POINTER_IS_ALIGNED(ptr,alignment) __assume_aligned(ptr,alignment)
#endif
#ifdef _MSC_VER
  //If any version is able to compile the library, it should already have this
  #define SQUIDS_COMPILER_ASSUME(axiom) __assume(axiom)
  #define SQUIDS_POINTER_IS_ALIGNED(ptr,alignment) \
    __assume((uintptr_t(ptr) & ((alignment) - 1)) == 0)
#endif
//final fallback
#ifndef SQUIDS_COMPILER_ASSUME
  #define SQUIDS_COMPILER_ASSUME(axiom) do{}while(0) //not available
  #define SQUIDS_POINTER_IS_ALIGNED(ptr,alignment) do{}while(0) //not available
#endif

#ifndef SQUIDS_USE_VECTOR_EXTENSIONS
  #if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
    #define SQUIDS_USE_VECTOR_EXTENSIONS 1
  #endif
#endif

namespace squids{
  
class SU_vector;

///This namespace contains implementation details
///which most users should not need to use directly
namespace detail{
    
#if SQUIDS_USE_VECTOR_EXTENSIONS
  //We will try to use SIMD vectors of 4 doubles.
  //This is both the best size for us (8 is not always applicable), and if size
  //4 vectors are not available the compiler should legalize them to size 2 or
  //scalars.
  #if defined(__GNUC__) || defined(__clang__)
    typedef double double_vector4 __attribute__((__vector_size__(32)));
  #endif
  #if defined(__INTEL_COMPILER)
    typedef double double_vector4 __attribute__((vector vectorlength(4)));
  #endif
#endif //SQUIDS_USE_VECTOR_EXTENSIONS
  
  ///All of the SU_vector operations in ../SU_inc/ are generated in terms
  ///of incmrement operations (+=), but we would like to be able to fuse
  ///them with plain assignments, increments, and decrements. Rather that
  ///creating three versions of those operations with tiny syntax differences
  ///this set of classes allows them to be used for all three operations
  ///with the same syntax, by means of the slightly underhanded trick of
  ///redefining the += operator
  
  struct AssignWrapper{
    constexpr static bool allowTargetResize=true;
    double* v;
    double& operator+=(double nv){
      *v=nv;
      return(*v);
    }
#if SQUIDS_USE_VECTOR_EXTENSIONS
    void apply4(double_vector4 rhs){
      *(double_vector4*)v=rhs;
    }
#endif //SQUIDS_USE_VECTOR_EXTENSIONS
    operator double() const{ return(*v); }
    template<typename T>
    static T& apply(T& target,const T& source){
      return(target=source);
    }
  };
  
  struct IncrementWrapper{
    constexpr static bool allowTargetResize=false;
    double* v;
    double& operator+=(double nv){
      *v+=nv;
      return(*v);
    }
#if SQUIDS_USE_VECTOR_EXTENSIONS
    void apply4(double_vector4 rhs){
      *(double_vector4*)v=*(double_vector4*)v + rhs;
    }
#endif //SQUIDS_USE_VECTOR_EXTENSIONS
    operator double() const{ return(*v); }
    template<typename T>
    static T& apply(T& target,const T& source){
      return(target+=source);
    }
  };
  
  struct DecrementWrapper{
    constexpr static bool allowTargetResize=false;
    double* v;
    double& operator+=(double nv){
      *v-=nv;
      return(*v);
    }
#if SQUIDS_USE_VECTOR_EXTENSIONS
    void apply4(double_vector4 rhs){
      *(double_vector4*)v=*(double_vector4*)v - rhs;
    }
#endif //SQUIDS_USE_VECTOR_EXTENSIONS
    operator double() const{ return(*v); }
    template<typename T>
    static T& apply(T& target,const T& source){
      return(target-=source);
    }
  };
  
  ///This class is used to adapt a plain array of values by wrapping
  ///each value with the given Wrapper type on demand.
  template<typename Wrapper>
  struct vector_wrapper{
    const unsigned int& dim;
    struct component_wrapper{
      double* components;
      
      component_wrapper(double* components):components(components){}
      Wrapper operator[](unsigned int i){
        return(Wrapper{components+i});
      }
    } components;
    vector_wrapper(const unsigned int& dim, double* components):
    dim(dim),components{components}{}
#if SQUIDS_USE_VECTOR_EXTENSIONS
    void apply4(double_vector4 rhs){ Wrapper{components.components}.apply4(rhs); };
#endif //SQUIDS_USE_VECTOR_EXTENSIONS
  };
  
  ///Constant used to indicate that the first argument of a proxy was an
  ///r-value reference and may be moved from
  constexpr static int Arg1Movable=1;
  ///Constant used to indicate that the second argument of a proxy was an
  ///r-value reference and may be moved from
  constexpr static int Arg2Movable=2;
  
  ///Used to look up properties of operations on SU_vectors
  template<typename Op>
  struct operation_traits{
    ///Indicates whether the calculation of each element of the result depends
    ///only on the corresponding elements of the operands, meaning that the
    ///calculation of the result elements may be reordered, or the inputs may
    ///be overwritten
    constexpr static bool elementwise=false;
    ///The number of SU_vector operands.
    ///If 1, suv2 is expected to be ignored
    constexpr static unsigned int vector_arity=2;
    ///The operands of this operation are guaranteed not to alias the location
    ///where the results will be stored
    constexpr static bool no_alias_target=false;
    ///The operands of this operation are guaranteed not to have the same size
    ///as the location where the results will be stored
    constexpr static bool equal_target_size=false;
    ///Both the operands of this operation the location where the results will
    ///be stored are guaranteed to have optimally aligned storage.
    constexpr static bool aligned_storgae=false;
  };
  
  ///The base class for objects representing arithmetic operations on SU_vectors.
  ///Rather than eagerly performing arithemtic, the calculation's type and operands
  ///are encoded in a proxy object which knows how to compute the operation once
  ///more is known about the context (such as where the result will be stored).
  ///This enables a number of operations which reduce copying of data.
  template<typename Op>
  struct EvaluationProxy{
    const SU_vector& suv1; ///the first operand
    const SU_vector& suv2; ///the second operand
    int flags;
    
    ///perform the stored operation, storing the results in target
    template<typename VW, bool Aligned=false>
    void compute(VW target) const{
      static_cast<const Op*>(this)->compute(target);
    }
    ///compute the stored operation into a temporary SU_vector
    operator SU_vector() const &;
    ///compute the stored operation, stealing memory if possible
    operator SU_vector() &&;
    
    ///multiplication by scalars
    SU_vector operator*(double a) const;
    ///addition with SU_vectors
    SU_vector operator+(const SU_vector& other) const;
    ///subtraction with SU_vectors
    SU_vector operator-(const SU_vector& other) const;
    ///time evolution according to an SU_vector operator
    SU_vector Evolve(const SU_vector& other, double t) const;
    
    ///negation
    SU_vector operator-() const &;
    SU_vector operator-() &&;
    
    ///scalar product between proxies
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
    double operator*(const ProxyType& other) const;
    ///addition of proxies
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
    SU_vector operator+(const ProxyType& other) const;
    ///subtraction of proxies
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
    SU_vector operator-(const ProxyType& other) const;
    ///time evolution of a proxy according to an operator given by another proxy
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
    SU_vector Evolve(const ProxyType& other ,double t) const;
    
    ///whether the result of the operation can be written directly
    ///into the storage of the first operand
    bool mayStealArg1() const{
      return(flags&detail::Arg1Movable && operation_traits<Op>::elementwise);
    }
  };
  
  //fwd decls for operation proxies
  
  ///The result of a time evolution operation
  struct EvolutionProxy : public EvaluationProxy<EvolutionProxy>{
    double t; ///the time over which the evolution is performed
    
    ///evolve suv1 according to the operator suv2 over time period t
    EvolutionProxy(const SU_vector& suv1,const SU_vector& suv2,double t):
    EvaluationProxy<EvolutionProxy>{suv1,suv2,0},t(t){}
    
    template<typename VW, bool Aligned=false>
    void compute(VW target) const;
  };
  
  ///The result of a time evolution operation
  struct FastEvolutionProxy : public EvaluationProxy<FastEvolutionProxy>{
    const double* coefficients;
    
    ///evolve suv1 according to the operator suv2 over time period t
    FastEvolutionProxy(const SU_vector& suv1,const double* c):
    EvaluationProxy<FastEvolutionProxy>{suv1,suv1,0},coefficients(c){}
    
    template<typename VW, bool Aligned=false>
    void compute(VW target) const;
  };
  
  template<>
  struct operation_traits<FastEvolutionProxy>{
    constexpr static bool elementwise=false;
    constexpr static unsigned int vector_arity=1;
    constexpr static bool no_alias_target=false;
    constexpr static bool equal_target_size=false;
    constexpr static bool aligned_storgae=false;
  };
  
  ///The result of adding two SU_vectors
  struct AdditionProxy : public EvaluationProxy<AdditionProxy>{
    ///The sum of suv1 and suv2
    AdditionProxy(const SU_vector& suv1,const SU_vector& suv2,int flags=0):
    EvaluationProxy<AdditionProxy>{suv1,suv2,flags}{}
    
    template<typename VW, bool Aligned=false>
    void compute(VW target) const;
  };
  
  template<>
  struct operation_traits<AdditionProxy>{
    constexpr static bool elementwise=true;
    constexpr static unsigned int vector_arity=2;
    constexpr static bool no_alias_target=false;
    constexpr static bool equal_target_size=false;
    constexpr static bool aligned_storgae=false;
  };
  
  ///The result of subtracting two SU_vectors
  struct SubtractionProxy : public EvaluationProxy<SubtractionProxy>{
    ///The difference of suv1 and suv2
    SubtractionProxy(const SU_vector& suv1,const SU_vector& suv2,int flags=0):
    EvaluationProxy<SubtractionProxy>{suv1,suv2,flags}{}
    
    template<typename VW, bool Aligned=false>
    void compute(VW target) const;
  };
  
  template<>
  struct operation_traits<SubtractionProxy>{
    constexpr static bool elementwise=true;
    constexpr static unsigned int vector_arity=2;
    constexpr static bool no_alias_target=false;
    constexpr static bool equal_target_size=false;
    constexpr static bool aligned_storgae=false;
  };
  
  ///The result of negating an SU_vector
  struct NegationProxy : public EvaluationProxy<NegationProxy>{
    ///The additive inverse of suv
    NegationProxy(const SU_vector& suv, int flags=0):
    EvaluationProxy<NegationProxy>{suv,suv,flags}{}
    
    template<typename VW, bool Aligned=false>
    void compute(VW target) const;
  };
  
  template<>
  struct operation_traits<NegationProxy>{
    constexpr static bool elementwise=true;
    constexpr static unsigned int vector_arity=1;
    constexpr static bool no_alias_target=false;
    constexpr static bool equal_target_size=false;
    constexpr static bool aligned_storgae=false;
  };
  
  ///The result of multiplying an SU_vector by a scalar
  struct MultiplicationProxy : public EvaluationProxy<MultiplicationProxy>{
    double a; ///scalar mutiplcation factor
    
    ///The product of suv1 and a
    MultiplicationProxy(const SU_vector& suv1,double a,int flags=0):
    EvaluationProxy<MultiplicationProxy>{suv1,suv1,flags},a(a){}
    
    template<typename VW, bool Aligned=false>
    void compute(VW target) const;
  };
  
  template<>
  struct operation_traits<MultiplicationProxy>{
    constexpr static bool elementwise=true;
    constexpr static unsigned int vector_arity=1;
    constexpr static bool no_alias_target=false;
    constexpr static bool equal_target_size=false;
    constexpr static bool aligned_storgae=false;
  };
  
  ///The result of i times the commutator of two SU_vectors
  struct iCommutatorProxy : public EvaluationProxy<iCommutatorProxy>{
    ///The commutator of suv1 and suv2
    iCommutatorProxy(const SU_vector& suv1,const SU_vector& suv2):
    EvaluationProxy<iCommutatorProxy>{suv1,suv2,0}{}
    
    template<typename VW, bool Aligned=false>
    void compute(VW suv_new) const;
  };
  
  ///The result of the anticommutator of two SU_vectors
  struct ACommutatorProxy : public EvaluationProxy<ACommutatorProxy>{
    ///The anticommutator of suv1 and suv2
    ACommutatorProxy(const SU_vector& suv1,const SU_vector& suv2):
    EvaluationProxy<ACommutatorProxy>{suv1,suv2,0}{}
    
    template<typename VW, bool Aligned=false>
    void compute(VW target) const;
  };
  
  /// Used to indicate that no SU_vector operand on the RHS of an assignment is
  /// the same as or shares storage with the LHS.
  constexpr static unsigned int NoAlias=1;
  /// Used to indicate that the SU_vector on the LHS of an assignment is known to
  /// already have the correct dimensions to match whatever is on the RHS.
  constexpr static unsigned int EqualSizes=2;
  /// Used to indicate that all SU_vector operands participating in an operation
  /// have ideally aligned storage. For this purpose storage is considered
  /// ideally aligned if either: the SU_vector has dimension 1, the SU_vector
  /// has an even number of components, and its first component is 32 byte
  /// aligned, or the SU_vector has an odd number of compoennts, larger than
  /// one, and the second component (index 1) is 32 byte aligned.
  constexpr static unsigned int AlignedStorage=4;
  
  /// Used to apply guarantees about operands or the destination storage of an
  /// EvaluationProxy.
  template<unsigned int Flags, typename WrappedType>
  struct GuaranteeWrapper : public WrappedType{
    constexpr static unsigned int Guarantees=Flags;
    
    GuaranteeWrapper(WrappedType op):WrappedType(op){}
    
    template<typename VW>
    void compute(VW target) const{
      WrappedType::template compute<VW,bool(Guarantees&AlignedStorage)>(target);
    }
  };
  
  ///\brief Convenience function for applying guarantees about operands or the
  /// destination storage to an EvaluationProxy.
  ///
  /// Typical use:
  ///
  ///     target = guarantee</*bitwise disjunction of guarantee flags*/>
  ///                       (Op(operand_1, operand_2));
  ///
  /// where target, operand_1, and operand_2 are SU_vectors, and Op is an
  /// operation (possibly an operator) returning and EvaluationProxy.
  ///
  /// Valid guarantee flags are NoAlias, EqualSizes, and AlignedStorage.
  /// Incorrect application of guarantees (guarantee which do not actually hold)
  /// must be considered to cause undefined behavior.
  template<unsigned int Flags, typename WrappedType>
  GuaranteeWrapper<Flags,WrappedType> guarantee(WrappedType w){
    return(GuaranteeWrapper<Flags,WrappedType>{w});
  }
  
  template<unsigned int Flags, typename WrappedType>
  struct operation_traits<GuaranteeWrapper<Flags,WrappedType>>{
    using base_traits=operation_traits<WrappedType>;
    constexpr static bool elementwise=base_traits::elementwise;
    constexpr static unsigned int vector_arity=base_traits::vector_arity;
    constexpr static bool no_alias_target=Flags&NoAlias;
    constexpr static bool equal_target_size=Flags&EqualSizes;
    constexpr static bool aligned_storgae=Flags&AlignedStorage;
  };
}
  
} //namespace SQuIDS

#endif //SQUIDS_DETAIL_PROXYFWD_H