#ifndef SQUIDS_DETAIL_PROXYFWD_H
#define SQUIDS_DETAIL_PROXYFWD_H

#define REQUIRE_EVALUATION_PROXY_CORE typename std::enable_if<std::is_base_of<detail::EvaluationProxy<ProxyType>,ProxyType>::value>
#define REQUIRE_EVALUATION_PROXY_TPARAM typename= REQUIRE_EVALUATION_PROXY_CORE
#define REQUIRE_EVALUATION_PROXY_FPARAM REQUIRE_EVALUATION_PROXY_CORE ::type* =nullptr

namespace squids{
  
class SU_vector;

///This namespace contains implementation details
///which most users should not need to use directly
namespace detail{
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
    private:
      double* components;
    public:
      component_wrapper(double* components):components(components){}
      Wrapper operator[](unsigned int i){
        return(Wrapper{components+i});
      }
    } components;
    vector_wrapper(const unsigned int& dim, double* components):
    dim(dim),components{components}{}
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
    template<typename VW>
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
    
    template<typename VW>
    void compute(VW target) const;
  };
  
  ///The result of adding two SU_vectors
  struct AdditionProxy : public EvaluationProxy<AdditionProxy>{
    ///The sum of suv1 and suv2
    AdditionProxy(const SU_vector& suv1,const SU_vector& suv2,int flags=0):
    EvaluationProxy<AdditionProxy>{suv1,suv2,flags}{}
    
    template<typename VW>
    void compute(VW target) const;
  };
  
  ///The result of subtracting two SU_vectors
  struct SubtractionProxy : public EvaluationProxy<SubtractionProxy>{
    ///The difference of suv1 and suv2
    SubtractionProxy(const SU_vector& suv1,const SU_vector& suv2,int flags=0):
    EvaluationProxy<SubtractionProxy>{suv1,suv2,flags}{}
    
    template<typename VW>
    void compute(VW target) const;
  };
  
  ///The result of negating an SU_vector
  struct NegationProxy : public EvaluationProxy<NegationProxy>{
    ///The additive inverse of suv
    NegationProxy(const SU_vector& suv, int flags=0):
    EvaluationProxy<NegationProxy>{suv,suv,flags}{}
    
    template<typename VW>
    void compute(VW target) const;
  };
  
  ///The result of multiplying an SU_vector by a scalar
  struct MultiplicationProxy : public EvaluationProxy<MultiplicationProxy>{
    double a; ///scalar mutiplcation factor
    
    ///The product of suv1 and a
    MultiplicationProxy(const SU_vector& suv1,double a,int flags=0):
    EvaluationProxy<MultiplicationProxy>{suv1,suv1,flags},a(a){}
    
    template<typename VW>
    void compute(VW target) const;
  };
  
  ///The result of i times the commutator of two SU_vectors
  struct iCommutatorProxy : public EvaluationProxy<iCommutatorProxy>{
    ///The commutator of suv1 and suv2
    iCommutatorProxy(const SU_vector& suv1,const SU_vector& suv2):
    EvaluationProxy<iCommutatorProxy>{suv1,suv2,0}{}
    
    template<typename VW>
    void compute(VW suv_new) const;
  };
  
  ///The result of the anticommutator of two SU_vectors
  struct ACommutatorProxy : public EvaluationProxy<ACommutatorProxy>{
    ///The anticommutator of suv1 and suv2
    ACommutatorProxy(const SU_vector& suv1,const SU_vector& suv2):
    EvaluationProxy<ACommutatorProxy>{suv1,suv2,0}{}
    
    template<typename VW>
    void compute(VW target) const;
  };
}
  
} //namespace SQuIDS

#endif //SQUIDS_DETAIL_PROXYFWD_H