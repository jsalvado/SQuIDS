#ifndef SQUIDS_DETAIL_PROXYIMPL_H
#define SQUIDS_DETAIL_PROXYIMPL_H

namespace detail{
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
    operator SU_vector() const &{
      SU_vector result(suv1.dim);
      compute(detail::vector_wrapper<detail::AssignWrapper>{result.dim,result.components});
      return(result);
    }
    ///compute the stored operation, stealing memory if possible
    operator SU_vector() &&{
      SU_vector result(static_cast<Op&&>(*this),nullptr);
      return(result);
    }
    
    ///multiplication by scalars
    SU_vector operator*(double a) const;
    ///addition with SU_vectors
    SU_vector operator+(const SU_vector& other) const;
    ///subtraction with SU_vectors
    SU_vector operator-(const SU_vector& other) const;
    ///time evolution according to an SU_vector operator
    SU_vector SUEvolve(const SU_vector& other, double t) const;
    
    ///negation
    SU_vector operator-() const &;
    SU_vector operator-() &&;
    
    ///scalar product between proxies
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
    double operator*(const ProxyType& other) const{
      return(SU_vector(*this)*SU_vector(other));
    }
    ///addition of proxies
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
    SU_vector operator+(const ProxyType& other) const;
    ///subtraction of proxies
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
    SU_vector operator-(const ProxyType& other) const;
    ///time evolution of a proxy according to an operator given by another proxy
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY_TPARAM>
    SU_vector SUEvolve(const ProxyType& other ,double t) const;
    
    ///whether the result of the operation can be written directly
    ///into the storage of the first operand
    bool mayStealArg1() const{
      return(flags&detail::Arg1Movable && operation_traits<Op>::elementwise);
    }
  };
  
  ///The result of a time evolution operation
  struct EvolutionProxy : public EvaluationProxy<EvolutionProxy>{
    double t; ///the time over which the evolution is performed
    
    ///evolve suv1 according to the operator suv2 over time period t
    EvolutionProxy(const SU_vector& suv1,const SU_vector& suv2,double t):
    EvaluationProxy<EvolutionProxy>{suv1,suv2,0},t(t){}
    
    template<typename VW>
    void compute(VW target) const{
      auto& suv_new=target; //alias for the name expected by generated code
#include "../SU_inc/EvolutionSelect.txt"
    }
  };
  
  ///The result of adding two SU_vectors
  struct AdditionProxy : public EvaluationProxy<AdditionProxy>{
    ///The sum of suv1 and suv2
    AdditionProxy(const SU_vector& suv1,const SU_vector& suv2,int flags=0):
    EvaluationProxy<AdditionProxy>{suv1,suv2,flags}{}
    
    template<typename VW>
    void compute(VW target) const{
      for(int i=0; i<target.dim*target.dim; i++)
        target.components[i] += suv1.components[i] + suv2.components[i];
    }
  };
  
  template<>
  struct operation_traits<AdditionProxy>{
    constexpr static bool elementwise=true;
  };
  
  ///The result of subtracting two SU_vectors
  struct SubtractionProxy : public EvaluationProxy<SubtractionProxy>{
    ///The difference of suv1 and suv2
    SubtractionProxy(const SU_vector& suv1,const SU_vector& suv2,int flags=0):
    EvaluationProxy<SubtractionProxy>{suv1,suv2,flags}{}
    
    template<typename VW>
    void compute(VW target) const{
      for(int i=0; i<target.dim*target.dim; i++)
        target.components[i] += suv1.components[i] - suv2.components[i];
    }
  };
  
  template<>
  struct operation_traits<SubtractionProxy>{
    constexpr static bool elementwise=true;
  };
  
  ///The result of negating an SU_vector
  struct NegationProxy : public EvaluationProxy<NegationProxy>{
    ///The additive inverse of suv
    NegationProxy(const SU_vector& suv, int flags=0):
    EvaluationProxy<NegationProxy>{suv,suv,flags}{}
    
    template<typename VW>
    void compute(VW target) const{
      for(int i=0; i<target.dim*target.dim; i++)
        target.components[i] += -suv1.components[i];
    }
  };
  
  template<>
  struct operation_traits<NegationProxy>{
    constexpr static bool elementwise=true;
  };
  
  ///The result of multiplying an SU_vector by a scalar
  struct MultiplicationProxy : public EvaluationProxy<MultiplicationProxy>{
    double a; ///scalar mutiplcation factor
    
    ///The product of suv1 and a
    MultiplicationProxy(const SU_vector& suv1,double a,int flags=0):
    EvaluationProxy<MultiplicationProxy>{suv1,suv1,flags},a(a){}
    
    template<typename VW>
    void compute(VW target) const{
      for(int i=0; i<target.dim*target.dim; i++)
        target.components[i] += a*suv1.components[i];
    }
  };
  
  template<>
  struct operation_traits<MultiplicationProxy>{
    constexpr static bool elementwise=true;
  };
  
  ///The result of i times the commutator of two SU_vectors
  struct iCommutatorProxy : public EvaluationProxy<iCommutatorProxy>{
    ///The commutator of suv1 and suv2
    iCommutatorProxy(const SU_vector& suv1,const SU_vector& suv2):
    EvaluationProxy<iCommutatorProxy>{suv1,suv2,0}{}
    
    template<typename VW>
    void compute(VW suv_new) const{
      suv_new.components[0]+=0;
#include "../SU_inc/iCommutatorSelect.txt"
    }
  };
  
  ///The result of the anticommutator of two SU_vectors
  struct ACommutatorProxy : public EvaluationProxy<ACommutatorProxy>{
    ///The anticommutator of suv1 and suv2
    ACommutatorProxy(const SU_vector& suv1,const SU_vector& suv2):
    EvaluationProxy<ACommutatorProxy>{suv1,suv2,0}{}
    
    template<typename VW>
    void compute(VW suv_new) const{
#include "../SU_inc/AnticommutatorSelect.txt"
    }
  };
  
  template<typename Op>
  SU_vector EvaluationProxy<Op>::operator*(double a) const{
    return(((SU_vector)*this) * a);
  }
  template<typename Op>
  SU_vector EvaluationProxy<Op>::operator+(const SU_vector& other) const{
    return((SU_vector)*this + other);
  }
  template<typename Op>
  SU_vector EvaluationProxy<Op>::operator-(const SU_vector& other) const{
    return((SU_vector)*this - other);
  }
  template<typename Op>
  SU_vector EvaluationProxy<Op>::SUEvolve(const SU_vector& other ,double t) const{
    return(((SU_vector)*this).SUEvolve(other,t));
  }
  
  template<typename Op>
  SU_vector EvaluationProxy<Op>::operator-() const &{
    return(-(SU_vector)*this);
  }
  template<typename Op>
  SU_vector EvaluationProxy<Op>::operator-() &&{
    return(-(SU_vector)std::move(*this));
  }
  
  template<typename Op>
  template<typename ProxyType, typename>
  SU_vector EvaluationProxy<Op>::operator+(const ProxyType& other) const{
    return(SU_vector(*this)+SU_vector(other));
  }
  template<typename Op>
  template<typename ProxyType, typename>
  SU_vector EvaluationProxy<Op>::operator-(const ProxyType& other) const{
    return(SU_vector(*this)-SU_vector(other));
  }
  template<typename Op>
  template<typename ProxyType, typename>
  SU_vector EvaluationProxy<Op>::SUEvolve(const ProxyType& other ,double t) const{
    return(SU_vector(*this).SUEvolve(other,t));
  }
} //namespace detail

#undef REQUIRE_EVALUATION_PROXY_CORE
#undef REQUIRE_EVALUATION_PROXY_TPARAM
#undef REQUIRE_EVALUATION_PROXY_FPARAM

#endif //SQUIDS_DETAIL_PROXYIMPL_H