#ifndef SQUIDS_DETAIL_PROXYIMPL_H
#define SQUIDS_DETAIL_PROXYIMPL_H

namespace detail{
  constexpr static int Arg1Movable=1;
  constexpr static int Arg2Movable=2;
  
  template<typename Op>
  struct operation_traits{
    constexpr static bool elementwise=false;
  };
  
  template<typename Op>
  struct EvaluationProxy{
    const SU_vector& suv1;
    const SU_vector& suv2;
    int flags;
    
    template<typename VW>
    void compute(VW target) const{
      static_cast<const Op*>(this)->compute(target);
    }
    operator SU_vector() const{
      SU_vector result(suv1.dim);
      compute(detail::vector_wrapper<detail::AssignWrapper>{result.dim,result.components});
      return(result);
    }
    
    SU_vector operator*(double a) const;
    SU_vector operator+(const SU_vector& other) const;
    SU_vector operator-(const SU_vector& other) const;
    SU_vector SUEvolve(const SU_vector& other ,double t) const;
    
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY>
    double operator*(const ProxyType& other) const{
      return(SU_vector(*this)*SU_vector(other));
    }
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY>
    SU_vector operator+(const ProxyType& other) const;
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY>
    SU_vector operator-(const ProxyType& other) const;
    template<typename ProxyType, REQUIRE_EVALUATION_PROXY>
    SU_vector SUEvolve(const ProxyType& other ,double t) const;
    
    bool mayStealArg1() const{
      return(flags&detail::Arg1Movable && operation_traits<Op>::elementwise);
    }
  };
  
  struct EvolutionProxy : public EvaluationProxy<EvolutionProxy>{
    double t;
    
    EvolutionProxy(const SU_vector& suv1,const SU_vector& suv2,double t):
    EvaluationProxy<EvolutionProxy>{suv1,suv2,0},t(t){}
    
    template<typename VW>
    void compute(VW target) const{
      auto& suv_new=target;
#include "../SU_inc/EvolutionSelect.txt"
    }
  };
  
  struct AdditionProxy : public EvaluationProxy<AdditionProxy>{
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
  
  struct SubtractionProxy : public EvaluationProxy<SubtractionProxy>{
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
  
  struct MultiplicationProxy : public EvaluationProxy<MultiplicationProxy>{
    double a;
    
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
  
  struct iCommutatorProxy : public EvaluationProxy<iCommutatorProxy>{
    iCommutatorProxy(const SU_vector& suv1,const SU_vector& suv2):
    EvaluationProxy<iCommutatorProxy>{suv1,suv2,0}{}
    
    template<typename VW>
    void compute(VW suv_new) const{
      suv_new.components[0]+=0;
#include "../SU_inc/iCommutatorSelect.txt"
    }
  };
  
  struct ACommutatorProxy : public EvaluationProxy<ACommutatorProxy>{
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

#undef REQUIRE_EVALUATION_PROXY

#endif //SQUIDS_DETAIL_PROXYIMPL_H