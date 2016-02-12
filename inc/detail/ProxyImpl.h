#ifndef SQUIDS_DETAIL_PROXYIMPL_H
#define SQUIDS_DETAIL_PROXYIMPL_H

namespace squids{

namespace detail{
  
  
  template<typename Op>
  ///compute the stored operation into a temporary SU_vector
  EvaluationProxy<Op>::operator SU_vector() const &{
    SU_vector result(suv1.dim);
    compute(detail::vector_wrapper<detail::AssignWrapper>{result.dim,result.components});
    return(result);
  }
  template<typename Op>
  ///compute the stored operation, stealing memory if possible
  EvaluationProxy<Op>::operator SU_vector() &&{
    SU_vector result(static_cast<Op&&>(*this),static_cast<Op*>(nullptr));
    return(result);
  }
  
  template<typename Op>
  template<typename ProxyType, typename>
  double EvaluationProxy<Op>::operator*(const ProxyType& other) const{
    return(SU_vector(*this)*SU_vector(other));
  }
    
  template<typename VW>
  void EvolutionProxy::compute(VW target) const{
    auto& suv_new=target; //alias for the name expected by generated code
#include "../SU_inc/EvolutionSelect.txt"
  }
  
  template<typename VW>
  void AdditionProxy::compute(VW target) const{
    for(unsigned int i=0; i<target.dim*target.dim; i++)
      target.components[i] += suv1.components[i] + suv2.components[i];
  }
  
  template<>
  struct operation_traits<AdditionProxy>{
    constexpr static bool elementwise=true;
  };
  
  template<typename VW>
  void SubtractionProxy::compute(VW target) const{
    for(unsigned int i=0; i<target.dim*target.dim; i++)
      target.components[i] += suv1.components[i] - suv2.components[i];
  }
  
  template<>
  struct operation_traits<SubtractionProxy>{
    constexpr static bool elementwise=true;
  };
  
  template<typename VW>
  void NegationProxy::compute(VW target) const{
    for(unsigned int i=0; i<target.dim*target.dim; i++)
      target.components[i] += -suv1.components[i];
  }
  
  template<>
  struct operation_traits<NegationProxy>{
    constexpr static bool elementwise=true;
  };
  
  template<typename VW>
  void MultiplicationProxy::compute(VW target) const{
      double* suv1c=suv1.components;
      for(unsigned int i=0; i<target.dim*target.dim; i++)
        target.components[i] += a*suv1c[i];
  }
  
  template<>
  struct operation_traits<MultiplicationProxy>{
    constexpr static bool elementwise=true;
  };
  
  template<typename VW>
  void iCommutatorProxy::compute(VW suv_new) const{
    suv_new.components[0]+=0;
#include "../SU_inc/iCommutatorSelect.txt"
  }
  
  template<typename VW>
  void ACommutatorProxy::compute(VW suv_new) const{
#include "../SU_inc/AnticommutatorSelect.txt"
  }
  
  template<typename Op>
  SU_vector EvaluationProxy<Op>::operator*(double a) const{
    return(static_cast<SU_vector>(*this) * a);
  }
  template<typename Op>
  SU_vector EvaluationProxy<Op>::operator+(const SU_vector& other) const{
    return(static_cast<SU_vector>(*this) + other);
  }
  template<typename Op>
  SU_vector EvaluationProxy<Op>::operator-(const SU_vector& other) const{
    return(static_cast<SU_vector>(*this) - other);
  }
  template<typename Op>
  SU_vector EvaluationProxy<Op>::Evolve(const SU_vector& other ,double t) const{
    return(static_cast<SU_vector>(*this).Evolve(other,t));
  }
  
  template<typename Op>
  SU_vector EvaluationProxy<Op>::operator-() const &{
    return(-static_cast<SU_vector>(*this));
  }
  template<typename Op>
  SU_vector EvaluationProxy<Op>::operator-() &&{
    return(-static_cast<SU_vector>(std::move(*this)));
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
  SU_vector EvaluationProxy<Op>::Evolve(const ProxyType& other ,double t) const{
    return(SU_vector(*this).Evolve(other,t));
  }
} //namespace detail

#undef REQUIRE_EVALUATION_PROXY_CORE
#undef REQUIRE_EVALUATION_PROXY_TPARAM
#undef REQUIRE_EVALUATION_PROXY_FPARAM

} //namespace squids
  
#endif //SQUIDS_DETAIL_PROXYIMPL_H