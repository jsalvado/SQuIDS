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
    
  template<typename VW, bool Aligned>
  void EvolutionProxy::compute(VW target) const{
    auto& suv_new=target; //alias for the name expected by generated code
#include "../SU_inc/EvolutionSelect.txt"
  }
  
  template<typename VW, bool Aligned>
  void FastEvolutionProxy::compute(VW target) const{
    auto& suv3=target; //alias for the name expected by generated code
    size_t offset=suv1.GetEvolveBufferSize()/2;
    const double* CX=coefficients;
    const double* SX=coefficients+offset;
#include "../SU_inc/FastEvolutionSelect.txt"
  }
  
  template<typename VW, bool Aligned>
  void AdditionProxy::compute(VW target) const{
    for(unsigned int i=0; i<target.dim*target.dim; i++)
      target.components[i] += suv1.components[i] + suv2.components[i];
  }
  
  template<typename VW, bool Aligned>
  void SubtractionProxy::compute(VW target) const{
    for(unsigned int i=0; i<target.dim*target.dim; i++)
      target.components[i] += suv1.components[i] - suv2.components[i];
  }
  
  template<typename VW, bool Aligned>
  void NegationProxy::compute(VW target) const{
    for(unsigned int i=0; i<target.dim*target.dim; i++)
      target.components[i] += -suv1.components[i];
  }
  
  template<typename VW, bool Aligned>
  void MultiplicationProxy::compute(VW target) const{
    //the rest of this is just a convoluted, and hopefully fast way of doing this:
    //for(unsigned int i=0; i<target.dim*target.dim; i++)
    //  target.components[i] += a*suv1.components[i];
    
    double* suv1c=suv1.components;
    auto size=suv1.size;
    SQUIDS_COMPILER_ASSUME(size>=1);
    //if the number of components is odd, handle the first and adjust pointers
    //to refer to the remainder.
    if(size%2==1){
      target.components[0] += a*suv1c[0];
      if(size==1)
        return;
      suv1c++;
      target.components.components++;
      size--;
    }
    //All even squares are multiples of four to begin with, and all odd squares
    //greater than one are one more than a multiple of four, so now size%4==0.
    SQUIDS_COMPILER_ASSUME(size%4==0);
    
#if SQUIDS_USE_VECTOR_EXTENSIONS
    //if the data is aligned, take the fast path
    if(Aligned || (!((intptr_t)suv1c%32) && !((intptr_t)target.components.components%32))){
      double_vector4 f={a,a,a,a};
      for(; size>0; size-=4, suv1c+=4, target.components.components+=4){
        //try to convince the compiler to use aligned loads and stores
        SQUIDS_POINTER_IS_ALIGNED(suv1c,32);
        SQUIDS_POINTER_IS_ALIGNED(target.components.components,32);
        target.apply4(f * *(double_vector4*)suv1c);
      }
    }
    else{ //slow path for unaligned data
#endif //SQUIDS_USE_VECTOR_EXTENSIONS
      for(unsigned int i=0; i<size; i++)
        target.components[i] += a*suv1c[i];
#if SQUIDS_USE_VECTOR_EXTENSIONS
    }
#endif //SQUIDS_USE_VECTOR_EXTENSIONS
  }
  
  template<typename VW, bool Aligned>
  void iCommutatorProxy::compute(VW suv_new) const{
    suv_new.components[0]+=0;
#include "../SU_inc/iCommutatorSelect.txt"
  }
  
  template<typename VW, bool Aligned>
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