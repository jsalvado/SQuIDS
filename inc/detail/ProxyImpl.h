#ifndef SQUIDS_DETAIL_PROXYIMPL_H
#define SQUIDS_DETAIL_PROXYIMPL_H

namespace detail{
  template<typename Op>
  struct EvaluationProxy{
    const SU_vector& suv1;
    const SU_vector& suv2;
    
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
  };
  
  struct EvolutionProxy : public EvaluationProxy<EvolutionProxy>{
    double t;
    
    EvolutionProxy(const SU_vector& suv1,const SU_vector& suv2,double t):
    EvaluationProxy<EvolutionProxy>{suv1,suv2},t(t){}
    
    template<typename VW>
    void compute(VW target) const{
      auto& suv_new=target;
      switch (suv_new.dim){
        case 2:
#include "../SU_inc/EvolutionSU2.txt"
          break;
        case 3:
#include "../SU_inc/EvolutionSU3.txt"
          break;
        case 4:
#include "../SU_inc/EvolutionSU4.txt"
          break;
        case 5:
#include "../SU_inc/EvolutionSU5.txt"
          break;
        case 6:
          throw std::runtime_error("SU(6) evolution support not yet implemented");
          break;
        default:
          throw std::runtime_error("SUEvolution :: Error : dim  ");
      };
    }
  };
  
  struct AdditionProxy : public EvaluationProxy<AdditionProxy>{
    AdditionProxy(const SU_vector& suv1,const SU_vector& suv2):
    EvaluationProxy<AdditionProxy>{suv1,suv2}{}
    
    template<typename VW>
    void compute(VW target) const{
      for(int i=0; i<target.dim*target.dim; i++)
        target.components[i] += suv1.components[i] + suv2.components[i];
    }
  };
  
  struct SubtractionProxy : public EvaluationProxy<SubtractionProxy>{
    SubtractionProxy(const SU_vector& suv1,const SU_vector& suv2):
    EvaluationProxy<SubtractionProxy>{suv1,suv2}{}
    
    template<typename VW>
    void compute(VW target) const{
      for(int i=0; i<target.dim*target.dim; i++)
        target.components[i] += suv1.components[i] - suv2.components[i];
    }
  };
  
  struct MultiplicationProxy : public EvaluationProxy<MultiplicationProxy>{
    double a;
    
    MultiplicationProxy(const SU_vector& suv1,double a):
    EvaluationProxy<MultiplicationProxy>{suv1,suv1},a(a){}
    
    template<typename VW>
    void compute(VW target) const{
      for(int i=0; i<target.dim*target.dim; i++)
        target.components[i] += a*suv1.components[i];
    }
  };
  
  struct iCommutatorProxy : public EvaluationProxy<iCommutatorProxy>{
    iCommutatorProxy(const SU_vector& suv1,const SU_vector& suv2):
    EvaluationProxy<iCommutatorProxy>{suv1,suv2}{}
    
    template<typename VW>
    void compute(VW target) const{
      auto& suv_new=target;
      suv_new.components[0]+=0;
      switch (suv_new.dim){
        case 2:
#include "../SU_inc/iConmutatorSU2.txt"
          break;
        case 3:
#include "../SU_inc/iConmutatorSU3.txt"
          break;
        case 4:
#include "../SU_inc/iConmutatorSU4.txt"
          break;
        case 5:
#include "../SU_inc/iConmutatorSU5.txt"
          break;
        case 6:
#include "../SU_inc/iConmutatorSU6.txt"
          break;
        default:
          throw std::runtime_error("SUiComutator :: Error.");
      }
    }
  };
  
  struct ACommutatorProxy : public EvaluationProxy<ACommutatorProxy>{
    ACommutatorProxy(const SU_vector& suv1,const SU_vector& suv2):
    EvaluationProxy<ACommutatorProxy>{suv1,suv2}{}
    
    template<typename VW>
    void compute(VW target) const{
      auto& suv_new=target;
      switch (suv_new.dim){
        case 2:
#include "../SU_inc/AnticonmutatorSU2.txt"
          break;
        case 3:
#include "../SU_inc/AnticonmutatorSU3.txt"
          break;
        case 4:
#include "../SU_inc/AnticonmutatorSU4.txt"
          break;
        case 5:
#include "../SU_inc/AnticonmutatorSU5.txt"
          break;
        case 6:
#include "../SU_inc/AnticonmutatorSU6.txt"
          break;
        default:
          throw std::runtime_error("SUAnticomutator :: Error.");
      }
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