#ifndef SQUIDS_DETAIL_PROXYFWD_H
#define SQUIDS_DETAIL_PROXYFWD_H

#define REQUIRE_EVALUATION_PROXY typename=typename std::enable_if<std::is_base_of<detail::EvaluationProxy<ProxyType>,ProxyType>::value>

namespace detail{
  
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
  
  //fwd decls for operation proxies
  template<typename Op>
  struct EvaluationProxy;
  struct EvolutionProxy;
  struct AdditionProxy;
  struct SubtractionProxy;
  struct MultiplicationProxy;
  struct iCommutatorProxy;
  struct ACommutatorProxy;
}

#endif //SQUIDS_DETAIL_PROXYFWD_H