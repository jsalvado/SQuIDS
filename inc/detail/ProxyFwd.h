#ifndef SQUIDS_DETAIL_PROXYFWD_H
#define SQUIDS_DETAIL_PROXYFWD_H

#define REQUIRE_EVALUATION_PROXY_CORE typename std::enable_if<std::is_base_of<detail::EvaluationProxy<ProxyType>,ProxyType>::value>
#define REQUIRE_EVALUATION_PROXY_TPARAM typename= REQUIRE_EVALUATION_PROXY_CORE
#define REQUIRE_EVALUATION_PROXY_FPARAM REQUIRE_EVALUATION_PROXY_CORE ::type* =nullptr

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