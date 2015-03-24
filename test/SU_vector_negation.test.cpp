#include <iostream>
#include <algorithm>
#include <SQuIDs/SUNalg.h>
#include "alloc_counting.h"

using squids::SU_vector;

void check_all_components_equal(const SU_vector& v, double expected){
  auto components=v.GetComponents();
  std::cout << "components " <<
  (std::all_of(components.begin(),components.end(),[=](double c){ return(c==expected); })?
   "are":"are not") << " correctly set\n";
}

void test_self_assign(unsigned int dim){
  const double value=61.5;
  std::cout << "Inverting vector with dimension " << dim << std::endl;
  SU_vector v(dim);
  v.SetAllComponents(value);
  try{
    v=-v;
    check_all_components_equal(v,-value);
  } catch(std::runtime_error& err){
    std::cout << "Exception: " << err.what() << '\n';
  }
}

void test_copy_assign(unsigned int dim){
  const double value=4.5;
  std::cout << "Assigning negated vector with dimension " << dim << std::endl;
  SU_vector v1(dim),v2;
  v1.SetAllComponents(value);
  try{
    v2=-std::move(v1);
    check_all_components_equal(v2,-value);
  } catch(std::runtime_error& err){
    std::cout << "Exception: " << err.what() << '\n';
  }
}

void test_move_assign(unsigned int dim){
  const double value=-3.2;
  std::cout << "Move assigning negated vector with dimension " << dim << std::endl;
  alloc_counting::reset_allocation_counters();
  size_t allocated;
  SU_vector v1(dim),v2;
  v1.SetAllComponents(value);
  try{
    v2=-std::move(v1);
    allocated=alloc_counting::mem_allocated;
    check_all_components_equal(v2,-value);
  } catch(std::runtime_error& err){
    std::cout << "Exception: " << err.what() << '\n';
  }
  std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
}

void test_negate_proxy(unsigned int dim){
  const double d1=4.4, d2=-9.6, d3=-(d1+d2);
  std::cout << "Constructing vector from negation of a proxy with dimension " << dim << std::endl;
  SU_vector v1(dim), v2(dim);
  v1.SetAllComponents(d1);
  v2.SetAllComponents(d2);
  alloc_counting::reset_allocation_counters();
  size_t allocated;
  try{
    SU_vector v3=-(std::move(v1)+std::move(v2)); //suppose we no longer care about v1 and v2
    allocated=alloc_counting::mem_allocated;
    check_all_components_equal(v3,d3);
  } catch(std::runtime_error& err){
    std::cout << "Exception: " << err.what() << '\n';
  }
  std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
}

int main(){
  alloc_counting::pattern_fill_allocs=true;
  alloc_counting::alloc_fill_pattern=0xFF;
  
  for(unsigned int i=2; i<=6; i++){
    test_self_assign(i);
    test_copy_assign(i);
    test_move_assign(i);
    test_negate_proxy(i);
    std::cout << std::endl;
  }
}