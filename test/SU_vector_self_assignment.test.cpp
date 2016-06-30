#include <iostream>
#include <algorithm>
#include <SQuIDs/SUNalg.h>
#include "alloc_counting.h"

using squids::SU_vector;

void check_all_components_equal(const SU_vector& v, double expected){
	const auto& components=v.GetComponents();
	if(!std::all_of(components.begin(),components.end(),[=](double c){ return(c==expected); }))
		std::cout << "components are not correctly set \n";
}

void check_allocation(const char* msg){
	size_t allocated=alloc_counting::mem_allocated;
	if(allocated!=0)
		std::cout << "Memory unnecesarily allocated by " << msg << '\n';
}

int main(){
	alloc_counting::pattern_fill_allocs=true;
	alloc_counting::alloc_fill_pattern=0xFF;
	
	for(unsigned int i=1; i<=6; i++){
		SU_vector v(i);
		
		v.SetAllComponents(5);
		alloc_counting::reset_allocation_counters();
		v=v;
		check_allocation("self assignment");
		check_all_components_equal(v,5);
		
		alloc_counting::reset_allocation_counters();
		v+=v;
		check_allocation("self increment");
		check_all_components_equal(v,10);
		
		alloc_counting::reset_allocation_counters();
		v-=v;
		check_allocation("self decrement");
		check_all_components_equal(v,0);
		
		v.SetAllComponents(6);
		alloc_counting::reset_allocation_counters();
		v*=4;
		size_t allocated=alloc_counting::mem_allocated;
		check_allocation("multiplication");
		check_all_components_equal(v,24);
		
		alloc_counting::reset_allocation_counters();
		v/=3;
		check_allocation("division");
		check_all_components_equal(v,8);
		
		alloc_counting::reset_allocation_counters();
		v=-v;
		check_allocation("self negation");
		check_all_components_equal(v,-8);
	}
}