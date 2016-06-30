#include <iostream>
#include <algorithm>
#include <SQuIDS/SUNalg.h>
#include "alloc_counting.h"

using squids::SU_vector;

void exercise_constructor(unsigned int dim){
	std::cout << "Dimension " << dim << '\n';
	const size_t size=dim*dim;
	const double initVal=3.14;
	
	SU_vector source1(dim);
	source1.SetAllComponents(initVal);
	
	std::unique_ptr<double[]> buffer(new double[size]);
	for(unsigned int i=0; i<size; i++)
		buffer[i]=initVal;
	SU_vector source2(dim,buffer.get());
	
	std::cout << "Internal storage\n";
	alloc_counting::reset_allocation_counters();
	SU_vector::clear_mem_cache();
	SU_vector dest1(source1);
	auto allocated=alloc_counting::mem_allocated;
	std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
	//check the number of components
	auto components=dest1.GetComponents();
	std::cout << components.size() << " components stored" << '\n';
	//check that all components are initialized to the correct value
	std::cout << "components " <<
		(std::all_of(components.begin(),components.end(),[=](double c){ return(c==initVal); })?
		 "are":"are not") << " correctly set\n";
	//check that memory is not aliased
	source1[0]=0;
	std::cout << "Memory aliasing: " << (dest1[0]==source1[0]) << '\n';
	
	std::cout << "External storage\n";
	alloc_counting::reset_allocation_counters();
	SU_vector::clear_mem_cache();
	SU_vector dest2(source2);
	allocated=alloc_counting::mem_allocated;
	std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
	//check the number of components
	components=dest2.GetComponents();
	std::cout << components.size() << " components stored" << '\n';
	//check that all components are initialized to the correct value
	std::cout << "components " <<
		(std::all_of(components.begin(),components.end(),[=](double c){ return(c==initVal); })?
		 "are":"are not") << " correctly set\n";
	//check that memory is not aliased
	source2[0]=0;
	std::cout << "Memory aliasing: " << (dest2[0]==source2[0]) << '\n';
	std::cout << '\n';
}

int main(){
	alloc_counting::pattern_fill_allocs=true;
	alloc_counting::alloc_fill_pattern=0xFF;
	for(unsigned int i=2; i<=6; i++)
		exercise_constructor(i);
}
