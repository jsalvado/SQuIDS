#include <iostream>
#include <algorithm>
#include <SQuIDS/SUNalg.h>
#include "alloc_counting.h"

void exercise_constructor(unsigned int dim){
	using squids::SU_vector;
	
	const size_t size=dim*dim;
	const double initVal=3.14;
	std::unique_ptr<double[]> buffer(new double[size]);
	for(unsigned int i=0; i<size; i++)
		buffer[i]=initVal;
	alloc_counting::reset_allocation_counters();
	SU_vector v(dim, buffer.get());
	auto allocated=alloc_counting::mem_allocated;
	std::cout << allocated << '\n';
	//check the number of components
	auto components=v.GetComponents();
	std::cout << components.size() << '\n';
	//check that all components are initialized to zero
	std::cout << std::all_of(components.begin(),components.end(),
	                         [=](double c){ return(c==initVal); }) << '\n';
}

int main(){
	for(unsigned int i=2; i<=6; i++)
		exercise_constructor(i);
}
