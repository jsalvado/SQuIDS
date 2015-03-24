#include <iostream>
#include <algorithm>
#include <SQuIDs/SUNalg.h>
#include "alloc_counting.h"

//cases:
// storage: none/internal/external
// size: matching/nonmatching

using squids::SU_vector;

void exercise_assignment(SU_vector& source, SU_vector& dest, double expectedVal){
	try{
		alloc_counting::reset_allocation_counters();
		dest=source;
		auto allocated=alloc_counting::mem_allocated;
		std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
		//check the number of components
		auto components=dest.GetComponents();
		std::cout << components.size() << " components stored" << '\n';
		//check that all components are initialized to the correct value
		std::cout << "components " <<
		(std::all_of(components.begin(),components.end(),[=](double c){ return(c==expectedVal); })?
		 "are":"are not") << " correctly set\n";
		//check that memory is not aliased
		source[0]=0;
		std::cout << "Memory aliasing: " << (dest[0]==source[0]) << '\n';
	}catch(std::runtime_error& e){
		std::cout << "Assignment failed: " << e.what() << '\n';
	}
}

void exercise_assignment_none(unsigned int dim){
	std::cout << "No storage\n";
	const size_t size=dim*dim;
	const double initVal=2.718;
	
	SU_vector dest;
	SU_vector source(dim);
	source.SetAllComponents(initVal);
	exercise_assignment(source,dest,initVal);
}

void exercise_assignment_internal(unsigned int dim, bool matching_size){
	std::cout << "Internal storage, " << (matching_size?"matching":"non-matching") << " size\n";
	const size_t size=dim*dim;
	const double initVal=2.718;
	
	const unsigned int src_dim=(matching_size?dim:(dim>2?dim-1:dim+1));
	const size_t src_size=src_dim*src_dim;
	
	SU_vector dest(dim);
	SU_vector source(src_dim);
	source.SetAllComponents(initVal);
	exercise_assignment(source,dest,initVal);
}

void exercise_assignment_external(unsigned int dim, bool matching_size){
	std::cout << "External storage, " << (matching_size?"matching":"non-matching") << " size\n";
	const size_t size=dim*dim;
	const double initVal=2.718;
	
	const unsigned int src_dim=(matching_size?dim:(dim>2?dim-1:dim+1));
	const size_t src_size=src_dim*src_dim;
	
	std::unique_ptr<double[]> buffer(new double[size]);
	SU_vector dest(dim,buffer.get());
	SU_vector source(src_dim);
	source.SetAllComponents(initVal);
	exercise_assignment(source,dest,initVal);
}

int main(){
	alloc_counting::pattern_fill_allocs=true;
	alloc_counting::alloc_fill_pattern=0xFF;
	for(unsigned int i=2; i<=6; i++){
		std::cout << "Dimension " << i << '\n';
		exercise_assignment_none(i);
		exercise_assignment_internal(i,true);
		exercise_assignment_internal(i,false);
		exercise_assignment_external(i,true);
		exercise_assignment_external(i,false);
		std::cout << '\n';
	}
}