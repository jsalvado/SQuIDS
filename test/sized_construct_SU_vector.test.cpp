#include <iostream>
#include <algorithm>
#include <SQuIDS/SUNalg.h>
#include "alloc_counting.h"

void exercise_constructor(unsigned int dim){
	squids::SU_vector v(dim);
	//check the number of components
	auto components=v.GetComponents();
	std::cout << components.size() << " components stored\n";
	//check that all components are initialized to zero
	std::cout << "components zero initialized: " <<
		std::all_of(components.begin(),components.end(),
	                [](double c){ return(c==0.0); }) << '\n';
}

int main(){
	alloc_counting::pattern_fill_allocs=true;
	alloc_counting::alloc_fill_pattern=0xFF;
	for(unsigned int i=2; i<=6; i++)
		exercise_constructor(i);
}
