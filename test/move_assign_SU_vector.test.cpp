#include <iostream>
#include <SQuIDs/SUNalg.h>
#include "alloc_counting.h"

//cases:
// source storage: none/internal/external
// dest storage: none/internal/external
// size: matching/nonmatching

void exercise_assignment(SU_vector& source, SU_vector& dest, double expectedVal){
	try{
		alloc_counting::reset_allocation_counters();
		dest=std::move(source);
		auto allocated=alloc_counting::mem_allocated;
		std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
		//check the number of components
		auto components=dest.GetComponents();
		std::cout << components.size() << " components stored" << '\n';
		//check that all components are initialized to the correct value
		std::cout << "components " <<
		(std::all_of(components.begin(),components.end(),[=](double c){ return(c==expectedVal); })?
		 "are":"are not") << " correctly set\n";
		std::cout << "source dimension: " << source.Dim() << '\n';
		//check that memory is not aliased
		if(source.Dim()){
			source[0]=5;
			std::cout << "Memory aliasing: " << (dest[0]==source[0]) << '\n';
		}
	}catch(std::runtime_error& e){
		std::cout << "Assignment failed: " << e.what() << '\n';
	}
}

void make_dest_none(unsigned int dim, SU_vector& source, double expectedVal, const std::string& sstr){
	std::cout << " -> no storage, " << sstr << "\n";
	SU_vector dest;
	exercise_assignment(source,dest,expectedVal);
}

void make_dest_internal(unsigned int dim, SU_vector& source, double expectedVal, const std::string& sstr){
	std::cout << " -> internal storage, " << sstr << "\n";
	SU_vector dest(dim);
	exercise_assignment(source,dest,expectedVal);
}

void make_dest_external(unsigned int dim, SU_vector& source, double expectedVal, const std::string& sstr){
	std::cout << " -> external storage, " << sstr << "\n";
	const size_t size=dim*dim;
	std::unique_ptr<double[]> buffer(new double[size]);
	SU_vector dest(dim,buffer.get());
	exercise_assignment(source,dest,expectedVal);
	//verify that memory has not shifted
	buffer[0]=22.5;
	std::cout << "Destination memeory has not moved: " << (dest[0]==buffer[0]) << '\n';
}

void make_src_none(unsigned int dim, bool matching_sizes, void(*make_dest)(unsigned int,SU_vector&,double,const std::string&)){
	std::cout << "no storage";
	SU_vector source;
	make_dest(dim,source,0,(matching_sizes?"matching sizes":"non-matching sizes"));
}

void make_src_internal(unsigned int dim, bool matching_sizes, void(*make_dest)(unsigned int,SU_vector&,double,const std::string&)){
	std::cout << "internal storage";
	const unsigned int src_dim=(matching_sizes?dim:(dim>2?dim-1:dim+1));
	const size_t src_size=src_dim*src_dim;
	const double initVal=.5772;
	
	SU_vector source(src_dim);
	source.SetAllComponents(initVal);
	
	make_dest(dim,source,initVal,(matching_sizes?"matching sizes":"non-matching sizes"));
}

void make_src_external(unsigned int dim, bool matching_sizes, void(*make_dest)(unsigned int,SU_vector&,double,const std::string&)){
	std::cout << "external storage";
	const unsigned int src_dim=(matching_sizes?dim:(dim>2?dim-1:dim+1));
	const size_t src_size=src_dim*src_dim;
	const double initVal=.5772;
	
	std::unique_ptr<double[]> buffer(new double[src_size]);
	SU_vector source(src_dim,buffer.get());
	source.SetAllComponents(initVal);
	
	make_dest(dim,source,initVal,(matching_sizes?"matching sizes":"non-matching sizes"));
}

int main(){
	alloc_counting::pattern_fill_allocs=true;
	alloc_counting::alloc_fill_pattern=0xFF;
	for(unsigned int i=2; i<=2; i++){
		std::cout << "Dimension " << i << '\n';
		make_src_none(i,true,make_dest_none);
		make_src_internal(i,true, make_dest_none);
		make_src_internal(i,true, make_dest_internal);
		make_src_internal(i,false,make_dest_internal);
		make_src_internal(i,true, make_dest_external);
		make_src_internal(i,false,make_dest_external);
		make_src_external(i,true, make_dest_none);
		make_src_external(i,true, make_dest_internal);
		make_src_external(i,false,make_dest_internal);
		make_src_external(i,true, make_dest_external);
		make_src_external(i,false,make_dest_external);
		std::cout << '\n';
	}
}