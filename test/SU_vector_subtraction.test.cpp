#include <iostream>
#include <algorithm>
#include <SQuIDs/SUNalg.h>
#include "alloc_counting.h"

//Subtraction is simple, so while we're here also thoroughly test arithmetic optimizations
//Things to test:
//subtracting SU_vectors without optimization
//subtracting with fused assignment
//	assignment to vector of correct size
//		owned storage -> no allocation
//		external storage -> no allocation
//	assignment to vector of incorrect size
//		no storage -> resize
//		owned storage -> resize
//		external storage -> exception
//	aliased assignment
//		same vectors
//		vectors using the same external storage
//subtracting an SU_vector and a proxy
//constructing or assigning a vector from a proxy which refences an r-value

using squids::SU_vector;

void check_all_components_equal(const SU_vector& v, double expected){
	auto components=v.GetComponents();
	std::cout << "components " <<
	(std::all_of(components.begin(),components.end(),[=](double c){ return(c==expected); })?
	 "are":"are not") << " correctly set\n";
}

void test_pure_subtract(unsigned int dim, bool matchingSizes){
	const double d1=6.36, d2=7.75, diff=d1-d2;
	SU_vector v1(dim), v2(matchingSizes?dim:(dim>2?dim-1:dim+1));
	std::cout << "Subtracting vectors with sizes " << v1.Dim() << " and " << v2.Dim() << '\n';
	v1.SetAllComponents(d1);
	v2.SetAllComponents(d2);
	try{
		SU_vector v3(v1-v2);
		check_all_components_equal(v3,diff);
	} catch(std::runtime_error& err){
		std::cout << "Exception: " << err.what() << '\n';
	}
}

void test_fused_assign_subtract(unsigned int dim, SU_vector& dest, const std::string& dStoreType){
	const double d1=4.97, d2=2.14, diff=d1-d2;
	SU_vector v1(dim), v2(dim);
	std::cout << "Assigning difference of vectors with sizes " << v1.Dim() << " and " << v2.Dim()
		<< " to vector with size " << dest.Dim() << " (" << dStoreType << ")\n";
	v1.SetAllComponents(d1);
	v2.SetAllComponents(d2);
	try{
		alloc_counting::reset_allocation_counters();
		dest=v1-v2;
		auto allocated=alloc_counting::mem_allocated;
		std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
		check_all_components_equal(dest,diff);
	} catch(std::runtime_error& err){
		std::cout << "Exception: " << err.what() << '\n';
	}
}

template<typename Callable>
void make_dest_none(Callable test){
	SU_vector dest;
	test(dest,"none");
}

template<typename Callable>
void make_dest_internal(unsigned int dim, Callable test){
	SU_vector dest(dim);
	test(dest,"internal");
}

template<typename Callable>
void make_dest_external(unsigned int dim, Callable test){
	const size_t size=dim*dim;
	std::unique_ptr<double[]> buffer(new double[size]);
	SU_vector dest(dim,buffer.get());
	test(dest,"external");
	//verify that memory has not shifted
	buffer[0]=22.5;
	std::cout << "Destination memory has not moved: " << (dest[0]==buffer[0]) << '\n';
}

//aliasing requires allocation!
void test_fused_assign_subtract_aliased_vectors(unsigned int dim){
	const double d1=5.28, d2=8.06, d3=1.83, diff1=d1-d2, diff2=d2-d3;
	SU_vector v1(dim), v2(dim), v3(dim);
	std::cout << "Aliased vector fused assign-subtract with size " << dim << '\n';
	v1.SetAllComponents(d1);
	v2.SetAllComponents(d2);
	v3.SetAllComponents(d3);
	try{
		alloc_counting::reset_allocation_counters();
		v1=v1-v2;
		auto allocated=alloc_counting::mem_allocated;
		std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
		check_all_components_equal(v1,diff1);
	} catch(std::runtime_error& err){
		std::cout << "Exception: " << err.what() << '\n';
	}
	try{
		alloc_counting::reset_allocation_counters();
		v3=v2-v3;
		auto allocated=alloc_counting::mem_allocated;
		std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
		check_all_components_equal(v3,diff2);
	} catch(std::runtime_error& err){
		std::cout << "Exception: " << err.what() << '\n';
	}
}

//aliasing requires allocation!
void test_fused_assign_subtract_aliased_storage(unsigned int dim){
	const size_t size=dim*dim;
	std::unique_ptr<double[]> b1(new double[size]), b2(new double[size]), b3(new double[size]);
	const double d1=3.70, d2=2.29, diff=d1-d2;
	std::cout << "Aliased storage fused assign-subtract with size " << dim << '\n';
	try{
		SU_vector v1(dim,b1.get()), v2(dim,b2.get()), v3(dim,b1.get());
		v1.SetAllComponents(d1);
		v2.SetAllComponents(d2);
		alloc_counting::reset_allocation_counters();
		v3=v1-v2;
		auto allocated=alloc_counting::mem_allocated;
		std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
		check_all_components_equal(v3,diff);
	} catch(std::runtime_error& err){
		std::cout << "Exception: " << err.what() << '\n';
	}
	try{
		SU_vector v1(dim,b1.get()), v2(dim,b2.get()), v3(dim,b2.get());
		v1.SetAllComponents(d1);
		v2.SetAllComponents(d2);
		alloc_counting::reset_allocation_counters();
		v3=v1-v2;
		auto allocated=alloc_counting::mem_allocated;
		std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
		check_all_components_equal(v3,diff);
	} catch(std::runtime_error& err){
		std::cout << "Exception: " << err.what() << '\n';
	}
}

//this test covers both performing another arithmetic operation on proxy
//(which forces the first proxy to be evaulated into a temporary) and
//also the construction of a new vector  and the assignment to a vector
//without storage allocated from a proxy which references an r-value
//from which memory can be taken
void test_subtract_vector_proxy(unsigned int dim){
	const double d1=7.31, d2=7.45, d3=2.46, diff=(d1-d2)-d3;
	SU_vector v1(dim), v2(dim), v3(dim);
	std::cout << "Subtracting 3 vectors with size " << dim << '\n';
	v1.SetAllComponents(d1);
	v2.SetAllComponents(d2);
	v3.SetAllComponents(d3);
	try{
		//test constructing a new vector
		alloc_counting::reset_allocation_counters();
		SU_vector v4((v1-v2)-v3);
		//a temporary SU_vector must be constructed for v1-v2
		//but this will then be an r-value whose memory v4 can take
		auto allocated=alloc_counting::mem_allocated;
		std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
		check_all_components_equal(v4,diff);
		
		//test assigning to a default constructed vector
		alloc_counting::reset_allocation_counters();
		SU_vector v5; //default construct; no backing storage yet
		v5=(v1-v2)-v3;
		allocated=alloc_counting::mem_allocated;
		std::cout << allocated/sizeof(double) << " entries allocated" << '\n';
		check_all_components_equal(v5,diff);
	} catch(std::runtime_error& err){
		std::cout << "Exception: " << err.what() << '\n';
	}
}

int main(){
	alloc_counting::pattern_fill_allocs=true;
	alloc_counting::alloc_fill_pattern=0xFF;
	
	for(unsigned int i=2; i<=6; i++){
		test_pure_subtract(i,true);
		test_pure_subtract(i,false);
		
		if(i<5){
		auto test_fused_assign_subtract=[i](SU_vector& dest, const std::string& dStoreType){
			::test_fused_assign_subtract(i,dest,dStoreType);
		};
		make_dest_none(test_fused_assign_subtract);
		make_dest_internal(i,test_fused_assign_subtract);
		make_dest_external(i,test_fused_assign_subtract);
		
		make_dest_internal((i>2?i-1:i+1),test_fused_assign_subtract);
		make_dest_external((i>2?i-1:i+1),test_fused_assign_subtract);
		}
		
		test_fused_assign_subtract_aliased_vectors(i);
		test_fused_assign_subtract_aliased_storage(i);
		
		test_subtract_vector_proxy(i);
		
		std::cout << '\n';
	}
}