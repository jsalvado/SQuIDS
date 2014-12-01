#include <iostream>
#include <algorithm>
#include "SQuIDS/SUNalg.h"
#include "alloc_counting.h"

//If we construct a new SU_vector from the result of an elementwise operation
//in which the first SU_vector operand is movable, the newly constructed
//object should take the movable object's backing storage, rather than
//allocating its own.

int main(){
	{
		SU_vector v(6);
		v.SetAllComponents(1);
		alloc_counting::reset_allocation_counters();
		SU_vector v2=2.*std::move(v);
		size_t mem_allocated=alloc_counting::mem_allocated/sizeof(double);
		std::cout << mem_allocated << std::endl;
		std::cout << v2.GetComponents()[0] << std::endl;
	}
	
	{
		SU_vector v1(6), v2(6);
		v1.SetAllComponents(1);
		v2.SetAllComponents(2);
		alloc_counting::reset_allocation_counters();
		SU_vector v3=std::move(v1)+std::move(v2);
		size_t mem_allocated=alloc_counting::mem_allocated/sizeof(double);
		std::cout << mem_allocated << std::endl;
		std::cout << v3.GetComponents()[0] << std::endl;
	}
	
	{
		SU_vector v1(6), v2(6);
		v1.SetAllComponents(1);
		v2.SetAllComponents(2);
		alloc_counting::reset_allocation_counters();
		SU_vector v3=std::move(v1)-std::move(v2);
		size_t mem_allocated=alloc_counting::mem_allocated/sizeof(double);
		std::cout << mem_allocated << std::endl;
		std::cout << v3.GetComponents()[0] << std::endl;
	}
	
	{ //this optimization does not apply for operations which mix components
		SU_vector v1(6), v2(6);
		v1.SetAllComponents(1);
		v2.SetAllComponents(2);
		alloc_counting::reset_allocation_counters();
		SU_vector v3=ACommutator(std::move(v1),std::move(v2));
		size_t mem_allocated=alloc_counting::mem_allocated/sizeof(double);
		std::cout << mem_allocated << std::endl;
	}
}
