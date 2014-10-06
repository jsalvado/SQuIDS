#include <iostream>
#include <SQuIDs/SUNalg.h>
#include "alloc_counting.h"

int main(){
	SU_vector v;
	auto allocated=alloc_counting::mem_allocated; //should be 0
	std::cout << allocated << " bytes allocated\n";
	std::cout << "vector dimension: " << v.Dim() << '\n';
	return(0);
}