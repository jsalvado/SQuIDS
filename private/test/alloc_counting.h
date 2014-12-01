// Simplistic shims for counting total memory allocated

#include <new>
#include <cstdlib>

namespace alloc_counting{
	size_t allocations=0;
	size_t mem_allocated=0;
	//if pattern_fill_allocs all newly allocated memory will be filled with alloc_fill_pattern
	bool pattern_fill_allocs=false;
	unsigned char alloc_fill_pattern=0;
	
	void reset_allocation_counters(){
		using namespace alloc_counting;
		allocations=0;
		mem_allocated=0;
	}
}

void* operator new(std::size_t size){
	using namespace alloc_counting;
	void* p=malloc(size);
	if(!p)
		throw std::bad_alloc();
	allocations++;
	mem_allocated+=size;
	if(pattern_fill_allocs)
		std::fill((unsigned char*)p,(unsigned char*)p+size,alloc_fill_pattern);
	return(p);
}

void* operator new[](std::size_t size){
	using namespace alloc_counting;
	void* p=malloc(size);
	if(!p)
		throw std::bad_alloc();
	allocations++;
	mem_allocated+=size;
	if(pattern_fill_allocs)
		std::fill((unsigned char*)p,(unsigned char*)p+size,alloc_fill_pattern);
	return(p);
}

void operator delete(void* p) noexcept{
	free(p);
}

void operator delete[](void* p) noexcept{
	free(p);
}