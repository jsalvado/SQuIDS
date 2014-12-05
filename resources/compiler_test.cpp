#include <memory>

template<typename T>
void sink(T&& t){}

int main(){
	std::unique_ptr<int> my_int(new int);
	sink(std::move(my_int));
}
