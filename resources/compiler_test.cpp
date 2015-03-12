#include <memory>

//need r-value references
//need defaulted template arguments for function templates
template<typename T, typename Dummy=void>
void sink(T&& t){}

struct type{
	int a,b;
	
	//need constexpr
	constexpr static bool flag=true;
	
	type():a(1),b(2){}
	//need delegating constructors
	type(int c):type(){
		a+=c;
		b-=c;
	}
	type(int a, int b):a(a),b(b){}
	
	//need defaulted functions
	~type()=default;
	
	//need ref qualified functions/r-value references for *this
	int foo() const &{ return(a); }
	int foo() &&{ return(b); }
};

int main(){
	//need unique_ptr
	//need nullptr
	std::unique_ptr<int> my_int(nullptr);
	sink(std::move(my_int));
	//need initializer lists
	type t{5,6};
	//need auto
	auto i=1;
}
