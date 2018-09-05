#include <iostream>
#include <algorithm>
#include <SQuIDS/SUNalg.h>

using squids::SU_vector;

void check_all_components_equal(const SU_vector& v, double expected, std::string context){
	auto components=v.GetComponents();
	if(!std::all_of(components.begin(),components.end(),[=](double c){
		return(std::abs(c-expected)<std::abs(expected*std::numeric_limits<double>::epsilon()));
	  })){
		std::cout << "components are not correctly set from " << context << '\n';
		std::cout.precision(20);
		std::cout << "  expected value: " << expected << " actual values: ";
		for(auto component : components)
			std::cout << component << ' ';
		std::cout << "relative error: " << std::abs(components.front()-expected)/std::numeric_limits<double>::epsilon();
		std::cout << std::endl;
	}
}

//Test whether operations which change code paths depending on the alignment
//of the SU_vector storage always give correct results

int main(){
	double a=1.7, b=2.2, c=3.4;
	std::string context;
	for(unsigned int i=2; i<=6; i++){
		std::cout << "size " << i << std::endl;
		//make buffers with some wiggle room
		std::unique_ptr<double[]> store1(new double[i*i+2*(32/sizeof(double))]);
		std::unique_ptr<double[]> store2(new double[i*i+2*(32/sizeof(double))]);
		std::unique_ptr<double[]> store3(new double[i*i+2*(32/sizeof(double))]);
		//always start from 32 byte alignement to make the test deterministic
		size_t base1=(32-((intptr_t)store1.get())%32)/sizeof(double);
		size_t base2=(32-((intptr_t)store2.get())%32)/sizeof(double);
		size_t base3=(32-((intptr_t)store3.get())%32)/sizeof(double);
		//cycle through all possible alignments
		for(unsigned int offset=0; offset<(32/sizeof(double)); offset++){
			std::cout << " offset " << offset << std::endl;
			SU_vector v1(i,store1.get()+offset+base1);
			SU_vector v2(i,store2.get()+offset+base2);
			SU_vector v3(i,store3.get()+offset+base3);
			
			//scalar multiplcation
			context="scalar multiplication and assignment";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v2=v1*c;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,a*c,context);
			
			context="scalar multiplication and increment";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v2+=v1*c;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,b+a*c,context);
			
			context="scalar multiplication and decrement";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v2-=v1*c;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,b-(a*c),context);
			
			//negation
			context="negation and assignment";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v2=-v1;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,-a,context);
			
			context="negation and increment";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v2+=-v1;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,b-a,context);
			
			context="negation and decrement";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v2-=-v1;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,b+a,context);
			
			//addition
			context="addition and assignment";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v3.SetAllComponents(c);
			v3=v1+v2;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,b,context);
			check_all_components_equal(v3,a+b,context);
			
			context="addition and increment";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v3.SetAllComponents(c);
			v3+=v1+v2;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,b,context);
			check_all_components_equal(v3,a+b+c,context);
			
			context="addition and decrement";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v3.SetAllComponents(c);
			v3-=v1+v2;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,b,context);
			check_all_components_equal(v3,c-(a+b),context);
			
			//subtraction
			context="subtraction and assignment";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v3.SetAllComponents(c);
			v3=v1-v2;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,b,context);
			check_all_components_equal(v3,a-b,context);
			
			context="subtraction and increment";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v3.SetAllComponents(c);
			v3+=v1-v2;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,b,context);
			check_all_components_equal(v3,c+a-b,context);
			
			context="subtraction and decrement";
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			v3.SetAllComponents(c);
			v3-=v1-v2;
			check_all_components_equal(v1,a,context);
			check_all_components_equal(v2,b,context);
			check_all_components_equal(v3,c-(a-b),context);
			
			//trace
			double expected=a*b*(i+2*i*i-2);
			v1.SetAllComponents(a);
			v2.SetAllComponents(b);
			double trace=v1*v2;
			//different order of operations can cause some rounding:
			if(std::abs(trace-expected)>1e-10)
				std::cout << "trace result not correct: " << expected << ' ' << trace << "\n";
		}
	}
}
