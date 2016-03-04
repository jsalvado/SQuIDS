#include <iostream>
#include <SQuIDs/SUNalg.h>

//The inner product of two SU_vectors is defined in terms of the SUTrace operation

int main(){
	using squids::SU_vector;
	for(unsigned int dim=2; dim<=SQUIDS_MAX_HILBERT_DIM; dim++){
		double prod=SU_vector::Generator(dim,0)*SU_vector::Generator(dim,0);
		if(prod!=dim)
			std::cout << "Failure: Product of two instances of generator 0 of dimension "
				<< dim << " should be " << dim << "; got " << prod << std::endl;
		
		for(unsigned int i=1; i<dim*dim; i++){
			prod=SU_vector::Generator(dim,i)*SU_vector::Generator(dim,i);
			if(prod!=2)
				std::cout << "Failure: Product of two instances of generator "
				<< i << " of dimension " << dim
				<< " should be 2; got " << prod << std::endl;
		}
		
		for(unsigned int i=0; i<dim*dim; i++){
			for(unsigned int j=0; j<dim*dim; j++){
				if(j==i) //this case was covered above
					continue;
				prod=SU_vector::Generator(dim,i)*SU_vector::Generator(dim,j);
				if(prod!=0)
					std::cout << "Failure: Product of generator " << i
					<< " with generator " << j << " of dimension " << dim
					<< " should be 0; got " << prod << std::endl;
			}
		}
		
		for(unsigned int i=0; i<dim*dim; i++){
			for(unsigned int j=0; j<dim*dim; j++){
				if(j==i) //we want pairs of different generators
					continue;
				prod=(SU_vector::Generator(dim,i)+SU_vector::Generator(dim,j))
					*(SU_vector::Generator(dim,i)+SU_vector::Generator(dim,j));
				double expected=(i==0 || j==0 ? dim : 2)+2;
				if(prod!=expected)
					std::cout << "Failure: Product of two instances of the sum of generator "
					<< i << " and generator " << j << " of dimension " << dim
					<< " should be " << expected << "; got " << prod << std::endl;
			}
		}
	}
}