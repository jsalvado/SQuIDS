#include <cmath>
#include <iostream>
#include <SQuIDs/SUNalg.h>

int main(){
	using namespace squids;

	const double pi=4*atan(1);
	const double base_tol=2e-16; //~machine epsilon
	Const transform;

	for(unsigned int dim=2; dim<=SQUIDS_MAX_HILBERT_DIM; dim++){
		//make sure the transformation between bases is full of non-trivial angles
		for(unsigned int i=0; i<dim-1; i++){
			double angle=2.*pi*(i+1.)/(dim+1.);
			transform.SetMixingAngle(i,dim-1,angle);
			transform.SetPhase(i,dim-1,-angle);
		}
		//test transformation for every generator
		for(unsigned int g=0; g<(dim*dim); g++){
			//get the next generator and make a copy
			SU_vector s=SU_vector::Generator(dim,g);
			SU_vector to=s;
			// rotate s using RotateToB1
			s.RotateToB1(transform);
			// rotate t using gsl_unitary_rotation
			auto U = transform.GetTransformationMatrix(dim);
			SU_vector tr = to.Rotate(U.get());
			//check that the result matches the original closely
			auto sc=s.GetComponents();
			auto tc=tr.GetComponents();
			//as the dimension increases more operations are necessary in the
			//rotation, so we allow for a slightly more error to accumulate
			for(unsigned int i=0; i<(dim*dim); i++){
				if(std::abs(sc[i]-tc[i])>2*dim*base_tol)
					std::cout << "Dimension " << dim << " generator " << g <<
					" component " << i << " has error " << (sc[i]-tc[i]) << std::endl;
			}
		}
	}
}
