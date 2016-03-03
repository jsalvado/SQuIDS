#include <cmath>
#include <complex>
#include <iostream>
#include <SQuIDS/SUNalg.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

int main(){
  using namespace squids;

  const double pi=4*atan(1);
  const double base_tol=2e-16; //~machine epsilon
  const auto unit=gsl_complex_rect(1,0);
  const auto zero=gsl_complex_rect(0,0);
  auto from_gsl=[](const gsl_complex& cg)->std::complex<double>{
    return(std::complex<double>(GSL_REAL(cg),GSL_IMAG(cg)));
  };
  Const transform;

  for(unsigned int dim=2; dim<=SQUIDS_MAX_HILBERT_DIM; dim++){
    //make sure the transformation between bases is full of non-trivial angles
    for(unsigned int i=0; i<dim-1; i++){
      double angle=2.*pi*(i+1.)/(dim+1.);
      transform.SetMixingAngle(i,dim-1,angle);
      transform.SetPhase(i,dim-1,-angle);
      transform.SetEnergyDifference(dim-1,1);
    }
    
    auto U=transform.GetTransformationMatrix(dim);
    
    gsl_matrix_complex* temp = gsl_matrix_complex_alloc(dim,dim);
    gsl_matrix_complex* result = gsl_matrix_complex_alloc(dim,dim);
    
    //test transformation for every generator
    for(unsigned int g=0; g<(dim*dim); g++){
      SU_vector v=SU_vector::Generator(dim,g);
      auto vg=v.GetGSLMatrix();
      
      //U^dagger * vg -> temp
      gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,unit,U.get(),vg.get(),zero,temp);
      //temp * U -> result
      gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,unit,temp,U.get(),zero,result);
      
      //do the transformation natively
      v.RotateToB1(transform);
      //fetch the reference result
      auto vgt=v.GetGSLMatrix();
		
      for(size_t i=0; i<dim; i++){
        for(size_t j=0; j<dim; j++){
          auto diff=from_gsl(gsl_matrix_complex_get(result,i,j))-from_gsl(gsl_matrix_complex_get(vgt.get(),i,j));
          if(std::abs(diff.real())>2*dim*base_tol)
            std::cout << "Dimension " << dim << " generator " << g <<
            " component " << i << ',' << j << "(real) has error " << diff.real() << std::endl;
          if(std::abs(diff.imag())>2*dim*base_tol)
            std::cout << "Dimension " << dim << " generator " << g <<
            " component " << i << ',' << j << "(imag) has error " << diff.imag() << std::endl;
        }
      }
    }
    gsl_matrix_complex_free(temp);
    gsl_matrix_complex_free(result);
  }
}
