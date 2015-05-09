#include <iostream>
#include <SQuIDs/SUNalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>

// This tests that conversion between GSL objects and SU_vectors
// works both ways.

double DC(double x){
  if (fabs(x) < 1.0e-15)
    return 0.;
  else
    return x;
}

int main(){
  using squids::SU_vector;
  
  for(unsigned int dim=2; dim<=6; dim++){
    gsl_matrix_complex * matrix_1 = gsl_matrix_complex_alloc(dim,dim);
    gsl_matrix_complex * matrix_2 = gsl_matrix_complex_alloc(dim,dim);
    
    //make some hermitian matrices
    for(int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++){
        gsl_matrix_complex_set(matrix_1,i,j,(gsl_complex){{static_cast<double>(i+j),static_cast<double>(i-j)}});
        gsl_matrix_complex_set(matrix_2,i,j,(gsl_complex){{static_cast<double>(i*j),static_cast<double>(2*(j-i))}});
      }
    }
    
    //convert them to SU vectors
    SU_vector m1(matrix_1);
    SU_vector m2(matrix_2);
    
    //compute the sum in the SU representation and convert back to a matrix
    auto matrix_sum = SU_vector(m1+m2).GetGSLMatrix();
    
    //print the results
    for(unsigned int i = 0; i < dim; i++){
      for(unsigned int j = 0; j < dim; j++){
        gsl_complex x = gsl_matrix_complex_get(matrix_sum.get(),i,j);
        std::cout << DC(GSL_REAL(x)) << "," << DC(GSL_IMAG(x)) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
	
    gsl_matrix_complex_free(matrix_1);
    gsl_matrix_complex_free(matrix_2);
  }
  
  return 0;
}
