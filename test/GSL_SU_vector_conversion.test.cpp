#include <iostream>
#include <SQuIDs/SUNalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>

// This test will proof that conversion between GSL objects and SU_vectors
// works both ways.

double DC(double x){
  if (fabs(x) < 1.0e-15)
    return 0.;
  else
    return x;
}

int main(){
  using squids::SU_vector;
  
  unsigned int dim = 3;
  gsl_matrix_complex * matrix_1 = gsl_matrix_complex_alloc(dim,dim);
  gsl_matrix_complex * matrix_2 = gsl_matrix_complex_alloc(dim,dim);

  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      gsl_matrix_complex_set(matrix_1,i,j,(gsl_complex){{static_cast<double>(i+j),static_cast<double>(i-j)}});
      gsl_matrix_complex_set(matrix_2,i,j,(gsl_complex){{static_cast<double>(i*j),static_cast<double>(j-i)}});
    }
  }

  SU_vector m1 = SU_vector(matrix_1);
  SU_vector m2 = SU_vector(matrix_2);

  gsl_matrix_complex_add(matrix_1,matrix_2);

  auto matrix_sum = SU_vector(m1+m2).GetGSLMatrix();

  for(unsigned int i = 0; i < dim; i++){
    for(unsigned int j = 0; j < dim; j++){
      gsl_complex x = gsl_matrix_complex_get(matrix_sum.get(),i,j);
      std::cout << DC(GSL_REAL(x)) << "," << DC(GSL_IMAG(x)) << " ";
    }
    std::cout << std::endl;
  }

  for(unsigned int i = 0; i < dim; i++){
    for(unsigned int j = 0; j < dim; j++){
      gsl_complex x = gsl_matrix_complex_get(matrix_1,i,j);
      std::cout << GSL_REAL(x) << "," << GSL_IMAG(x) << " ";
    }
    std::cout << std::endl;
  }
	
  gsl_matrix_complex_free(matrix_1);
  gsl_matrix_complex_free(matrix_2);

  return 0;
}
