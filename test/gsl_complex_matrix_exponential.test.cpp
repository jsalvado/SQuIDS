#include <iostream>
#include <iomanip>
#include <fstream>
#include <SQuIDs/SUNalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace squids;

void simple_print(const gsl_matrix_complex * matrix){
  for(unsigned int i=0; i < matrix->size1; i++){
    for(unsigned int j=0; j < matrix->size2; j++){
      auto c = gsl_matrix_complex_get(matrix,i,j);
      std::cout << GSL_REAL(c) << " " << GSL_IMAG(c) << " ";
    }
    std::cout << std::endl;
  }
}

int main(){
  unsigned int dim = 3;
  gsl_matrix_complex * M = gsl_matrix_complex_alloc(dim,dim);
  gsl_matrix_complex * eM = gsl_matrix_complex_alloc(dim,dim);

  // set printing precision
  std::cout << std::scientific;
  std::cout << std::setprecision(8);
  // first matrix
  gsl_matrix_complex_set(M,0,0,gsl_complex_rect(1,2));
  gsl_matrix_complex_set(M,0,1,gsl_complex_rect(3,-2));
  gsl_matrix_complex_set(M,0,2,gsl_complex_rect(5,-7));

  gsl_matrix_complex_set(M,1,0,gsl_complex_rect(3,9));
  gsl_matrix_complex_set(M,1,1,gsl_complex_rect(4,-2));
  gsl_matrix_complex_set(M,1,2,gsl_complex_rect(9,-7));

  gsl_matrix_complex_set(M,2,0,gsl_complex_rect(5,7));
  gsl_matrix_complex_set(M,2,1,gsl_complex_rect(1,-1));
  gsl_matrix_complex_set(M,2,2,gsl_complex_rect(5,-3));

  // calculating matrix exponential
  squids::math_detail::matrix_exponential(eM,M);
  // simple print matrix
  simple_print(eM);
  std::cout << std::endl;

  // second matrix
  gsl_matrix_complex_set(M,0,0,gsl_complex_rect(0,0.1));
  gsl_matrix_complex_set(M,0,1,gsl_complex_rect(1,20));
  gsl_matrix_complex_set(M,0,2,gsl_complex_rect(50,-1));

  gsl_matrix_complex_set(M,1,0,gsl_complex_rect(4,5));
  gsl_matrix_complex_set(M,1,1,gsl_complex_rect(0,0));
  gsl_matrix_complex_set(M,1,2,gsl_complex_rect(0.9,0.8));

  gsl_matrix_complex_set(M,2,0,gsl_complex_rect(50,70));
  gsl_matrix_complex_set(M,2,1,gsl_complex_rect(0.1,-0.1));
  gsl_matrix_complex_set(M,2,2,gsl_complex_rect(5,-5));

  squids::math_detail::matrix_exponential(eM,M);
  simple_print(eM);
  std::cout << std::endl;

  // third matrix
  gsl_matrix_complex_set(M,0,0,gsl_complex_rect(0.1,0.01));
  gsl_matrix_complex_set(M,0,1,gsl_complex_rect(0.02,0.2));
  gsl_matrix_complex_set(M,0,2,gsl_complex_rect(0.3,0.3));

  gsl_matrix_complex_set(M,1,0,gsl_complex_rect(-1,-1));
  gsl_matrix_complex_set(M,1,1,gsl_complex_rect(-5,1));
  gsl_matrix_complex_set(M,1,2,gsl_complex_rect(0.1,0.));

  gsl_matrix_complex_set(M,2,0,gsl_complex_rect(0,0));
  gsl_matrix_complex_set(M,2,1,gsl_complex_rect(0.7,0.3));
  gsl_matrix_complex_set(M,2,2,gsl_complex_rect(1,1));

  squids::math_detail::matrix_exponential(eM,M);
  simple_print(eM);
  std::cout << std::endl;

  //freeing matrices
  gsl_matrix_complex_free(M);
  gsl_matrix_complex_free(eM);

  return 0;
}
