#include <iostream>
#include <SQuIDs/SUNalg.h>
#include <iostream>


int main(){
  using squids::SU_vector;

  gsl_matrix_complex *mb=gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex *mb2=gsl_matrix_complex_alloc(2,2);

  gsl_complex zt00, zt01, zt10, zt11;
  GSL_SET_COMPLEX(&zt00, 1, 0); GSL_SET_COMPLEX(&zt01, 1, 1);
  GSL_SET_COMPLEX(&zt10, 1, -1); GSL_SET_COMPLEX(&zt11, 1, 0);

  gsl_matrix_complex_set(mb,0,0,zt00); gsl_matrix_complex_set(mb,0,1,zt01);
  gsl_matrix_complex_set(mb,1,0,zt10); gsl_matrix_complex_set(mb,1,1,zt11);


  SU_vector v3(mb);
  auto m3=v3.GetGSLMatrix();

  for (size_t i = 0; i < mb->size1; i++) {
    for (size_t j = 0; j < mb->size2; j++) {
      if( GSL_REAL(gsl_matrix_complex_get(mb, i, j)) -
	  GSL_REAL(gsl_matrix_complex_get(m3.get(), i, j))>1e-15){
	std::cout << " comp: "<< i << "  " << j << " RE:  " << GSL_REAL(gsl_matrix_complex_get(mb, i, j)) -
	  GSL_REAL(gsl_matrix_complex_get(m3.get(), i, j)) << std::endl;
      }

      if( GSL_IMAG(gsl_matrix_complex_get(mb, i, j)) -
	  GSL_IMAG(gsl_matrix_complex_get(m3.get(), i, j))>1e-15){
	std::cout <<  " comp: "<< i << "  " << j << " IM:  " << GSL_IMAG(gsl_matrix_complex_get(mb, i, j)) -
	  GSL_IMAG(gsl_matrix_complex_get(m3.get(), i, j)) << std::endl;
      }

    }
  }

  return(0);
}
