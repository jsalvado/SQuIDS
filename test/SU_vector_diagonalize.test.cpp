#include <cmath>
#include <iostream>
#include <iomanip>
#include <SQuIDs/SUNalg.h>
#include <chrono>
#include <gsl/gsl_complex_math.h>


int main(){
	using namespace squids;

  SU_vector v(std::vector<double>{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9});

  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();
  // first we use the analitic formulaes
  auto eigen_system = v.GetEigenSystem();
  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();

  std::clock_t c_start_2 = std::clock();
  auto t_start_2 = std::chrono::high_resolution_clock::now();
  // then we use the GSL numerical solution
  unsigned int dim = v.Dim();
  gsl_vector * eigenvalues = gsl_vector_alloc(dim);
  gsl_matrix_complex * eigenvectors = gsl_matrix_complex_alloc(dim,dim);
  auto matrix=v.GetGSLMatrix();
  gsl_eigen_hermv_workspace * ws = gsl_eigen_hermv_alloc(dim);
  gsl_eigen_hermv(matrix.get(),eigenvalues,eigenvectors,ws);
  gsl_eigen_hermv_free(ws);
  gsl_eigen_hermv_sort(eigenvalues,eigenvectors,GSL_EIGEN_SORT_VAL_ASC);
  std::clock_t c_end_2 = std::clock();
  auto t_end_2 = std::chrono::high_resolution_clock::now();

  std::cout << std::setprecision(10);
  // checking eigenvalues agreement
  for(unsigned int i = 0; i < dim; i++){
    double eps = gsl_vector_get(eigen_system.first.get(),i) - gsl_vector_get(eigenvalues,i);
    if (eps > 1.0e-10)
      std::cout << gsl_vector_get(eigen_system.first.get(),i) << " " << gsl_vector_get(eigenvalues,i) << std::endl;
  }

  // checking eigenvectors
  // the eigenvectors can be off by a phase, but the absolute value must match
  for(unsigned int i = 0; i < eigen_system.second.get()->size1; i++){
    for(unsigned int j = 0; j < eigen_system.second.get()->size2; j++){
      double eps = gsl_complex_abs(gsl_matrix_complex_get(eigen_system.second.get(),j,i)) - gsl_complex_abs(gsl_matrix_complex_get(eigenvectors,j,i));
      if (eps > 1.0e-10) {
        std::cout << gsl_complex_abs(gsl_matrix_complex_get(eigen_system.second.get(),j,i)) << std::endl;
        std::cout << gsl_complex_abs(gsl_matrix_complex_get(eigenvectors,j,i)) << std::endl;
      }
    }
  }

  /*
  std::cout << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << std::endl;
  std::cout << 1000.0 * (c_end_2-c_start_2) / CLOCKS_PER_SEC << std::endl;
  std::cout << std::chrono::duration<double, std::milli>(t_end-t_start).count() << std::endl;
  std::cout << std::chrono::duration<double, std::milli>(t_end_2-t_start_2).count() << std::endl;
  */

  return 0;
}
