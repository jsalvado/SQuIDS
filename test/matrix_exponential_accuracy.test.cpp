// exp(-iK) rho exp(iK) from SU_vector::UTransform and detail::HermitianEvolver
// against a reference built from many small steps, for |K| spanning every
// Pade degree of the scaling-and-squaring exponential (the degree 9 branch,
// 0.95 < |K| <= 2.1, used to drop the A^8 terms). Prints nothing on success.
#include <cmath>
#include <iostream>
#include <SQuIDS/SUNalg.h>
#include <SQuIDS/detail/hermitian_evolver.h>

using squids::SU_vector;

int main(){
  for(unsigned int dim = 2; dim <= 6; dim++){
    squids::detail::HermitianEvolver ev(dim);
    SU_vector K(dim), rho(dim), ref(dim), u(dim), e(dim);
    for(unsigned int k = 0; k < dim*dim; k++){
      K[k] = std::cos(1.3*k + 0.2*dim);
      rho[k] = std::sin(0.7*k + 1.1);
    }
    double knorm = 0;
    for(unsigned int k = 0; k < dim*dim; k++) knorm = std::max(knorm, std::abs(K[k]));
    for(double scale = 0.01; scale < 300; scale *= 1.3){
      const double s = scale/knorm;
      // reference: 2^12 steps, each with |K|/4096 < 0.1 (degree 3 or 5 branch)
      ref = rho;
      for(int i = 0; i < 4096; i++)
        ref = ref.UTransform(K, gsl_complex_rect(0., s/4096.));
      u = rho.UTransform(K, gsl_complex_rect(0., s));
      SU_vector Ks = K*s;
      ev.Evolve(Ks, rho, e);
      double du = 0, de = 0, nr = 0;
      for(unsigned int k = 0; k < dim*dim; k++){
        du = std::max(du, std::abs(u[k] - ref[k]));
        de = std::max(de, std::abs(e[k] - ref[k]));
        nr = std::max(nr, std::abs(ref[k]));
      }
      if(du > 1e-10*nr)
        std::cout << "UTransform SU(" << dim << ") |K|~" << scale << ": rel error " << du/nr << std::endl;
      if(de > 1e-10*nr)
        std::cout << "HermitianEvolver SU(" << dim << ") |K|~" << scale << ": rel error " << de/nr << std::endl;
    }
  }
  return 0;
}
