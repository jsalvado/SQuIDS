#ifndef SQUIDS_DETAIL_MATRIXEXP_H
#define SQUIDS_DETAIL_MATRIXEXP_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>

namespace squids{
namespace math_detail{

struct gsl_matrix_complex_holder{
  gsl_matrix_complex* m;
  gsl_matrix_complex_holder():m(nullptr){}
  ~gsl_matrix_complex_holder(){ gsl_matrix_complex_free(m); }
  void reset(unsigned int s1, unsigned int s2){
    if(!m || s1!=m->size1 || s2!=m->size2){
      if(m)
        gsl_matrix_complex_free(m);
      m = gsl_matrix_complex_alloc(s1,s2);
    }
  }
  operator gsl_matrix_complex*(){ return m; }
  gsl_matrix_complex* operator->(){ return m; }
};

void gsl_matrix_complex_mul(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s);
void gsl_matrix_complex_add(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s);
void gsl_matrix_complex_sub(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s);

///Compute the matrix exponential
///\param eA matrix into which the result will be written
///\param A the matrix to be exponentiated
void matrix_exponential(gsl_matrix_complex * eA, const gsl_matrix_complex *A);

} // close math_detail namespace
} // close squids namespace

#endif
