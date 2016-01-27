#ifndef SQUIDS_DETAIL_MATRIXEXP_H
#define SQUIDS_DETAIL_MATRIXEXP_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>

namespace squids{
namespace math_detail{

void gsl_matrix_complex_mul(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s);
void gsl_matrix_complex_add(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s);
void gsl_matrix_complex_sub(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s);

void matrix_exponential(gsl_matrix_complex * eA, const gsl_matrix_complex *A);

} // close math_detail namespace
} // close squids namespace

#endif
