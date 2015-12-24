#ifndef SQUIDS_DETAIL_MATRIXEXP_H
#define SQUIDS_DETAIL_MATRIXEXP_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>

namespace squids{

namespace math_detail{

void gsl_matrix_complex_mul(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s);
void gsl_matrix_complex_add(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s);
void gsl_matrix_complex_sub(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s);

double exact_1_norm(const gsl_matrix_complex *M);
double one_normest_matrix_power(const gsl_matrix_complex * M, unsigned int p);
double one_normest_core(const gsl_matrix_complex *A, unsigned int t = 2, unsigned int itmax = 5);
int ell(const gsl_matrix_complex * A,unsigned int m);

void pade3(const gsl_matrix_complex *A, const gsl_matrix_complex *id,
            const gsl_matrix_complex *A2,
            gsl_matrix_complex *U, gsl_matrix_complex *V);
void pade5(const gsl_matrix_complex *A, const gsl_matrix_complex *id,
            const gsl_matrix_complex *A2, const gsl_matrix_complex *A4,
            gsl_matrix_complex *U, gsl_matrix_complex *V);
void pade7(const gsl_matrix_complex *A, const gsl_matrix_complex *id,
            const gsl_matrix_complex *A2, const gsl_matrix_complex *A4, const gsl_matrix_complex *A6,
            gsl_matrix_complex *U, gsl_matrix_complex *V);
void pade9(const gsl_matrix_complex *A, const gsl_matrix_complex *id,
            const gsl_matrix_complex *A2, const gsl_matrix_complex *A4, const gsl_matrix_complex *A6, const gsl_matrix_complex *A8,
            gsl_matrix_complex *U, gsl_matrix_complex *V);
void pade13(const gsl_matrix_complex *A, const gsl_matrix_complex *id,
            const gsl_matrix_complex *A2, const gsl_matrix_complex *A4, const gsl_matrix_complex *A6,
            gsl_matrix_complex *U, gsl_matrix_complex *V);

void solve_P_Q(const gsl_matrix_complex * U, const gsl_matrix_complex *V,
                gsl_matrix_complex * SPQ);

void matrix_exponential(gsl_matrix_complex * eA, const gsl_matrix_complex *A);

} // close math_detal namespace
} // close squids namespace

#endif
