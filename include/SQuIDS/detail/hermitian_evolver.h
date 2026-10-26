#ifndef SQUIDS_DETAIL_HERMITIAN_EVOLVER_H
#define SQUIDS_DETAIL_HERMITIAN_EVOLVER_H

// exp(-iK) rho exp(iK) for Hermitian K, through an eigendecomposition of K.
//
// SU_vector::UTransform evaluates the same thing with a general complex matrix
// exponential (scaling and squaring), whose cost grows with the norm of K and
// which allocates GSL matrices. Exponential integrators call this once or twice
// per step, so this version keeps every buffer and the eigensolver workspace,
// and its cost does not depend on |K|. The matrix conversions and the
// eigensolver are public for integrators that work in matrix form.

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include "../SUNalg.h"

namespace squids{
namespace detail{

class HermitianEvolver {
public:
  typedef std::complex<double> cplx;

  explicit HermitianEvolver(unsigned int dim):
  N(dim), n2(dim*dim), basis(n2*n2), gram_inv(n2*n2, 0.0), proj(n2, 0.0),
  Km(n2), Rm(n2), Tm(n2), Wm(n2), lam(dim),
  gK(gsl_matrix_complex_alloc(dim,dim)), gW(gsl_matrix_complex_alloc(dim,dim)),
  geval(gsl_vector_alloc(dim)), ws(gsl_eigen_hermv_alloc(dim))
  {
    // Read the SU_vector basis off SQuIDS itself: B_k is the matrix of the
    // vector whose only nonzero component is k.
    SU_vector e(N);
    for(unsigned int k=0; k<n2; k++){
      e.SetAllComponents(0.0);
      e[k] = 1.0;
      e.GetGSLMatrix(gK);
      for(unsigned int i=0; i<N; i++)
        for(unsigned int j=0; j<N; j++){
          gsl_complex z = gsl_matrix_complex_get(gK,i,j);
          basis[k*n2+i*N+j] = cplx(GSL_REAL(z), GSL_IMAG(z));
        }
    }
    // components from a Hermitian matrix M: c = G^-1 b, b_j = Tr(B_j M),
    // G_jk = Tr(B_j B_k) (real because every B is Hermitian)
    gsl_matrix* G = gsl_matrix_alloc(n2,n2);
    for(unsigned int a=0; a<n2; a++)
      for(unsigned int b=0; b<n2; b++)
        gsl_matrix_set(G,a,b,TraceProduct(&basis[a*n2],&basis[b*n2]));
    gsl_permutation* p = gsl_permutation_alloc(n2);
    gsl_matrix* Ginv = gsl_matrix_alloc(n2,n2);
    int sign;
    gsl_linalg_LU_decomp(G,p,&sign);
    gsl_linalg_LU_invert(G,p,Ginv);
    for(unsigned int a=0; a<n2; a++)
      for(unsigned int b=0; b<n2; b++)
        gram_inv[a*n2+b] = gsl_matrix_get(Ginv,a,b);
    gsl_matrix_free(Ginv);
    gsl_permutation_free(p);
    gsl_matrix_free(G);
  }

  ~HermitianEvolver(){
    gsl_eigen_hermv_free(ws);
    gsl_vector_free(geval);
    gsl_matrix_complex_free(gW);
    gsl_matrix_complex_free(gK);
  }
  HermitianEvolver(const HermitianEvolver&) = delete;
  HermitianEvolver& operator=(const HermitianEvolver&) = delete;

  unsigned int Dimension() const { return N; }

  ///\brief out = exp(-i K) rho exp(i K); out may alias rho.
  void Evolve(const SU_vector& K, const SU_vector& rho, SU_vector& out){
    ToMatrix(K, &Km[0]);
    ToMatrix(rho, &Rm[0]);
    Eigen(&Km[0], &lam[0], &Wm[0]);
    Adjoint(&Wm[0], &Rm[0], &Tm[0]);   // T = W^dagger R W
    PhaseRotate(&lam[0], 1.0, &Tm[0]); // T_ij *= exp(-i (l_i - l_j))
    Conjugate(&Wm[0], &Tm[0], &Rm[0]); // R = W T W^dagger
    FromMatrix(&Rm[0], out);
  }

  ///\brief M = matrix of v (row major, N x N).
  void ToMatrix(const SU_vector& v, cplx* M) const{
    for(unsigned int i=0; i<n2; i++)
      M[i] = 0;
    for(unsigned int k=0; k<n2; k++){
      const double c = v[k];
      if(c == 0) continue;
      const cplx* B = &basis[k*n2];
      for(unsigned int i=0; i<n2; i++)
        M[i] += c*B[i];
    }
  }

  ///\brief v = SU_vector of the Hermitian matrix M.
  void FromMatrix(const cplx* M, SU_vector& v){
    for(unsigned int a=0; a<n2; a++)
      proj[a] = TraceProduct(&basis[a*n2], M);
    for(unsigned int a=0; a<n2; a++){
      double c = 0;
      for(unsigned int b=0; b<n2; b++)
        c += gram_inv[a*n2+b]*proj[b];
      v[a] = c;
    }
  }

  ///\brief Hermitian M = W diag(lam) W^dagger; M is not modified.
  void Eigen(const cplx* M, double* l, cplx* W){
    for(unsigned int i=0; i<N; i++)
      for(unsigned int j=0; j<N; j++)
        gsl_matrix_complex_set(gK,i,j,gsl_complex_rect(M[i*N+j].real(), M[i*N+j].imag()));
    gsl_eigen_hermv(gK, geval, gW, ws);
    for(unsigned int i=0; i<N; i++){
      l[i] = gsl_vector_get(geval,i);
      for(unsigned int j=0; j<N; j++){
        gsl_complex z = gsl_matrix_complex_get(gW,i,j);
        W[i*N+j] = cplx(GSL_REAL(z), GSL_IMAG(z));
      }
    }
  }

  ///\brief T = W^dagger R W.
  void Adjoint(const cplx* W, const cplx* R, cplx* T){
    for(unsigned int i=0; i<N; i++)
      for(unsigned int j=0; j<N; j++){
        cplx s = 0;
        for(unsigned int k=0; k<N; k++)
          s += std::conj(W[k*N+i])*R[k*N+j];
        tmp[i*N+j] = s;
      }
    Product(&tmp[0], W, T);
  }

  ///\brief R = W T W^dagger.
  void Conjugate(const cplx* W, const cplx* T, cplx* R){
    Product(W, T, &tmp[0]);
    for(unsigned int i=0; i<N; i++)
      for(unsigned int j=0; j<N; j++){
        cplx s = 0;
        for(unsigned int k=0; k<N; k++)
          s += tmp[i*N+k]*std::conj(W[j*N+k]);
        R[i*N+j] = s;
      }
  }

  ///\brief T_ij *= exp(-i t (l_i - l_j)): conjugation by exp(-i t diag(l)).
  void PhaseRotate(const double* l, double t, cplx* T) const{
    for(unsigned int i=0; i<N; i++)
      for(unsigned int j=0; j<N; j++)
        if(i != j)
          T[i*N+j] *= std::polar(1.0, -t*(l[i]-l[j]));
  }

private:
  unsigned int N, n2;
  std::vector<cplx> basis;      // n2 matrices of N x N
  std::vector<double> gram_inv; // n2 x n2
  std::vector<double> proj;
  std::vector<cplx> Km, Rm, Tm, Wm;
  std::vector<cplx> tmp = std::vector<cplx>(36);
  std::vector<double> lam;
  gsl_matrix_complex *gK, *gW;
  gsl_vector* geval;
  gsl_eigen_hermv_workspace* ws;

  // Re Tr(X Y) for N x N complex matrices
  double TraceProduct(const cplx* X, const cplx* Y) const{
    double s = 0;
    for(unsigned int i=0; i<N; i++)
      for(unsigned int j=0; j<N; j++)
        s += (X[i*N+j]*Y[j*N+i]).real();
    return s;
  }

  void Product(const cplx* A, const cplx* B, cplx* C) const{
    for(unsigned int i=0; i<N; i++)
      for(unsigned int j=0; j<N; j++){
        cplx s = 0;
        for(unsigned int k=0; k<N; k++)
          s += A[i*N+k]*B[k*N+j];
        C[i*N+j] = s;
      }
  }
};

} // namespace detail
} // namespace squids

#endif
