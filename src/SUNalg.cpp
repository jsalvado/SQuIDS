 /******************************************************************************
 *    This program is free software: you can redistribute it and/or modify     *
 *   it under the terms of the GNU General Public License as published by      *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   This program is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                             *
 *   Authors:                                                                  *
 *      Carlos Arguelles (University of Wisconsin Madison)                     *
 *         carguelles@icecube.wisc.edu                                         *
 *      Christopher Weaver (University of Wisconsin Madison)                   *
 *         chris.weaver@icecube.wisc.edu                                       *
 *      Jordi Salvado (University of Wisconsin Madison)                        *
 *         jsalvado@icecube.wisc.edu                                           *
 ******************************************************************************/

#include "SUNalg.h"

#include <ostream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <iostream>
#include <complex>
#define KRONECKER(i,j)  ( (i)==(j) ? 1 : 0 )
#define SQ(x)  (( x ) * ( x ))

#include <iomanip>

/*
-----------------------------------------------------------------------
Auxiliary Functions
-----------------------------------------------------------------------
*/

namespace{
gsl_complex to_gsl(std::complex<double> c){
  gsl_complex g;
  GSL_SET_COMPLEX(&g,c.real(),c.imag());
  return(g);
}

void gsl_matrix_complex_normalize( gsl_matrix_complex * m) {
  for ( size_t i = 0; i < m->size1; i++ ){
    double norm = 0;
    for ( size_t j = 0; j < m->size2; j++ ){
      norm += gsl_complex_abs2(gsl_matrix_complex_get(m,j,i));
    }
    norm = sqrt(norm);
    for ( size_t j = 0; j < m->size2; j++ ){
      gsl_matrix_complex_set(m,j,i,gsl_complex_div_real(gsl_matrix_complex_get(m,j,i),norm));
    }
  }
}
}

namespace squids{
/*
-----------------------------------------------------------------------
SU_vector implementation
-----------------------------------------------------------------------
*/


/*
-----------------------------------------------------------------------
Constructors
-----------------------------------------------------------------------
*/


SU_vector::SU_vector():
dim(0),
size(0),
components(nullptr),
isinit(false),
isinit_d(false){}

SU_vector::SU_vector(const SU_vector& V):
dim(V.dim),
size(V.size),
components(new double[size]),
isinit(true),
isinit_d(false){
  std::copy(V.components,V.components+size,components);
}

SU_vector::SU_vector(SU_vector&& V):
dim(V.dim),
size(V.size),
components(V.components),
isinit(V.isinit),
isinit_d(V.isinit_d){
  if(V.isinit){
    V.dim=0;
    V.size=0;
    V.components=nullptr;
    V.isinit=false;
  }
}

SU_vector::SU_vector(unsigned int d,double* comp):
dim(d),
size(dim*dim),
components(comp),
isinit(false),
isinit_d(true)
{
  if(dim>SQUIDS_MAX_HILBERT_DIM)
    throw std::runtime_error("SU_vector::SU_vector(unsigned int, double*): Invalid size: only up to SU(" SQUIDS_MAX_HILBERT_DIM_STR ") is supported");
};

SU_vector::SU_vector(unsigned int d):
dim(d),
size(dim*dim),
components(new double[size]),
isinit(true),
isinit_d(false)
{
  if(dim>SQUIDS_MAX_HILBERT_DIM)
    throw std::runtime_error("SU_vector::SU_vector(unsigned int): Invalid size: only up to SU(" SQUIDS_MAX_HILBERT_DIM_STR ") is supported");
  std::fill(components,components+size,0.0);
};

SU_vector::SU_vector(const std::vector<double>& comp):
dim(sqrt(comp.size())),
size(comp.size()),
components(new double[size]),
isinit(true),
isinit_d(false)
{
  if(dim*dim!=size)
    throw std::runtime_error("SU_vector::SU_vector(std::vector<double>): Vector size must be a square");
  if(dim>SQUIDS_MAX_HILBERT_DIM)
    throw std::runtime_error("SU_vector::SU_vector(std::vector<double>): Invalid size: only up to SU(" SQUIDS_MAX_HILBERT_DIM_STR ") is supported");
  std::copy(comp.begin(),comp.end(),components);
};

SU_vector::SU_vector(const gsl_matrix_complex* m):
dim(m->size1),
size(dim*dim),
components(new double[size]),
isinit(true),
isinit_d(false)
{
  if(m->size1!=m->size2)
    throw std::runtime_error("SU_vector::SU_vector(gsl_matrix_complex*): Matrix must be square");
  if(dim>SQUIDS_MAX_HILBERT_DIM)
    throw std::runtime_error("SU_vector::SU_vector(gsl_matrix_complex*): Invalid size: only up to SU(" SQUIDS_MAX_HILBERT_DIM_STR ") is supported");

  std::fill(components,components+size,0.0);

  double m_real[dim][dim]; double m_imag[dim][dim];
  for(unsigned int i=0; i<dim; i++){
    for(unsigned int j=0; j<dim; j++){
      m_real[i][j] = GSL_REAL(gsl_matrix_complex_get(m,i,j));
      m_imag[i][j] = GSL_IMAG(gsl_matrix_complex_get(m,i,j));
    }
  }
  // the following rules are valid ONLY when the initial
  // matrix is hermitian

#include "MatrixToSUSelect.txt"
};

namespace{
  struct sq_array_2D{
    unsigned int d;
    double* data;

    struct index_proxy{
      double* data;
      double operator[](unsigned int j) const{
        return(data[j]);
      }
    };
    index_proxy operator[](unsigned int i) const{
      return(index_proxy{data+d*i});
    }
  };

void ComponentsFromMatrices(double* components, unsigned int dim, const sq_array_2D& m_real, const sq_array_2D& m_imag){
#include "MatrixToSUSelect.txt"
}
} // close unnamed namespace

SU_vector SU_vector::Projector(unsigned int d, unsigned int ii){
  if(d>SQUIDS_MAX_HILBERT_DIM)
    throw std::runtime_error("SU_vector::Projector(unsigned int, unsigned int): Invalid size: only up to SU(" SQUIDS_MAX_HILBERT_DIM_STR ") is supported");
  if(ii>=d)
    throw std::runtime_error("SU_vector::Projector(unsigned int, unsigned int): Invalid component: must be smaller than dimension");

  double m_real[d][d]; double m_imag[d][d];
  for(unsigned int i=0; i<d; i++){
    for(unsigned int j=0; j<d; j++){
      m_real[i][j] = KRONECKER(i,j)*KRONECKER(i,ii);
      m_imag[i][j] = 0.0;
    }
  }

  SU_vector v(d);
  ComponentsFromMatrices(v.components,d,sq_array_2D{d,m_real[0]},sq_array_2D{d,m_imag[0]});
  return(v);
}

SU_vector SU_vector::Identity(unsigned int d){
  if(d>SQUIDS_MAX_HILBERT_DIM)
    throw std::runtime_error("SU_vector::Identity(unsigned int): Invalid size: only up to SU(" SQUIDS_MAX_HILBERT_DIM_STR ") is supported");

  double m_real[d][d]; double m_imag[d][d];
  for(unsigned int i=0; i<d; i++){
    for(unsigned int j=0; j<d; j++){
      m_real[i][j] = KRONECKER(i,j);
      m_imag[i][j] = 0.0;
    }
  }

  SU_vector v(d);
  ComponentsFromMatrices(v.components,d,sq_array_2D{d,m_real[0]},sq_array_2D{d,m_imag[0]});
  return(v);
}

SU_vector SU_vector::PosProjector(unsigned int d, unsigned int ii){
  if(d>SQUIDS_MAX_HILBERT_DIM)
    throw std::runtime_error("SU_vector::PosProjector(unsigned int, unsigned int): Invalid size: only up to SU(" SQUIDS_MAX_HILBERT_DIM_STR ") is supported");
  if(ii>=d)
    throw std::runtime_error("SU_vector::PosProjector(unsigned int, unsigned int): Invalid component: must be smaller than dimension");

  double m_real[d][d]; double m_imag[d][d];
  for(unsigned int i=0; i<d; i++){
    for(unsigned int j=0; j<d; j++){
      if(i<ii)
        m_real[i][j] = KRONECKER(i,j);
      else
        m_real[i][j] = 0.0;
      m_imag[i][j] = 0.0;
    }
  }

  SU_vector v(d);
  ComponentsFromMatrices(v.components,d,sq_array_2D{d,m_real[0]},sq_array_2D{d,m_imag[0]});
  return(v);
}

SU_vector SU_vector::NegProjector(unsigned int d, unsigned int ii){
  if(d>SQUIDS_MAX_HILBERT_DIM)
    throw std::runtime_error("SU_vector::NegProjector(unsigned int, unsigned int): Invalid size: only up to SU(" SQUIDS_MAX_HILBERT_DIM_STR ") is supported");
  if(ii>=d)
    throw std::runtime_error("SU_vector::NegProjector(unsigned int, unsigned int): Invalid component: must be smaller than dimension");

  double m_real[d][d]; double m_imag[d][d];
  for(unsigned int i=0; i<d; i++){
    for(unsigned int j=0; j<d; j++){
      if(d-i<ii)
        m_real[i][j] = KRONECKER(i,j);
      else
        m_real[i][j] = 0.0;
      m_imag[i][j] = 0.0;
    }
  }

  SU_vector v(d);
  ComponentsFromMatrices(v.components,d,sq_array_2D{d,m_real[0]},sq_array_2D{d,m_imag[0]});
  return(v);
}

SU_vector SU_vector::Generator(unsigned int d, unsigned int ii){
  if(d>SQUIDS_MAX_HILBERT_DIM)
    throw std::runtime_error("SU_vector::Component(unsigned int, unsigned int): Invalid size: only up to SU(" SQUIDS_MAX_HILBERT_DIM_STR ") is supported");
  if(ii>=d*d)
    throw std::runtime_error("SU_vector::Component(unsigned int, unsigned int): Invalid component: must be smaller than dimension");
  SU_vector v(d);
  v.components[ii] = 1.0;
  return(v);
}


/*
-----------------------------------------------------------------------
Operations
-----------------------------------------------------------------------
*/

void SU_vector::SetAllComponents(double x){
  for(unsigned int i = 0; i < size; i++)
    components[i] = x;
}

std::vector<double> SU_vector::GetComponents() const{
  std::vector<double> x ( dim*dim );
  for (unsigned int i = 0; i < dim*dim ; i++)
    x[i] = components[i];
  return x;
}

std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)>
SU_vector::GetGSLMatrix() const {
  if( !isinit and !isinit_d)
    throw std::runtime_error("SU_vector::GetGSLMatrix(): SU_vector not initialized.");
  gsl_matrix_complex * matrix = gsl_matrix_complex_alloc(dim,dim);
#include "SUToMatrixSelect.txt"
  return std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)>(matrix,gsl_matrix_complex_free);
}

void gsl_complex_matrix_exponential(gsl_matrix_complex *eA, const gsl_matrix_complex *A, unsigned int dimx){
  int j,k=0;
  gsl_complex temp;
  gsl_matrix* matreal =gsl_matrix_alloc(2*dimx,2*dimx);
  gsl_matrix* expmatreal =gsl_matrix_alloc(2*dimx,2*dimx);
  //Converting the complex matrix into real one using A=[Areal, Aimag;-Aimag,Areal]
  for (j = 0; j < dimx;j++)
    for (k = 0; k < dimx;k++){
      temp=gsl_matrix_complex_get(A,j,k);
      gsl_matrix_set(matreal,j,k,-GSL_IMAG(temp));
      gsl_matrix_set(matreal,dimx+j,dimx+k,-GSL_IMAG(temp));
      gsl_matrix_set(matreal,j,dimx+k,GSL_REAL(temp));
      gsl_matrix_set(matreal,dimx+j,k,-GSL_REAL(temp));
    }

  gsl_linalg_exponential_ss(matreal,expmatreal,GSL_PREC_DOUBLE);

  double realp;
  double imagp;
  for (j = 0; j < dimx;j++)
    for (k = 0; k < dimx;k++){
      realp=gsl_matrix_get(expmatreal,j,k);
      imagp=gsl_matrix_get(expmatreal,j,dimx+k);
      gsl_matrix_complex_set(eA,j,k,gsl_complex_rect(realp,imagp));
    }
  gsl_matrix_free(matreal);
  gsl_matrix_free(expmatreal);
}

void gsl_matrix_complex_change_basis_UCMU(const gsl_matrix_complex* U, gsl_matrix_complex* M){
  int numneu = U->size1;
  gsl_matrix_complex* U1 = gsl_matrix_complex_alloc(numneu,numneu);
  gsl_matrix_complex* U2 = gsl_matrix_complex_alloc(numneu,numneu);
  gsl_matrix_complex_memcpy(U1,U);
  gsl_matrix_complex_memcpy(U2,U);

  gsl_matrix_complex* T1 = gsl_matrix_complex_alloc(numneu,numneu);

  // doing : U M U^dagger
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,
                 gsl_complex_rect(1.0,0.0),M,
                 U1,gsl_complex_rect(0.0,0.0),T1);
  gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,
                 gsl_complex_rect(1.0,0.0),U2,
                 T1,gsl_complex_rect(0.0,0.0),M);

  gsl_matrix_complex_free(U1);
  gsl_matrix_complex_free(U2);
  gsl_matrix_complex_free(T1);
}

SU_vector SU_vector::UTransform(const SU_vector& v, gsl_complex scale) const{
  auto mv=v.GetGSLMatrix();
  auto mu=(*this).GetGSLMatrix();
  //rescale the matrix
  gsl_matrix_complex_scale(mv.get(),scale);

  //declare and calculate matrix exponential
  gsl_matrix_complex* em = gsl_matrix_complex_alloc (dim, dim);
  math_detail::matrix_exponential(em,mv.get());
  // sandwich
  gsl_matrix_complex_change_basis_UCMU(em, mu.get());
  gsl_matrix_complex_free(em);

  return SU_vector(mu.get());
}

std::pair<std::unique_ptr<gsl_vector,void (*)(gsl_vector*)>,
std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)>>
SU_vector::GetEigenSystem(bool order) const{
  gsl_vector * eigenvalues = gsl_vector_alloc(dim);
  gsl_matrix_complex * eigenvectors = gsl_matrix_complex_alloc(dim,dim);
  switch (dim) {
    case 3:
          {
#include "EigenSystemSU3.txt"
          }
          break;
    default:
      auto matrix=(*this).GetGSLMatrix();
      gsl_eigen_hermv_workspace * ws = gsl_eigen_hermv_alloc(dim);
      gsl_eigen_hermv(matrix.get(),eigenvalues,eigenvectors,ws);
      gsl_eigen_hermv_free(ws);
  }
  // sorting eigenvalues
  if (order)
    gsl_eigen_hermv_sort(eigenvalues,eigenvectors,GSL_EIGEN_SORT_VAL_ASC);
  return std::make_pair(
    std::unique_ptr<gsl_vector,void (*)(gsl_vector*)>(eigenvalues,gsl_vector_free),
    std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)>(eigenvectors,gsl_matrix_complex_free));
}

/*
SU_vector SU_vector::Rotate(const gsl_matrix_complex* m) const{
  const SU_vector& suv1=*this;
  SU_vector suv_new(dim);
  if ( m->size1 != dim or m->size2 != dim )
    throw std::runtime_error("SU_vector::Rotate(gsl_matrix_complex): matrix dimensions and SU_vector dimensions do not match.");
#include "unitary_rotation_switcher.h"
  return suv_new;
}
*/

SU_vector SU_vector::Rotate(const gsl_matrix_complex* U) const{
  if ( U->size1 != dim or U->size2 != dim )
    throw std::runtime_error("SU_vector::Rotate(gsl_matrix_complex): matrix dimensions and SU_vector dimensions do not match.");

  const SU_vector& suv1=*this;
  auto m = suv1.GetGSLMatrix();
  gsl_matrix_complex_change_basis_UCMU(U,m.get());
  return SU_vector(std::move(m));
}

SU_vector SU_vector::Rotate(unsigned int ii, unsigned int jj, double th, double del) const{
  const SU_vector& suv=*this;
  SU_vector suv_rot(dim);
  unsigned int i=ii+1, j=jj+1; //convert to 1 based indices to interface with Mathematica generated code

  assert(i<j && "Components selected for rotation must be in ascending order");

#include "rotation_switcher.h"

  return suv_rot;
}

void SU_vector::RotateToB0(const Const& param){
  for(unsigned int j=1; j<dim; j++){
    for(unsigned int i=0; i<j; i++)
      *this=Rotate(i,j,-param.GetMixingAngle(i,j),param.GetPhase(i,j)); //note negated angle
  }
}

void SU_vector::RotateToB1(const Const& param){
  for(unsigned int j=dim-1; j>0; j--){
    for(unsigned int i=j-1; i<dim; i--)
      *this=Rotate(i,j,param.GetMixingAngle(i,j),param.GetPhase(i,j));
  }
}

double SUTrace(const SU_vector& suv1,const SU_vector& suv2){
  double gen_trace = 0.0;
  double id_trace = 0.0;

  for(unsigned int i=1; i < suv1.Size(); i++)
    gen_trace += (suv1.components[i])*(suv2.components[i]);

  id_trace = (suv1.components[0])*(suv2.components[0])*double(suv1.dim);
  return id_trace+2.0*gen_trace;
}

bool SU_vector::operator==(const SU_vector& other) const{
  if(dim != other.dim) //vectors of different sizes cannot be equal
    return false;
  //if one vector contains data and the other does not they cannot be equal
  if((isinit || isinit_d) != (other.isinit || other.isinit_d))
    return false;
  //if both vectors are empty, consider them equal
  if(!isinit && !isinit_d)
    return true;
  //test whether all components are equal
  for(unsigned int i=0; i < size; i++){
    if(components[i] != other.components[i])
      return false;
  }
  return true;
}

double SU_vector::operator*(const SU_vector &other) const{
  if(size!=other.size)
    throw std::runtime_error("Non-matching dimensions in SU_vector inner product");
  return SUTrace(*this,other);
}

SU_vector& SU_vector::operator=(const SU_vector& other){
  if(this==&other)
    return(*this);
  if(size!=other.size){
    if(isinit_d) //can't resize
      throw std::runtime_error("Non-matching dimensions in assignment to SU_vector with external storage");
    //can resize
    if(isinit)
      delete[] components;
    dim=other.dim;
    size=other.size;
    components=new double[size];
    isinit=true;
  }

  std::copy(other.components,other.components+size,components);
  return *this;
}

SU_vector& SU_vector::operator=(SU_vector&& other){
  if(this==&other)
    return(*this);
  if(isinit){ //this vector owns storage
    std::swap(dim,other.dim);
    std::swap(size,other.size);
    std::swap(components,other.components);
    std::swap(isinit,other.isinit);
    std::swap(isinit_d,other.isinit_d);
  }
  else if(isinit_d){ //this vector has non-owned storage
    //conservatively assume that we cannot shift our storage, so fall back on an actual copy
    if(size!=other.size)
      throw std::runtime_error("Non-matching dimensions in assignment to SU_vector with external storage");
    std::copy(other.components,other.components+size,components);
  }
  else{
    //this vector has no storage, so just take everything from other
    dim = other.dim;
    size = other.size;
    components = other.components;
    isinit = other.isinit;
    isinit_d = other.isinit_d;
    if(other.isinit){ //other must relinquish ownership
      other.dim=0;
      other.size=0;
      other.components=nullptr;
      other.isinit=false;
    }
  }
  return(*this);
}

SU_vector& SU_vector::operator+=(const SU_vector &other){
  if(size!=other.size)
    throw std::runtime_error("Non-matching dimensions in SU_vector increment");
  for(unsigned int i=0; i < size; i++)
    components[i] += other.components[i];
  return *this;
}

SU_vector& SU_vector::operator-=(const SU_vector &other){
  if(size!=other.size)
    throw std::runtime_error("Non-matching dimensions in SU_vector decrement");
  for(unsigned int i=0; i < size; i++)
    components[i] -= other.components[i];
  return *this;
}

SU_vector& SU_vector::operator *=(double x){
  for(unsigned int i = 0; i < size; i++)
    components[i] *= x;
  return *this;
}

SU_vector& SU_vector::operator /=(double x){
  for(unsigned int i = 0; i < size; i++)
    components[i] /= x;
  return *this;
}

std::ostream& operator<<(std::ostream& os, const SU_vector& V){
  for(unsigned int i=0; i< V.size-1; i++)
    os << V.components[i] << "  ";
  os << V.components[V.size-1];
  return os;
}
	
detail::iCommutatorProxy iCommutator(const SU_vector& suv1,const SU_vector& suv2){
	if(suv1.dim!=suv2.dim)
		throw std::runtime_error("Commutator error, non-matching dimensions ");
	return(detail::iCommutatorProxy{suv1,suv2});
}

detail::ACommutatorProxy ACommutator(const SU_vector& suv1,const SU_vector& suv2){
	if(suv1.dim!=suv2.dim)
		throw std::runtime_error("Anti Commutator error: non-matching dimensions ");
	return(detail::ACommutatorProxy{suv1,suv2});
}

} //namespace squids
