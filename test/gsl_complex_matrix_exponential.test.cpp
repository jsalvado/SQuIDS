#include <iostream>
#include <iomanip>
#include <fstream>
#include <SQuIDS/SUNalg.h>
#include <gsl/gsl_blas.h>
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

void fill_matrix(gsl_matrix_complex* m, 
                 const std::initializer_list<std::initializer_list<std::initializer_list<double>>>& data){
  unsigned int size2=data.size();
  unsigned int size1=0;
  if(size2){
    size1=data.begin()->size();
    for(const auto& row : data)
      assert(row.size()==size1);
  }
  if(m->size1!=size1 || m->size2!=size2){
    gsl_matrix_complex_free(m);
    m=gsl_matrix_complex_alloc(size1,size2);
  }
  unsigned int i=0;
  for(const auto& row : data){
    unsigned int j=0;
    for(const auto& entry : row){
      auto it=entry.begin();
      double re=*it++;
      double im=*it;
      gsl_matrix_complex_set(m,i,j++,gsl_complex_rect(re, im));
    }
    i++;
  }
}

///given M and e=exp(M) from a trusted calculation, check that eM==e to within tol
bool check_exponential(gsl_matrix_complex* M, gsl_matrix_complex* e, gsl_matrix_complex* eM, double tol){
  double maxErr=0;
  auto close_enough=[tol,&maxErr](double observed, double expected){
    double err;
    if(!expected)
      err=std::abs(observed);
    else
      err=std::abs(observed-expected)/expected;
    if(err>maxErr)
      maxErr=err;
    return err<tol;
  };
  unsigned int dim=M->size1; //assume M is square
  bool good=true;
  for(unsigned int i=0; i<dim; i++){
    for(unsigned int j=0; j<dim; j++){
      auto cc = gsl_matrix_complex_get(eM,i,j);
      auto ce = gsl_matrix_complex_get(e,i,j);
      if(!close_enough(GSL_REAL(cc),GSL_REAL(ce)) || !close_enough(GSL_IMAG(cc),GSL_IMAG(ce)))
        good=false;
    }
  }
  if(!good){
    std::cout << "Error:\nInputMatrix:\n";
    simple_print(M);
    std::cout << "Expected Exponential:\n";
    simple_print(e);
    std::cout << "Computed Exponential:\n";
    simple_print(eM);
    std::cout << "Max. Err.: " << maxErr << '\n';
  }
  return good;
}

int main(){
  unsigned int dim = 3;
  gsl_matrix_complex * M = gsl_matrix_complex_alloc(dim,dim);
  gsl_matrix_complex * eM = gsl_matrix_complex_alloc(dim,dim);

  // set printing precision
  std::cout << std::scientific;
  std::cout << std::setprecision(8);
  // first matrix
  fill_matrix(M,{
    {{1,2},{3,-2},{5,-7}},
    {{3,9},{4,-2},{9,-7}},
    {{5,7},{1,-1},{5,-3}},
  });
  // calculating matrix exponential
  squids::math_detail::matrix_exponential(eM,M);
  // simple print matrix
  simple_print(eM);
  std::cout << std::endl;

  // second matrix
  fill_matrix(M,{
    {{ 0,0.1},{  1,  20},{ 50, -1}},
    {{ 4,  5},{  0,   0},{0.9,0.8}},
    {{50, 70},{0.1,-0.1},{  5, -5}},
  });
  squids::math_detail::matrix_exponential(eM,M);
  simple_print(eM);
  std::cout << std::endl;

  // third matrix
  fill_matrix(M,{
    {{0.1,0.01},{0.02,0.2},{0.3,0.3}},
    {{ -1,  -1},{  -5,  1},{0.1,  0}},
    {{  0,   0},{ 0.7,0.3},{  1,  1}},
  });
  squids::math_detail::matrix_exponential(eM,M);
  simple_print(eM);
  std::cout << std::endl;
  
  //Use the identity that, for an invertable matrix S, 
  //exp(S*M*S^-1) = S*exp(M)*S^-1
  //to test exponentiating a variety of (real) matrices whose exponentials
  //we can simply calculate directly, having choosen M to be a diagonal matrix.
  //Adjust the scale of the elements on the diagonal of M to test all of the
  //polynomial approximations used internally in the exponentiation implementation.
  auto S = gsl_matrix_complex_alloc(dim,dim);
  auto Si = gsl_matrix_complex_alloc(dim,dim);
  auto p = gsl_matrix_complex_alloc(dim,dim); //scratch space for products
  auto e = gsl_matrix_complex_alloc(dim,dim); //expected exponential
  
  //    [ 1 4 7 ]
  //S = [ 2 5 8 ]
  //    [ 3 6 1 ]
  fill_matrix(S,{
    {{1,0},{4,0},{7,0}},
    {{2,0},{5,0},{8,0}},
    {{3,0},{6,0},{1,0}},
  });
  
  //              [ -43  38 -3 ]
  //S^-1 = (1/24) [  22 -20  6 ]
  //              [  -3   6 -3 ]
  fill_matrix(Si,{
    {{-43,0},{ 38,0},{-3,0}},
    {{ 22,0},{-20,0},{ 6,0}},
    {{ -3,0},{  6,0},{-3,0}},
  });
  gsl_matrix_complex_scale(Si, gsl_complex_rect(1./24, 0));
  
  double tol=1e-12;
  
  for(double scale=.001; scale<2; scale*=sqrt(2)){
    //construct the matrix to exponentiate
    gsl_matrix_complex_set_all(M,GSL_COMPLEX_ZERO);
    gsl_matrix_complex_set(M,0,0,gsl_complex_rect(  -scale,0));
    gsl_matrix_complex_set(M,1,1,gsl_complex_rect(-2*scale,0));
    gsl_matrix_complex_set(M,2,2,gsl_complex_rect(-4*scale,0));
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,S,M,GSL_COMPLEX_ZERO,p);
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,p,Si,GSL_COMPLEX_ZERO,M);
    
    squids::math_detail::matrix_exponential(eM,M);
    
    //pseudo-analytically contruct the correct exponential
    gsl_matrix_complex_set_all(e,GSL_COMPLEX_ZERO);
    gsl_matrix_complex_set(e,0,0,gsl_complex_rect(exp(  -scale),0));
    gsl_matrix_complex_set(e,1,1,gsl_complex_rect(exp(-2*scale),0));
    gsl_matrix_complex_set(e,2,2,gsl_complex_rect(exp(-4*scale),0));
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,S,e,GSL_COMPLEX_ZERO,p);
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,p,Si,GSL_COMPLEX_ZERO,e);
    
    check_exponential(M,e,eM,tol);
  }
  gsl_matrix_complex_free(M);
  gsl_matrix_complex_free(eM);
  gsl_matrix_complex_free(S);
  gsl_matrix_complex_free(Si);
  gsl_matrix_complex_free(p);
  gsl_matrix_complex_free(e);

  return 0;
}
