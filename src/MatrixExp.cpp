#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_blas.h>

#include "SQuIDS/detail/MatrixExp.h"
#include "SQuIDS/detail/ProxyFwd.h"

namespace squids{

namespace math_detail{

void gsl_matrix_complex_print(const gsl_matrix_complex* matrix){
    for(unsigned int i = 0; i < matrix->size1; i++){
        for(unsigned int j = 0; j < matrix->size2; j++){
            std::cout << gsl_matrix_complex_get(matrix,i,j).dat[0] <<
            "+i" << gsl_matrix_complex_get(matrix,i,j).dat[1] << " ";
        }
        std::cout << std::endl;
    }
}

void gsl_vector_complex_print(const gsl_vector_complex* vec){
    for(unsigned int i = 0; i < vec->size; i++){
        std::cout << gsl_vector_complex_get(vec,i).dat[0] <<
        "+i" << gsl_vector_complex_get(vec,i).dat[1] << " ";
    }
    std::cout << std::endl;
}

void gsl_vector_print(const gsl_vector* vec){
    for(unsigned int i = 0; i < vec->size; i++){
        std::cout << gsl_vector_get(vec,i) << ' ';
    }
    std::cout << std::endl;
}

void gsl_matrix_complex_mul(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s){
  assert(O->size1 == I->size1 and O->size2 == I->size2);
  for(unsigned int i = 0; i < I->size1; i++){
    for(unsigned int j = 0; j < I->size2; j++){
      gsl_matrix_complex_set(O,i,j,gsl_complex_mul(gsl_matrix_complex_get(I,i,j),s));
    }
  }
}

void gsl_matrix_complex_add(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s){
  assert(O->size1 == I->size1 and O->size2 == I->size2);
  for(unsigned int i = 0; i < I->size1; i++){
    for(unsigned int j = 0; j < I->size2; j++){
      gsl_matrix_complex_set(O,i,j,gsl_complex_add(gsl_matrix_complex_get(O,i,j),
                                                   gsl_complex_mul(gsl_matrix_complex_get(I,i,j),s)));
    }
  }
}

void gsl_matrix_complex_sub(gsl_matrix_complex *O,const gsl_matrix_complex *I,gsl_complex s){
  assert(O->size1 == I->size1 and O->size2 == I->size2);
  for(unsigned int i = 0; i < I->size1; i++){
    for(unsigned int j = 0; j < I->size2; j++){
      gsl_matrix_complex_set(O,i,j,gsl_complex_sub(gsl_matrix_complex_get(O,i,j),
                                                   gsl_complex_mul(gsl_matrix_complex_get(I,i,j),s)));
    }
  }
}

double exact_1_norm(const gsl_matrix_complex *M){
  double norm = -1;
  for(unsigned int j = 0; j < M->size2; j++){
    double tmp = 0;
    for(unsigned int i = 0; i < M->size1; i++){
      tmp+= gsl_complex_abs(gsl_matrix_complex_get(M,i,j));
    }
    if (tmp > norm)
      norm = tmp;
  }
  return norm;
}
	
double one_normest_core(const gsl_matrix_complex *A, unsigned int t = 2, unsigned int itmax = 5);

double one_normest_matrix_power(const gsl_matrix_complex * M, unsigned int p){
  assert(M->size1 == M->size2);
  if (p == 1)
    return one_normest_core(M);

  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder MP;
  MP.reset(M->size1,M->size2);
  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder tmp;
  tmp.reset(M->size1,M->size2);
  
  for(unsigned int i=0; i<p-1; i++){
      if(i==0){
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,M,M,GSL_COMPLEX_ZERO,MP);
        //gsl_matrix_complex_print(MP);
      } else if (i%2 != 0){
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,M,MP,GSL_COMPLEX_ZERO,tmp);
        //gsl_matrix_complex_print(tmp);
      }else{
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,M,tmp,GSL_COMPLEX_ZERO,MP);
        //gsl_matrix_complex_print(MP);
      }
  }

  if(p%2!=0 and p>2)
    gsl_matrix_complex_memcpy(MP,tmp);

  double onc = one_normest_core(MP);
  return onc;
}

void resample_column(unsigned int i, gsl_matrix *X,const gsl_rng* rng){
  for(unsigned int j=0; j < X->size2; j++)
    gsl_matrix_set(X,j,i,gsl_rng_uniform_int(rng,2)*2.-1.);
}

void resample_column(unsigned int i, gsl_matrix_complex *X,const gsl_rng* rng){
  for(unsigned int j=0; j < X->size2; j++){
    //std::cout << "random guy " << gsl_rng_uniform_int(rng,2) << std::endl;
    gsl_matrix_complex_set(X,j,i,gsl_complex_rect(gsl_rng_uniform_int(rng,2)*2.-1.,0.));
  }
}

bool are_parallel(const gsl_vector *a,const gsl_vector *b){
  assert(a->size == b->size);
  double dp;
  gsl_blas_ddot(a,b,&dp);
  return (std::abs(dp-a->size) < 1.0e-6);
}

bool are_parallel(const gsl_vector_complex *a,const gsl_vector_complex *b){
  assert(a->size == b->size);
  gsl_complex dp;
  /*
  std::cout << "va" << std::endl;
  gsl_vector_complex_print(a);
  std::cout << "vb" << std::endl;
  gsl_vector_complex_print(b);
  */

  gsl_blas_zdotc(a,b,&dp);
  //std::cout << "dp " << gsl_complex_abs(dp) << " " << a->size << std::endl;

  return (std::abs(gsl_complex_abs(dp)-a->size) < 1.0e-6);
}

bool every_col_is_parallel(const gsl_matrix *X, const gsl_matrix *Y){
  for(unsigned int i=0; i < X->size2; i++){
    bool is_there_a_parallel_col = false;
    gsl_vector_const_view vx = gsl_matrix_const_column(X,i);
    for(unsigned int j=0; j < Y->size2; j++){
      gsl_vector_const_view vy = gsl_matrix_const_column(Y,j);
      if(are_parallel(&vx.vector,&vy.vector)){
        is_there_a_parallel_col = true;
        break;
      }
    }
    if(not is_there_a_parallel_col)
      return false;
  }
  return true;
}

bool every_col_is_parallel(const gsl_matrix_complex *X, const gsl_matrix_complex *Y){
  for(unsigned int i=0; i < X->size2; i++){
    bool is_there_a_parallel_col = false;
    gsl_vector_complex_const_view vx = gsl_matrix_complex_const_column(X,i);
    for(unsigned int j=0; j < Y->size2; j++){
      gsl_vector_complex_const_view vy = gsl_matrix_complex_const_column(Y,j);
      if(are_parallel(&vx.vector,&vy.vector)){
        is_there_a_parallel_col = true;
        break;
      }
    }
    if(not is_there_a_parallel_col)
      return false;
  }
  return true;
}

bool column_needs_resampling(unsigned int i, const gsl_matrix *X, const gsl_matrix *Y = nullptr){
  gsl_vector_const_view v = gsl_matrix_const_column(X,i);
  for(unsigned int j = 0; j < i; j++){
    gsl_vector_const_view cXj = gsl_matrix_const_column(X,j);
    if(are_parallel(&v.vector,&cXj.vector))
        return true;
  }

  if(Y!=nullptr){
    for(unsigned int j=0; j<Y->size2; j++){
      gsl_vector_const_view cYj = gsl_matrix_const_column(Y,j);
      if(are_parallel(&v.vector,&cYj.vector))
          return true;
    }
  }

  return false;
}

bool column_needs_resampling(unsigned int i, const gsl_matrix_complex *X, const gsl_matrix_complex *Y = nullptr){
  gsl_vector_complex_const_view v = gsl_matrix_complex_const_column(X,i);
  for(unsigned int j = 0; j < i; j++){
    gsl_vector_complex_const_view cXj = gsl_matrix_complex_const_column(X,j);
    if(are_parallel(&v.vector,&cXj.vector))
        return true;
  }

  if(Y!=nullptr){
    for(unsigned int j=0; j<Y->size2; j++){
      gsl_vector_complex_const_view cYj = gsl_matrix_complex_const_column(Y,j);
      if(are_parallel(&v.vector,&cYj.vector))
          return true;
    }
  }

  return false;
}

double norm_1d_1(const gsl_vector* v){
  double norm =0;
  for(unsigned int i=0; i<v->size; i++){
    norm += std::abs(gsl_vector_get(v,i));
  }
  return norm;
}

double norm_1d_1(const gsl_vector_complex* v){
  double norm =0;
  for(unsigned int i=0; i<v->size; i++){
    norm += gsl_complex_abs(gsl_vector_complex_get(v,i));
  }
  return norm;
}

double norm_1d_inf(const gsl_vector* v){
  double norm =-1;
  for(unsigned int i=0; i<v->size; i++){
    double tmp = std::abs(gsl_vector_get(v,i));
    if(tmp > norm)
      norm = tmp;
  }
  return norm;
}

double norm_1d_inf(const gsl_vector_complex* v){
  double norm =-1;
  for(unsigned int i=0; i<v->size; i++){
    double tmp = gsl_complex_abs(gsl_vector_complex_get(v,i));
    if(tmp > norm)
      norm = tmp;
  }
  return norm;
}

void norm_1d_1(gsl_vector * M, const gsl_matrix_complex * Y){
  assert(M->size <= Y->size2);
  for(unsigned int i=0; i<M->size; i++){
    gsl_vector_complex_const_view cz = gsl_matrix_complex_const_column(Y,i);
    gsl_vector_set(M,i,norm_1d_1(&cz.vector));
  }
}

void norm_1d_inf(gsl_vector * H,  const gsl_matrix_complex * Z){
  assert(H->size <= Z->size1);
  for(unsigned int i=0; i<H->size; i++){
    gsl_vector_complex_const_view cz = gsl_matrix_complex_const_row(Z,i);
    gsl_vector_set(H,i,norm_1d_inf(&cz.vector));
  }
}

void gsl_matrix_elementary_set(gsl_matrix *X, std::vector<unsigned int>& ind){
  assert(X->size2 <= ind.size());
  for(unsigned int j=0; j<X->size2; j++){
    gsl_vector_view cx = gsl_matrix_column(X,j);
    gsl_vector_set_zero(&cx.vector);
    gsl_vector_set(&cx.vector,ind[j],1.);
  }
}

void gsl_matrix_elementary_set(gsl_matrix_complex *X, std::vector<unsigned int>& ind){
  assert(X->size2 <= ind.size());
  for(unsigned int j=0; j<X->size2; j++){
    gsl_vector_complex_view cx = gsl_matrix_complex_column(X,j);
    gsl_vector_complex_set_zero(&cx.vector);
    gsl_vector_complex_set(&cx.vector,ind[j],GSL_COMPLEX_ONE);
  }
}

void sign_round_up(gsl_matrix * S, const gsl_matrix_complex * Y){
  for(unsigned int i=0; i<Y->size1; i++)
    for(unsigned int j=0; j<Y->size2; j++)
      gsl_matrix_set(S,i,j,(GSL_REAL(gsl_matrix_complex_get(Y,i,j)) >=0) ? 1 : -1);
}

void sign_round_up(gsl_matrix_complex * S, const gsl_matrix_complex * Y){
  for(unsigned int i=0; i<Y->size1; i++)
    for(unsigned int j=0; j<Y->size2; j++)
      gsl_matrix_complex_set(S,i,j,(GSL_REAL(gsl_matrix_complex_get(Y,i,j)) >=0) ? GSL_COMPLEX_ONE : gsl_complex_rect(-1,0.));
}

double one_normest_product(const gsl_matrix_complex *A,const gsl_matrix_complex *B){
  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder product;
  product.reset(A->size1,A->size2);
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,B,A,GSL_COMPLEX_ZERO,product);
  double normest = one_normest_core(product);
  return normest;
}
  
struct gsl_rng_holder{
  const gsl_rng_type* T;
  gsl_rng* r;
  gsl_rng_holder():T(gsl_rng_env_setup()),r(gsl_rng_alloc(T)){}
  ~gsl_rng_holder(){ gsl_rng_free(r); }
  operator gsl_rng*(){ return r; }
};
  
struct one_normest_core_scratch_space{
  gsl_matrix_complex* X;
  gsl_matrix_complex* Y;
  gsl_matrix_complex* Z;
  gsl_matrix_complex* S;
  gsl_matrix_complex* S_old;
  gsl_vector_complex* w;
  gsl_vector* h;
  gsl_vector* mags;
  std::vector<unsigned int> ind;
  std::vector<unsigned int> ind_hist;
  std::vector<unsigned int> used_enties;
  std::vector<unsigned int> unused_enties;
  
  one_normest_core_scratch_space():
  X(nullptr),Y(nullptr),Z(nullptr),S(nullptr),S_old(nullptr),
  w(nullptr),h(nullptr),mags(nullptr){}
  
  ~one_normest_core_scratch_space(){
    gsl_matrix_complex_free(X);
    gsl_matrix_complex_free(Y);
    gsl_matrix_complex_free(Z);
    gsl_matrix_complex_free(S);
    gsl_matrix_complex_free(S_old);
    gsl_vector_complex_free(w);
    gsl_vector_free(h);
    gsl_vector_free(mags);
  }
  
  void reset(unsigned int n, unsigned int t){
    //matrices are all the same size and can be grouped together
    if(!X || n!=X->size1 || t!=X->size2){
      if(X){
        gsl_matrix_complex_free(X);
        gsl_matrix_complex_free(Y);
        gsl_matrix_complex_free(Z);
        gsl_matrix_complex_free(S);
        gsl_matrix_complex_free(S_old);
      }
      X = gsl_matrix_complex_alloc(n,t);
      Y = gsl_matrix_complex_alloc(n,t);
      Z = gsl_matrix_complex_alloc(n,t);
      S = gsl_matrix_complex_alloc(n,t);
      S_old = gsl_matrix_complex_alloc(n,t);
    }
    //w and h are the same size
    if(!w || n!=w->size){
      if(w){
        gsl_vector_complex_free(w);
        gsl_vector_free(h);
      }
      w = gsl_vector_complex_alloc(n);
      h = gsl_vector_alloc(n);
    }
    if(!mags || t!=mags->size){
      if(mags)
        gsl_vector_free(mags);
      mags = gsl_vector_alloc(t);
    }
    ind.resize(n);
    ind_hist.clear();
    ind_hist.reserve(n);
    used_enties.reserve(n);
    unused_enties.reserve(n);
  }
};

double one_normest_core(const gsl_matrix_complex *A, unsigned int t, unsigned int itmax){
  assert(A->size1 == A->size2);
  if ( itmax < 2 )
    throw std::runtime_error("At least two iterations are needed.");
  if (t < 1 )
    throw std::runtime_error("At least one column is needed.");
  if (t >= A->size1)
    throw std::runtime_error("t should be smaller than the order of the matrix.");

  unsigned int n = A->size1;
  unsigned int nmults = 0;
  unsigned int nresamples = 0;
  
  SQUIDS_THREAD_LOCAL gsl_rng_holder r;
  SQUIDS_THREAD_LOCAL one_normest_core_scratch_space scratch;
  scratch.reset(n,t);
  gsl_matrix_complex* X=scratch.X;
  gsl_matrix_complex* Y=scratch.Y;
  gsl_matrix_complex* Z=scratch.Z;
  gsl_matrix_complex* S=scratch.S;
  gsl_matrix_complex* S_old=scratch.S_old;
  gsl_vector_complex* w=scratch.w;
  gsl_vector* h=scratch.h;
  gsl_vector* mags=scratch.mags;
  std::vector<unsigned int>& ind=scratch.ind;
  std::vector<unsigned int>& ind_hist=scratch.ind_hist;
  std::vector<unsigned int>& used_enties=scratch.used_enties;
  std::vector<unsigned int>& unused_enties=scratch.unused_enties;

  gsl_matrix_complex_set_all(X,GSL_COMPLEX_ONE);
  gsl_matrix_complex_set_all(S,GSL_COMPLEX_ZERO);

  // this produces random X guys
  if (t > 1){
      //gsl_matrix_complex_print(X);
    for(unsigned int i = 1; i < t; i++){
      resample_column(i,X,r);
      //gsl_matrix_complex_print(X);
    }

    for(unsigned int i = 0; i < t; i++){
      while (column_needs_resampling(i,X)) {
        resample_column(i,X,r);
        nresamples++;
      }
    }
  }
  //gsl_matrix_scale(X,1./n);
  gsl_matrix_complex_scale(X,gsl_complex_rect(1./n,0.));
  double est_old=0;
  unsigned int k=0;
  unsigned int ind_best=0;
  for(unsigned int i=0; i<n; i++)
    ind[i] = i;

  double est; unsigned int best_j;
  //std::cout << "start loop of death" << std::endl;
  while (true){
    //std::cout << "k " << k << std::endl;
    //std::cout << "being iteration" << std::endl;
    //gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,X,0.0,Y);
    //std::cout << "A" << std::endl;
    //gsl_matrix_complex_print(A);
    //std::cout << "X" << std::endl;
    //gsl_matrix_complex_print(X);
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,A,X,GSL_COMPLEX_ZERO,Y);
    //std::cout << "Y" << std::endl;
    //gsl_matrix_complex_print(Y);
    nmults++;
    norm_1d_1(mags,Y);
    best_j = gsl_vector_max_index(mags);
    est = gsl_vector_get(mags,best_j);

    //(1)
    if(est > est_old or k == 2){
      if (k>=2)
        ind_best = ind[best_j];
      gsl_vector_complex_const_view cy = gsl_matrix_complex_const_column(Y,best_j);
      //gsl_vector_memcpy(w,&cy.vector);
      gsl_vector_complex_memcpy(w,&cy.vector);
    }

    //std::cout << "est " << est << " best_j " << best_j << std::endl;

    if(est <= est_old and k >=2){
      est = est_old;
      break;
    }

    if (k>itmax){
      break;
    }
    // save previous step
    est_old = est;
    //gsl_matrix_memcpy(S_old,S);
    gsl_matrix_complex_memcpy(S_old,S);
    // new S
    sign_round_up(S,Y);
    //(2)
    //std::cout << "S_old" << std::endl;
    //gsl_matrix_complex_print(S_old);
    //std::cout << "S" << std::endl;
    //gsl_matrix_complex_print(S);
    //std::cout << (every_col_is_parallel(S,S_old) ? "true" : "false") << std::endl;
    if (every_col_is_parallel(S,S_old)){
      break;
    }
    //std::cout << "end e c p c" << std::endl;

    if(t>1){
      for(unsigned int i=0; i<t; i++){
        while(column_needs_resampling(i,S,S_old)){
          resample_column(i,S,r);
          nresamples++;
        }
      }
    }

    //(3)
    //gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,S,0.0,Z);
    gsl_blas_zgemm(CblasTrans,CblasNoTrans,GSL_COMPLEX_ONE,A,S,GSL_COMPLEX_ZERO,Z);
    nmults++;
    norm_1d_inf(h,Z);
    //std::cout << "h" << std::endl;
    //gsl_vector_print(h);
    //std::cout << "Z" << std::endl;
    //gsl_matrix_complex_print(Z);

    //(4)
    if ( k>=2 and gsl_vector_max_index(h) == ind_best )
      break;

    // reorder options
    // this might be ineficietn, but lets try it
    std::vector<std::pair<double,unsigned int>> pares(n);
    for(unsigned int ii = 0; ii < n; ii++)
      pares[ii] = std::pair<double,unsigned int>{gsl_vector_get(h,ii),ind[ii]};
    std::sort(pares.begin(),pares.end());
    std::reverse(pares.begin(),pares.end());
    for(unsigned int ii = 0; ii < n; ii++){
      gsl_vector_set(h,ii,pares[ii].first);
      ind[ii] = pares[ii].second;
    }

    /*
    std::cout << "print " << n << std::endl;
    for(auto tt : pares)
        std::cout << tt.first << " " << tt.second << " ";
    std::cout << std::endl;
    */

    if (t>1){

      //if(condition_magica que depende de ind e ind_hist)
      //  break;

      //(5)
      used_enties.clear();
      unused_enties.clear();
      for(unsigned int ii : ind){
        bool ii_is_in_indhist = false;
        for(unsigned int jj : ind_hist){
          if(ii==jj)
            ii_is_in_indhist = true;
        }
        if (ii_is_in_indhist)
          used_enties.push_back(ii);
        else
          unused_enties.push_back(ii);
      }
      //std::cout << "unused entries size" << std::endl;
      //std::cout << unused_enties.size() << std::endl;

      ind.clear();
      for(auto ii : unused_enties){
        ind.push_back(ii);
      }
      for(auto ii : used_enties){
        ind.push_back(ii);
      }
    }

    /*
    std::cout << "ind" << std::endl;
    for(auto ii:ind)
      std::cout << ii << " ";
    std::cout << std::endl;
    */

    gsl_matrix_elementary_set(X,ind);
    for(unsigned int ii=0; ii < t; ii++){
      ind_hist.push_back(ind[ii]);
    }
    k++;
    //std::cout << "X" << std::endl;
    //gsl_matrix_complex_print(X);
  }

  //std::cout << "est: " << est << std::endl;

  //std::cout << "END CRAZY Function" << std::endl;
  return est;
}

int ell(const gsl_matrix_complex * A,unsigned int m){
 unsigned int p = 2*m+1;
 double abs_c_recip = gsl_sf_choose(2*p,p)*gsl_sf_fact(2*p+1);
 //std::cout << abs_c_recip << std::endl;
 // round up threshold
 double u = pow(2.,-53);

  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder absA;
  absA.reset(A->size1,A->size2);
 gsl_matrix_complex_memcpy(absA,A);
 for(unsigned int i = 0; i < absA->size1; i++){
   for(unsigned int j = 0; j < absA->size2; j++){
      gsl_matrix_complex_set(absA,i,j,
          gsl_complex_rect(gsl_complex_abs(gsl_matrix_complex_get(absA,i,j)),0.));
   }
 }

 double est = one_normest_matrix_power(absA,p);

 if (not est)
   return 0;
 double alpha = est / (exact_1_norm(A) * abs_c_recip);
 double log2_alpha_div_u = log(alpha/u)/M_LN2;
 int value = static_cast<int>(std::ceil(log2_alpha_div_u / (2 * m)));
 return std::max(value,0);
}

void pade3(const gsl_matrix_complex *A, const gsl_matrix_complex *id,
            const gsl_matrix_complex *A2,
            gsl_matrix_complex *U, gsl_matrix_complex *V){
  std::vector<double> b {120., 60., 12., 1.};

  // here we will use V as a temp to calculate U
  gsl_matrix_complex_mul(V,id,gsl_complex_rect(b[1],0.));
  gsl_matrix_complex_add(V,A2,gsl_complex_rect(b[3],0.));
  // this function sets the U matrix value
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,A,V,GSL_COMPLEX_ZERO,U);

  // this lines set the value of V
  gsl_matrix_complex_mul(V,A2,gsl_complex_rect(b[2],0.));
  gsl_matrix_complex_add(V,id,gsl_complex_rect(b[0],0.));

  return;
}

void pade5(const gsl_matrix_complex *A, const gsl_matrix_complex *id,
            const gsl_matrix_complex *A2, const gsl_matrix_complex *A4,
            gsl_matrix_complex *U, gsl_matrix_complex *V){
  std::vector<double> b {30240., 15120., 3360., 420., 30., 1.};

  // here we will use V as a temp to calculate U
  gsl_matrix_complex_mul(V,id,gsl_complex_rect(b[1],0.));
  gsl_matrix_complex_add(V,A2,gsl_complex_rect(b[3],0.));
  gsl_matrix_complex_add(V,A4,gsl_complex_rect(b[5],0.));
  // this function sets the U matrix value
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,A,V,GSL_COMPLEX_ZERO,U);

  // this lines set the value of V
  gsl_matrix_complex_mul(V,id,gsl_complex_rect(b[0],0.));
  gsl_matrix_complex_add(V,A2,gsl_complex_rect(b[2],0.));
  gsl_matrix_complex_add(V,A4,gsl_complex_rect(b[4],0.));

  return;
}

void pade7(const gsl_matrix_complex *A, const gsl_matrix_complex *id,
            const gsl_matrix_complex *A2, const gsl_matrix_complex *A4, const gsl_matrix_complex *A6,
            gsl_matrix_complex *U, gsl_matrix_complex *V){
  std::vector<double> b {17297280., 8648640., 1995840.,
                         277200., 25200., 1512., 56., 1};

  // here we will use V as a temp to calculate U
  gsl_matrix_complex_mul(V,id,gsl_complex_rect(b[1],0.));
  gsl_matrix_complex_add(V,A2,gsl_complex_rect(b[3],0.));
  gsl_matrix_complex_add(V,A4,gsl_complex_rect(b[5],0.));
  gsl_matrix_complex_add(V,A6,gsl_complex_rect(b[7],0.));
  // this function sets the U matrix value
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,A,V,GSL_COMPLEX_ZERO,U);

  // this lines set the value of V
  gsl_matrix_complex_mul(V,id,gsl_complex_rect(b[0],0.));
  gsl_matrix_complex_add(V,A2,gsl_complex_rect(b[2],0.));
  gsl_matrix_complex_add(V,A4,gsl_complex_rect(b[4],0.));
  gsl_matrix_complex_add(V,A6,gsl_complex_rect(b[6],0.));

  return;
}

void pade9(const gsl_matrix_complex *A, const gsl_matrix_complex *id,
            const gsl_matrix_complex *A2, const gsl_matrix_complex *A4, const gsl_matrix_complex *A6,
            gsl_matrix_complex *U, gsl_matrix_complex *V){
  std::vector<double> b {17643225600., 8821612800., 2075673600., 302702400.,
                         30270240.,2162160., 110880., 3960., 90., 1.};

  // here we will use V as a temp to calculate U
  gsl_matrix_complex_mul(V,id,gsl_complex_rect(b[1],0.));
  gsl_matrix_complex_add(V,A2,gsl_complex_rect(b[3],0.));
  gsl_matrix_complex_add(V,A4,gsl_complex_rect(b[5],0.));
  gsl_matrix_complex_add(V,A6,gsl_complex_rect(b[7],0.));
  // this function sets the U matrix value
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,A,V,GSL_COMPLEX_ZERO,U);

  // this lines set the value of V
  gsl_matrix_complex_mul(V,id,gsl_complex_rect(b[0],0.));
  gsl_matrix_complex_add(V,A2,gsl_complex_rect(b[2],0.));
  gsl_matrix_complex_add(V,A4,gsl_complex_rect(b[4],0.));
  gsl_matrix_complex_add(V,A6,gsl_complex_rect(b[6],0.));

  return;
}

void pade13(const gsl_matrix_complex *A, const gsl_matrix_complex *id,
            const gsl_matrix_complex *A2, const gsl_matrix_complex *A4, const gsl_matrix_complex *A6,
            gsl_matrix_complex *U, gsl_matrix_complex *V){
  std::vector<double> b {64764752532480000., 32382376266240000., 7771770303897600.,
                        1187353796428800., 129060195264000., 10559470521600., 670442572800.,
                        33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.};

  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder tmp;
  tmp.reset(A->size1,A->size2);

  // here we will use V as a temp to calculate U
  gsl_matrix_complex_mul(V,A2,gsl_complex_rect(b[9],0.));
  gsl_matrix_complex_add(V,A4,gsl_complex_rect(b[11],0.));
  gsl_matrix_complex_add(V,A6,gsl_complex_rect(b[13],0.));
  // this function sets the U matrix value
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,A6,V,GSL_COMPLEX_ZERO,tmp);

  gsl_matrix_complex_add(tmp,id,gsl_complex_rect(b[1],0.));
  gsl_matrix_complex_add(tmp,A2,gsl_complex_rect(b[3],0.));
  gsl_matrix_complex_add(tmp,A4,gsl_complex_rect(b[5],0.));
  gsl_matrix_complex_add(tmp,A6,gsl_complex_rect(b[7],0.));

  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,A,tmp,GSL_COMPLEX_ZERO,U);

  // this lines set the value of V
  gsl_matrix_complex_mul(V,A2,gsl_complex_rect(b[8],0.));
  gsl_matrix_complex_add(V,A4,gsl_complex_rect(b[10],0.));
  gsl_matrix_complex_add(V,A6,gsl_complex_rect(b[12],0.));

  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,A6,V,GSL_COMPLEX_ZERO,tmp);

  gsl_matrix_complex_mul(V,id,gsl_complex_rect(b[0],0.));
  gsl_matrix_complex_add(V,A2,gsl_complex_rect(b[2],0.));
  gsl_matrix_complex_add(V,A4,gsl_complex_rect(b[4],0.));
  gsl_matrix_complex_add(V,A6,gsl_complex_rect(b[6],0.));
  gsl_matrix_complex_add(V,tmp,GSL_COMPLEX_ONE);

//  gsl_matrix_complex_free(tmp);
  return;
}

struct gsl_permutation_holder{
  gsl_permutation* p;
  gsl_permutation_holder():p(nullptr){}
  ~gsl_permutation_holder(){ gsl_permutation_free(p); }
  void reset(unsigned int n){
    if(!p || n!=p->size){
      if(p)
        gsl_permutation_free(p);
      p = gsl_permutation_alloc(n);
    }
  }
  operator gsl_permutation*(){ return p; }
  gsl_permutation* operator->(){ return p; }
};

void solve_P_Q(const gsl_matrix_complex * U, const gsl_matrix_complex *V,
                gsl_matrix_complex * SPQ){
  // this functions solves the problem (V-U) X = (U+V)
  // not sure what do to here at the moment
  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder P;
  P.reset(U->size1,U->size2);
  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder Q;
  Q.reset(U->size1,U->size2);
  SQUIDS_THREAD_LOCAL gsl_permutation_holder per;
  per.reset(U->size1);
  /*
  std::cout << "V" << std::endl;
  gsl_matrix_complex_print(V);
  std::cout << "U" << std::endl;
  gsl_matrix_complex_print(U);
  */
  gsl_matrix_complex_memcpy(P,V);
  gsl_matrix_complex_memcpy(Q,V);
  gsl_matrix_complex_add(P,U,GSL_COMPLEX_ONE);
  gsl_matrix_complex_sub(Q,U,GSL_COMPLEX_ONE);

  int sign;
  /*
  std::cout << "Q" << std::endl;
  gsl_matrix_complex_print(Q);
  std::cout << "P" << std::endl;
  gsl_matrix_complex_print(P);
  */

  gsl_linalg_complex_LU_decomp(Q,per,&sign);
  for (unsigned int col = 0; col < U->size1; col++){
    gsl_vector_complex_view in_column = gsl_matrix_complex_column(P,col);
    gsl_vector_complex_view out_column = gsl_matrix_complex_column(SPQ,col);
    gsl_linalg_complex_LU_solve(Q,per,&in_column.vector,&out_column.vector);
  }
}

void matrix_exponential(gsl_matrix_complex* eA, const gsl_matrix_complex *A){
  //gsl_matrix_complex_print(A); std::cout << '\n';
  bool isDiagonal=true;
  for(unsigned int i = 0; i<A->size1; i++){
    for(unsigned int j = 0; j<A->size2; j++){
      if(i!=j){
        if(gsl_matrix_complex_get(A,i,j).dat[0] || gsl_matrix_complex_get(A,i,j).dat[1])
          isDiagonal=false;
      }
    }
  }
  if(isDiagonal){
    gsl_matrix_complex_set_all(eA,GSL_COMPLEX_ZERO);
    for(unsigned int i = 0; i<A->size1; i++)
      gsl_matrix_complex_set(eA,i,i,gsl_complex_exp(gsl_matrix_complex_get(A,i,i)));
    return; //done!
  }
  
  // construct the identity matrix
  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder id;
  id.reset(A->size1,A->size2);
  gsl_matrix_complex_set_identity(id);
  // these matrices are always needed
  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder U;
  U.reset(A->size1,A->size2);
  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder V;
  V.reset(A->size1,A->size2);
  // try Pade order 3.
  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder A2;
  A2.reset(A->size1,A->size2);
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,A,A,GSL_COMPLEX_ZERO,A2);
  //gsl_matrix_complex_print(A2);
  double d6 = pow(one_normest_matrix_power(A2,3),1./6.);
  double eta_1 = std::max(pow(one_normest_matrix_power(A2,2),1./4.),d6);
  if (eta_1 < 1.495585217958292e-002 and ell(A, 3) == 0){
    pade3(A,id,A2,U,V);
    solve_P_Q(U,V,eA);
    return;
  }
  // try Pade order 5
  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder A4;
  A4.reset(A->size1,A->size2);
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,A2,A2,GSL_COMPLEX_ZERO,A4);
  double d4 = pow(exact_1_norm(A4),1./4.);
  double eta_2 = std::max(d4,d6);
  if (eta_2 < 2.539398330063230e-001 and ell(A, 5) == 0){
    pade5(A,id,A2,A4,U,V);
    solve_P_Q(U,V,eA);
    return;
  }
  // try Pade order 7 and 9
  //gsl_matrix_complex * A6 = gsl_matrix_complex_alloc(A->size1,A->size2);
  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder A6;
  A6.reset(A->size1,A->size2);
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,A2,A4,GSL_COMPLEX_ZERO,A6);
  d6 = pow(exact_1_norm(A6),1./6.);
  double d8 = pow(one_normest_matrix_power(A4,2),1./8.);
  double eta_3 = std::max(d6,d8);

  if( eta_3 < 9.504178996162932e-001 and ell(A, 7) == 0 ){
    pade7(A,id,A2,A4,A6,U,V);
    solve_P_Q(U,V,eA);
    // free allocated matrices
    return;
  }

  if( eta_3 < 2.097847961257068e+000 and ell(A, 9) == 0 ){
    pade9(A,id,A2,A4,A6,U,V);
    solve_P_Q(U,V,eA);
    // free allocated matrices
    return;
  }

  // use pade order 13
  double d10 = pow(one_normest_product(A4,A6),1./10.);
  double eta_4 = std::max(d8,d10);
  double eta_5 = std::min(eta_3,eta_4);
  double theta_13 = 4.25;
  int u = std::ceil(log(eta_5/theta_13)/M_LN2);
  unsigned int s = std::max(u,0);
  SQUIDS_THREAD_LOCAL gsl_matrix_complex_holder B;
  B.reset(A->size1,A->size2);
  gsl_matrix_complex_memcpy(B,A);
  gsl_matrix_complex_scale(B,gsl_complex_rect(pow(2.,s*(-1.)),0));
  s += ell(B,13);
  //std::cout << "s " << s << std::endl;
  // rescale all matrices
  /*
  std::cout << "A2" << std::endl;
  gsl_matrix_complex_print(A2);
  */
  gsl_matrix_complex_scale(A2,gsl_complex_rect(pow(2.,-2.*s),0.));
  gsl_matrix_complex_scale(A4,gsl_complex_rect(pow(2.,-4.*s),0.));
  gsl_matrix_complex_scale(A6,gsl_complex_rect(pow(2.,-6.*s),0.));

  /*
  std::cout << "B" << std::endl;
  gsl_matrix_complex_print(B);
  std::cout << "A2" << std::endl;
  gsl_matrix_complex_print(A2);
  std::cout << "A4" << std::endl;
  gsl_matrix_complex_print(A4);
  std::cout << "A6" << std::endl;
  gsl_matrix_complex_print(A6);
  */

  pade13(B,id,A2,A4,A6,U,V);
  /*
  std::cout << "U" << std::endl;
  gsl_matrix_complex_print(U);
  std::cout << "V" << std::endl;
  gsl_matrix_complex_print(V);
  */
  solve_P_Q(U,V,eA);

  //std::cout << "eA" << std::endl;
  //gsl_matrix_complex_print(eA);
  // now do repeated squaring
  for(unsigned int i=0; i<s; i++){
    if (i%2 == 0){
      gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,eA,eA,GSL_COMPLEX_ZERO,B);
    } else {
      gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,GSL_COMPLEX_ONE,B,B,GSL_COMPLEX_ZERO,eA);
    }
  }

  if(s%2 != 0)
    gsl_matrix_complex_memcpy(eA,B);

  /*
  std::cout << "eA" << std::endl;
  gsl_matrix_complex_print(eA);
  std::cout << "B" << std::endl;
  gsl_matrix_complex_print(B);
  std::cout << std::endl;
  */

  return;
}

} // close math_detail namespace
} // close squids namespace
