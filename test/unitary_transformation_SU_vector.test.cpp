#include <iostream>
#include <SQuIDs/SUNalg.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <iostream>


void gsl_matrix_complex_change_basis_UCMU(gsl_matrix_complex* U, gsl_matrix_complex* M){
  int numneu = U->size1;
  gsl_matrix_complex *U1 = gsl_matrix_complex_alloc(numneu,numneu);
  gsl_matrix_complex *U2 = gsl_matrix_complex_alloc(numneu,numneu);
  gsl_matrix_complex_memcpy(U1,U);
  gsl_matrix_complex_memcpy(U2,U);
  
  gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(numneu,numneu);
  
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


    void gsl_complex_matrix_exponential(gsl_matrix_complex *eA, gsl_matrix_complex *A, int dimx)
  {
    int j,k=0;
    gsl_complex temp;
    gsl_matrix *matreal =gsl_matrix_alloc(2*dimx,2*dimx);
    gsl_matrix *expmatreal =gsl_matrix_alloc(2*dimx,2*dimx);
    //Converting the complex matrix into real one using A=[Areal, Aimag;-Aimag,Areal]
    for (j = 0; j < dimx;j++)
      for (k = 0; k < dimx;k++)
        {
	  temp=gsl_matrix_complex_get(A,j,k);
	  gsl_matrix_set(matreal,j,k,GSL_REAL(temp));
	  gsl_matrix_set(matreal,dimx+j,dimx+k,GSL_REAL(temp));
	  gsl_matrix_set(matreal,j,dimx+k,GSL_IMAG(temp));
	  gsl_matrix_set(matreal,dimx+j,k,-GSL_IMAG(temp));
        }

    gsl_linalg_exponential_ss(matreal,expmatreal,GSL_PREC_DOUBLE);

    double realp;
    double imagp;
    for (j = 0; j < dimx;j++)
      for (k = 0; k < dimx;k++)
        {
	  realp=gsl_matrix_get(expmatreal,j,k);
	  imagp=gsl_matrix_get(expmatreal,j,dimx+k);
	  gsl_matrix_complex_set(eA,j,k,gsl_complex_rect(realp,imagp));
        }
    gsl_matrix_free(matreal);
    gsl_matrix_free(expmatreal);
  }


int main(){
  using squids::SU_vector;
  
  int dim=2;
  int Ngenerators=dim*dim;
  double th=1.456;
  double dl=1.231;
  double t=3.456;

  //Check that the time evolution is unitary
  for(unsigned int i=0;i<Ngenerators;i++){
    for(unsigned int j=0;j<Ngenerators;j++){
      SU_vector v1=SU_vector::Generator(dim,i);
      SU_vector v2=SU_vector::Generator(dim,j);
      double sprod=v1*v2;
      for(unsigned int it=0;it<dim;it++){
        SU_vector v3=SU_vector::Projector(dim,it);
        double out=sprod-v1.Evolve(v3,t)*v2.Evolve(v3,t);
        if(fabs(out)>1e-15){
          std::cout << "Generators: " << i << "  " << j << "  Evolution Operator: " << it <<
            "  \t" << out  << std::endl;
        }
      }
    }
  }

  //Check that the rotations are unitary transformations
  for(unsigned int i=0;i<Ngenerators;i++){
    for(unsigned int j=0;j<Ngenerators;j++){
      SU_vector v1=SU_vector::Generator(dim,i);
      SU_vector v2=SU_vector::Generator(dim,j);
      double sprod=v1*v2;
      for(unsigned int ir=1;ir<dim;ir++){
        for(unsigned int jr=ir+1;jr<dim;jr++){
          double out=sprod-v1.Rotate(ir,jr,th,dl)*v2.Rotate(ir,jr,th,dl);
          if(fabs(out)>1e-15){
            std::cout << "Generators: " << i << "  " << j << "  Rotation: " << ir << "  " << jr <<
              "  \t" << out  << std::endl;
          }
        }
      }
    }
  }



    
  //gsl_matrix_complex *m=gsl_matrix_complex_alloc(2,2);
  //gsl_matrix_complex *em=gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex *mb=gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex *mb2=gsl_matrix_complex_alloc(2,2);

  gsl_complex z00, z01, z10, z11;
  GSL_SET_COMPLEX(&z00, 0, 0); GSL_SET_COMPLEX(&z01, 0, 1);
  GSL_SET_COMPLEX(&z10, 0, 1
		  ); GSL_SET_COMPLEX(&z11, 1, 0);

  gsl_complex zt00, zt01, zt10, zt11;
  GSL_SET_COMPLEX(&zt00, 1, 0); GSL_SET_COMPLEX(&zt01, 0, 1);
  GSL_SET_COMPLEX(&zt10, 1, 0); GSL_SET_COMPLEX(&zt11, 1, 0);

  
  //gsl_matrix_complex_set(m,0,0,z00); gsl_matrix_complex_set(m,0,1,z01);
  //gsl_matrix_complex_set(m,1,0,z10); gsl_matrix_complex_set(m,1,1,z11);

  gsl_matrix_complex_set(mb,0,0,zt00); gsl_matrix_complex_set(mb,0,1,zt01);
  gsl_matrix_complex_set(mb,1,0,zt10); gsl_matrix_complex_set(mb,1,1,zt11);

  //gsl_matrix_complex_set(mb2,0,0,zt00); gsl_matrix_complex_set(mb2,0,1,zt01);
  //gsl_matrix_complex_set(mb2,1,0,zt10); gsl_matrix_complex_set(mb2,1,1,zt11);

  
  //SU_vector v1(mb);
  //SU_vector v2(m);

  
  //gsl_complex_matrix_exponential(em, m, 2);
  
  
  //gsl_matrix_complex_change_basis_UCMU(em, mb);
  
  for (size_t i = 0; i < mb->size1; i++) {
    for (size_t j = 0; j < mb->size2; j++) {
      std::cout << GSL_REAL(gsl_matrix_complex_get(mb, i, j))
		<< "  " << GSL_IMAG(gsl_matrix_complex_get(mb, i, j))  << "  \t "; 
    }
    std::cout << std::endl;
  }

  std::cout << "------------------" << std::endl;
  
  //SU_vector v3=v1.UTransform(0.58*v2);
  

  //gsl_matrix_complex * outmat = gsl_matrix_complex_alloc (dim, dim);
  //gsl_matrix_complex * em2 = gsl_matrix_complex_alloc (dim, dim);
    
  //gsl_complex_matrix_exponential(em2,m,2);    
  //gsl_matrix_complex_change_basis_UCMU(em2, mb2);

  
  SU_vector v3(mb);
  auto m3=v3.GetGSLMatrix();
  
  
  for (size_t i = 0; i < mb->size1; i++) {
    for (size_t j = 0; j < mb->size2; j++) {
      std::cout << GSL_REAL(gsl_matrix_complex_get(m3.get(), i, j))
		<< "  " << GSL_IMAG(gsl_matrix_complex_get(m3.get(), i, j))  << "  \t "; 
    }
    std::cout << std::endl;
  }
  
  
  
  
  // for(unsigned int i=0;i<Ngenerators;i++){
  //   for(unsigned int j=0;j<Ngenerators;j++){
  //     for(unsigned int k=0;k<Ngenerators;k++){
  // 	if(i!=j && j!=k){
  // 	  SU_vector v1=(SU_vector::Generator(dim,i)+SU_vector::Generator(dim,j));
  // 	  SU_vector v2=SU_vector::Generator(dim,k);
  // 	  SU_vector v3=v2.UTransform(v1);
  // 	  SU_vector v4=(v2-v3.UTransform(-v1));
  // 	  auto out=v4.GetComponents();
  // 	  for(int di=0;di<dim*dim;di++){
  // 	    if(fabs(out[di])>1e-13)
  // 	      std::cout <<i << "  " << j << "  " << k  <<"  Component: "<<di<<" --> " << out[di] << std::endl;
  // 	  }
  // 	}
  //     }
  //   }
  // }
  

  return(0);
}
