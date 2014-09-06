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
 *      Jordi Salvado (University of Wisconsin Madison)                        *
 *         jsalvado@icecube.wisc.edu                                           *
 ******************************************************************************/



#include "SUNalg.h"

#define SQR(x)      ((x)*(x))                        // x^2
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  // |x|^2
#define POW10(x)    (exp(M_LN10*(x)))                // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( (b) > 0.0 ? (fabs(a)) : (-fabs(a)) )
#define KRONECKER(i,j)  ( (i)==(j) ? 1 : 0 )


/*
 * We implement the SU(N) semialgebraic solution of the
 * diffusion problem. We wil asumme that the problem
 * will be solve in the mass-interaction basis
*/

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


SU_vector::SU_vector(){
  dim = 0;
  size = 0;
  isinit = false;
  isinit_d = false;
};

SU_vector::SU_vector( const SU_vector& V){
  dim=V.dim;
  size=V.size;
  components=new double[size];
  for(int i=0;i<size;i++)
    components[i]=V.components[i];
}

SU_vector::SU_vector(int d,double* comp){
  isinit = false;
  isinit_d = false;
  InitSU_vector(d,comp);
};

SU_vector::SU_vector(int d){
  isinit = false;
  isinit_d = false;
  InitSU_vector(d);
};

SU_vector::SU_vector(std::vector<double> comp){
  isinit = false;
  isinit_d = false;
  InitSU_vector(comp);
};

SU_vector::SU_vector(string Type,int ii,int d){
  isinit = false;
  isinit_d = false; 
  InitSU_vector(Type,ii,d);
};

SU_vector::SU_vector(gsl_matrix_complex* m){
  isinit = false;
  isinit_d = false;  
  InitSU_vector(m);
};


SU_vector::~SU_vector(){
  if(isinit){
    delete [] components;
    isinit=false;
  }
};

/*
-----------------------------------------------------------------------
Initializers
-----------------------------------------------------------------------
*/

void SU_vector::InitSU_vector(int d){
  size = SQR(d);
  dim = d;
  if(isinit_d){
    throw std::runtime_error( "SU_vector::InitSU_vector : Initialization not allowed, vector already initialized by reference");
  }

  if(MAXSIZE<size || size<0){
    throw std::runtime_error("SU_vector::InitSU_vector : Not allowed size, only up to SU(6)");
  }
  components = new double[size];
  isinit=true;
  memset(components, 0.0, sizeof(double)*size);
};

void SU_vector::InitSU_vector(const SU_vector& V){
  dim=V.dim;
  size=V.size;
  if(isinit){
    delete[] components;
  }
  components=new double[size];
  for(int i=0;i<size;i++){
    components[i]=V.components[i];
  }

  isinit = true;
};

void SU_vector::InitSU_vector(int d,double* comp){
  dim = d;
  size = SQR(d);
  if(MAXSIZE<size || size<0){
    throw std::runtime_error("SU_vector::InitSU_vector : Not allowed size, only up to SU(6)");
  }
  if(isinit){
    delete[] components;
  }
  components = comp;
  isinit=false;
  isinit_d=true;
  
};

void SU_vector::InitSU_vector(gsl_matrix_complex* m){
  size = m->size1*m->size1;
  dim=m->size1;
  if(isinit_d){
    throw std::runtime_error( "SU_vector::InitSU_vector : Initialization not allowed, vector already initialized by reference");
  }
  if(isinit){
    delete[] components;
  }

  if(MAXSIZE<size || size<0){
    throw std::runtime_error("SU_vector::InitSU_vector : Not allowed size, only up to SU(6)");
  }
  
  if((dim%(int)sqrt(size))!=0){
    throw std::runtime_error("SU_vector::InitSU_vector : error wrong matrix.");
  }
  components=new double[size];
  isinit=true;
  memset(components, 0.0, sizeof(double)*size);
  
  double m_real[dim][dim]; double m_imag[dim][dim];
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      m_real[i][j] = gsl_matrix_complex_get(m,i,j).dat[0];
      m_imag[i][j] = gsl_matrix_complex_get(m,i,j).dat[1];
    }
  }
  // the following rules are valid ONLY when the initial
  // matrix is hermitian
    
  switch (dim){
  case 2:
#include "MatrixToSU2.txt"
    break;
  case 3:
#include "MatrixToSU3.txt"
    break;
  case 4:
#include "MatrixToSU4.txt"
    break;
  case 5:
#include "MatrixToSU5.txt"
    break;
  default:
    throw std::runtime_error("GLS_MATRIX_COMPLEX to SU_vector :: Error. ");
  }
};


void SU_vector::InitSU_vector(std::vector<double> comp){
  size = (int)comp.size();
  dim = (int)sqrt(comp.size());
  if(isinit_d){
    throw std::runtime_error( "SU_vector::InitSU_vector : Initialization not allowed, vector already initialized by reference");
  }

  if(isinit){
    delete[] components;
  }

  if(MAXSIZE<size || size<0){
    throw std::runtime_error("SU_vector::InitSU_vector : Not allowed size, only up to SU(6)");
  }


  if(dim*dim!=size){
    throw std::runtime_error("SU_vector::InitSU_vector : Error whrong vector size");
  }

  components=new double[size];
  isinit=true;
  for(int i = 0; i< size; i++)
    components[i] = comp[i];
};

void SU_vector::InitSU_vector(string Type,int ii,int d){
  dim = d;
  size = SQR(dim);
  if(isinit_d){
    throw std::runtime_error( "SU_vector::InitSU_vector : Initialization not allowed, vector already initialized by reference");
  }
  if(isinit){
    delete[] components;
    isinit=false;
  }


  if(MAXSIZE<size || size<0){
    throw std::runtime_error("SU_vector::InitSU_vector : Not allowed size, only up to SU(6)");
  }

  if(!isinit){
    components=new double[size];
    memset(components, 0.0, sizeof(double)*size);
    isinit=true;
  }
  
  double m_real[dim][dim]; double m_imag[dim][dim];

  if(Type=="Proj"||Type=="proj"||Type=="Projector"||Type=="projector"){
    if(ii>dim|| ii<0){
      throw std::runtime_error("SU_vector::SU_vector : input parameter out of range for projectors");
    }

    for(int i=0; i<dim; i++){
      for(int j=0; j<dim; j++){
	m_real[i][j] = KRONECKER(i,j)*KRONECKER(i,ii);
	m_imag[i][j] = 0.0;
      }
    }
  }else if(Type == "Identity"){
    for(int i=0; i<dim; i++){
      for(int j=0; j<dim; j++){
	m_real[i][j] = KRONECKER(i,j);
	m_imag[i][j] = 0.0;
      }
    }
  }else if(Type =="PosProjector"){
    if(ii>dim|| ii<0){
      throw std::runtime_error("SU_vector::SU_vector : input parameter out of range for positive projectors ");
      exit(1);
    }
    for(int i=0; i<dim; i++){
      for(int j=0; j<dim; j++){
	if(i<ii)
	  m_real[i][j] = KRONECKER(i,j);
	else
	  m_real[i][j] = 0.0;
	
	m_imag[i][j] = 0.0;
      }
    }
  }else if(Type == "NegProjector"){
    if(ii>dim|| ii<0){
      throw std::runtime_error("SU_vector::SU_vector : input parameter out of range for negative projectors ");
    }
    for(int i=0; i<dim; i++){
      for(int j=0; j<dim; j++){
	if(dim-i<ii)
	  m_real[i][j] = KRONECKER(i,j);
	else
	  m_real[i][j] = 0.0;
    	
	m_imag[i][j] = 0.0;
      }
    }
  }else if(Type == "Component" || Type == "component" || Type  == "Comp" || Type == "comp"
	   || Type == "Generator" || Type == "Gen" || Type == "generator" || Type == "gen"){
    if(ii>size|| ii<0){
      throw std::runtime_error("SU_vector::SU_vector : input parameter out of range for generators");
    }
    components[ii] = 1.0;
    return;
  }else{
    throw std::runtime_error("SU_vector::SU_vector :  Option not available");
  }

  // the following rules are valid ONLY when the initial
  // matrix is hermitian
  switch (dim){
  case 2:
#include "MatrixToSU2.txt"
    break;
  case 3:
#include "MatrixToSU3.txt"
    break;
  case 4:
#include "MatrixToSU4.txt"
    break;
  case 5:
#include "MatrixToSU5.txt"
    break;
  default:
    throw std::runtime_error("GLS_MATRIX_COMPLEX to SU_vector :: Error. ");
  }    
};


/*
-----------------------------------------------------------------------
Operations
-----------------------------------------------------------------------
*/


vector<double> SU_vector::GetComponents(){
  vector<double> x ( dim*dim );
  for ( int i = 0; i < dim*dim ; i ++ )
    x[i] = components[i];
  return x;
}



SU_vector SU_vector::Rescale(double x){
  for(int i=0; i< size; i++){
    components[i] = x*components[i];
  };
  return *this;
};

SU_vector SU_vector::Rotate(int i, int j, double th, double del){
  SU_vector suv=*this;
  SU_vector suv_rot(dim);

  switch (dim){
  case 2:
    if (i == 1 and j == 2){
#include "RotationSU2_12.txt"
    };
    break;
  case 3:
    switch (i){
    case 1:
      switch (j){
      case 2:
#include "RotationSU3_12.txt"
	break;
      case 3:
#include "RotationSU3_13.txt"
	break;
      };
      break;
    case 2:
      switch (j){
      case 3:
#include "RotationSU3_23.txt"
	break;
      };
    };
    break;
  case 4:
    switch (i){
    case 1:
      switch (j){
      case 2:
#include "RotationSU4_12.txt"
	break;
      case 3:
#include "RotationSU4_13.txt"
	break;
      case 4:
#include "RotationSU4_14.txt"
	break;
      };
      break;
    case 2:
      switch (j){
      case 3:
#include "RotationSU4_23.txt"
	break;
      case 4:
#include "RotationSU4_24.txt"
	break;
      };
      break;  
    case 3:
      switch (j){
      case 4:
#include "RotationSU4_34.txt"
	break;
      };
      break;
    }
    break;
  case 5:
    switch (i){
    case 1:
      switch (j){
      case 2:
#include "RotationSU5_12.txt"
	break;
      case 3:
#include "RotationSU5_13.txt"
	break;
      case 4:
#include "RotationSU5_14.txt"
	break;
      case 5:
#include "RotationSU5_15.txt"
	break;                            
      };
      break;
    case 2:
      switch (j){
      case 3:
#include "RotationSU5_23.txt"
	break;
      case 4:
#include "RotationSU5_24.txt"
	break;
      case 5:
#include "RotationSU5_25.txt"
	break;
      };
      break;  
    case 3:
      switch (j){
      case 4:
#include "RotationSU5_34.txt"
	break;
      case 5:
#include "RotationSU5_35.txt"
	break;
      };
      break;
    case 4:
      switch (j){
      case 5:
#include "RotationSU5_45.txt"
	break;
      };
      break;
    }
    break;
  default:
    printf("SUN_rotation error. \n");
    exit(0);
  };
    
  for(int i=0; i < size; i++){
    components[i] = suv_rot.components[i];
  };
    
  return *this;
};

void SU_vector::RotateToB0(const Const* param){
  Const param2;
  param2.th12= -param->th12;
  param2.th13=  -param->th13;
  param2.th23=  -param->th23;
  param2.th14=  -param->th14;
  param2.th24=  -param->th24;
  param2.th34=  -param->th34;
  param2.th15=  -param->th15;
  param2.th25=  -param->th25;
  param2.th35=  -param->th35;
  param2.th45=  -param->th45;
  param2.th16=  -param->th16;
  param2.th26=  -param->th26;
  param2.th36=  -param->th36;
  param2.th46= -param->th46; 
  param2.th56= -param->th56; 
  RotateToB1(&param2);
}

void SU_vector::RotateToB1(const Const* param){
  SU_vector suv(dim);
  SU_vector suv_rot(dim);

  for(int i=0; i < size; i++){
    suv.components[i] = components[i];
  };


  double th12,th13,th23,del12,del13,del23;
  double th14,th24,th34,del14,del24,del34;
  //double th15,th25,th35,th45,del15,del25,del35,del45;
    
  switch (dim){
  case 2:
    th12 = param->th12;
    del12 = 0.0;
#include "RotationSU2.txt"

    break;
  case 3:
    th12 = param->th12;
    del12 = 0.0;
    th13 = param->th13;
    del13 = param->delta1;
    th23 = param->th23;
    del23 = 0.0;
#include "RotationSU3.txt"
    break;
  case 4:
    suv.Rotate(3,4,param->th34,0.0);
    suv.Rotate(2,4,param->th24,0.0);
    suv.Rotate(1,4,param->th14,param->delta2);
    suv.Rotate(2,3,param->th23,0.0);
    suv.Rotate(1,3,param->th13,param->delta1);
    suv.Rotate(1,2,param->th12,0.0);
            
   suv_rot = suv;


    break;
  case 5:

    break;
  default:
    printf("SUN_rotation error. \n");
    exit(1);
  };
  
  for(int i=0; i < SQR(dim); i++){
    components[i] = suv_rot.components[i];
  };

};


double SU_vector::SUTrace(SU_vector* suv1,SU_vector* suv2){
  double gen_trace = 0.0;
  double id_trace = 0.0;
  for(int i=1; i < SQR(suv1->dim); i++){
    gen_trace += (suv1->components[i])*(suv2->components[i]);
  };
  id_trace = (suv1->components[0])*(suv2->components[0])*double(suv1->dim);

  return id_trace+2.0*gen_trace;
};

double SU_vector::SUTrace(SU_vector& suv1,const SU_vector& suv2){
  double gen_trace = 0.0;
  double id_trace = 0.0;

  for(int i=1; i < SQR(suv1.dim); i++){
    gen_trace += (suv1.components[i])*(suv2.components[i]);
  };

  id_trace = (suv1.components[0])*(suv2.components[0])*double(suv1.dim);
  return id_trace+2.0*gen_trace;
};

void SU_vector::SUPrint(void) const {
  for(int i=0; i< SQR(dim); i++){
    std::cout << components[i] << " ";
  };
  std::cout << std::endl;
};


bool SU_vector::operator==(const SU_vector &other){
  bool equal = true;
  for(int i=0; i < SQR(dim); i++){
    if (components[i] != other.components[i]){
      equal = false;
      break;
    }
  };
  if (dim != other.dim)
    equal = false;
    
  return equal;
};

double SU_vector::operator*(const SU_vector &other){
  return SUTrace(*this,other);
};

SU_vector SU_vector::operator*(const double x){
  SU_vector su_new(dim);
  for(int i=0; i < SQR(dim); i++){
    su_new.components[i] = x*components[i];
  };
  return su_new;
};

SU_vector operator*(const double x, const SU_vector &other){
  int dim = other.Dim();
  SU_vector su_new(dim);
  for(int i=0; i < SQR(dim); i++){
    su_new[i] = x*other[i];
  };
  return su_new;
};

SU_vector operator*(const SU_vector &other, const double x){
  return x*other;
};


SU_vector & SU_vector::operator=(const SU_vector &other){
  if(size!=other.size){
    throw std::runtime_error(" SU_vector::operator= Error not equal size ");
  }
  if(this==&other){
    throw std::runtime_error(" SU_vector::operator= invalid self-assignment ");
  }

  dim = other.dim;
  size = other.size;
  
  for(int i=0; i < size; i++){
    components[i] = other.components[i];
  };
  return *this;
};

SU_vector & SU_vector::operator+=(const SU_vector &other){
  for(int i=0; i < SQR(dim); i++){
    components[i] += other.components[i];
  };
  return *this;
};

SU_vector & SU_vector::operator-=(const SU_vector &other){
  for(int i=0; i < SQR(dim); i++){
    components[i] -= other.components[i];
  };
  return *this;
};

SU_vector SU_vector::operator+(const SU_vector &other){
  SU_vector su_new(dim);
  for(int i=0; i < SQR(dim); i++){
    su_new.components[i] = components[i] + other.components[i];
  };
  return su_new;
};

SU_vector SU_vector::operator-(const SU_vector &other){
  SU_vector su_new(dim);
 
  for(int i=0; i < SQR(dim); i++){
    su_new.components[i] = components[i] - other.components[i];
  };
 
  return su_new;
};


SU_vector SU_vector::SUEvolve(SU_vector & suv1, double t){
  SU_vector suv2(dim);
  for(int i=0; i < size; i++){
    suv2.components[i] = components[i];
  };
  SU_vector suv_new(dim);

  // here we will asumme that suv1 is the
  // evolution, suv1 should be diagonal(linear combination of the SU(n) weights). 
  // The given expressions are only valid when suv1 is a projector-like SU-vector.
  
  switch (dim){
  case 2:
#include "EvolutionSU2.txt"
    break;
  case 3:
#include "EvolutionSU3.txt"
    break;
  case 4:
#include "EvolutionSU4.txt"
    break;
  case 5:
#include "EvolutionSU5.txt"
    break;
  default:
    throw std::runtime_error("SUEvolution :: Error : dim  ");
  };
  return suv_new;
}


ostream& operator<<(ostream& os, const SU_vector& V){
  for(int i=0; i< V.size-1; i++)
    os << V.components[i] << "  ";
  os <<V.components[V.size-1];
  return os;
}


SU_vector SU_alg::iCommutator(const SU_vector& suv1,const SU_vector& suv2){
  if(suv1.dim!=dim || suv2.dim!=dim){ 
    throw std::runtime_error("SU_alg: Commutator error, not right dimensions ");
  }
  SU_vector suv_new(dim);

  switch (dim){
  case 2:
#include "iConmutatorSU2.txt"
    break;
  case 3:
#include "iConmutatorSU3.txt"
    break;
  case 4:
#include "iConmutatorSU4.txt"
    break;
  case 5:
#include "iConmutatorSU5.txt"
    break;
  default:
    throw std::runtime_error("SUiConmutator :: Error.");
    exit(1);
  };
  return suv_new;
};

SU_vector SU_alg::ACommutator(const SU_vector& suv1,const SU_vector& suv2){
  if(suv1.dim!=dim || suv2.dim!=dim){ 
    throw std::runtime_error("SU_alg: Anti Commutator error, not right dimensions ");
  }
  SU_vector suv_new(dim);
  switch (dim){
  case 2:
#include "AnticonmutatorSU2.txt"
    break;
  case 3:
#include "AnticonmutatorSU3.txt"
    break;
  case 4:
#include "AnticonmutatorSU4.txt"
    break;
  case 5:
#include "AnticonmutatorSU5.txt"
    break;
  default:
    throw std::runtime_error("SUAnticonmutator :: Error.");
    exit(1);
  };
  return suv_new;
};
