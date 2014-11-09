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
  if(size>MAXSIZE)
    throw std::runtime_error("SU_vector::SU_vector(unsigned int, double*): Invalid size: only up to SU(6) is supported");
};

SU_vector::SU_vector(unsigned int d):
dim(d),
size(dim*dim),
components(new double[size]),
isinit(true),
isinit_d(false)
{
  if(size>MAXSIZE)
    throw std::runtime_error("SU_vector::SU_vector(unsigned int): Invalid size: only up to SU(6) is supported");
  std::fill(components,components+size,0.0);
};

SU_vector::SU_vector(std::vector<double> comp):
dim(sqrt(comp.size())),
size(comp.size()),
components(new double[size]),
isinit(true),
isinit_d(false)
{
  if(dim*dim!=size)
    throw std::runtime_error("SU_vector::SU_vector(std::vector<double>): Vector size must be a square");
  if(size>MAXSIZE)
    throw std::runtime_error("SU_vector::SU_vector(std::vector<double>): Invalid size: only up to SU(6) is supported");
  std::copy(comp.begin(),comp.end(),components);
};

SU_vector::SU_vector(gsl_matrix_complex* m):
dim(m->size1),
size(dim*dim),
components(new double[size]),
isinit(true),
isinit_d(false)
{
  if(m->size1!=m->size2)
    throw std::runtime_error("SU_vector::SU_vector(gsl_matrix_complex*): Matrix must be square");
  if(size>MAXSIZE)
    throw std::runtime_error("SU_vector::SU_vector(gsl_matrix_complex*): Invalid size: only up to SU(6) is supported");
  
  std::fill(components,components+size,0.0);
  
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

SU_vector::~SU_vector(){
  if(isinit){
    delete[] components;
    isinit=false;
  }
};



void SU_vector::SetBackingStore(double* comp){
  if(isinit)
    delete[] components;
  components = comp;
  isinit=false;
  isinit_d=true;
}

namespace{
  struct sq_array_2D{
    unsigned int d;
    double* data;
    //sq_array_2D(unsigned int d):d(d){}
    
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
    case 6:
#include "MatrixToSU6.txt"
      break;
    default:
      throw std::runtime_error("GLS_MATRIX_COMPLEX to SU_vector :: Error. ");
  }
}
}

SU_vector SU_vector::Projector(unsigned int d, unsigned int ii){
  if(d>sqrt(MAXSIZE))
    throw std::runtime_error("SU_vector::Projector(unsigned int, unsigned int): Invalid size: only up to SU(6) is supported");
  if(ii>=d)
    throw std::runtime_error("SU_vector::Projector(unsigned int, unsigned int): Invalid component: must be smaller than dimension");
  
  double m_real[d][d]; double m_imag[d][d];
  for(int i=0; i<d; i++){
    for(int j=0; j<d; j++){
      m_real[i][j] = KRONECKER(i,j)*KRONECKER(i,ii);
      m_imag[i][j] = 0.0;
    }
  }
  
  SU_vector v(d);
  ComponentsFromMatrices(v.components,d,sq_array_2D{d,m_real[0]},sq_array_2D{d,m_imag[0]});
  return(v);
}

SU_vector SU_vector::Identity(unsigned int d){
  if(d>sqrt(MAXSIZE))
    throw std::runtime_error("SU_vector::Identity(unsigned int): Invalid size: only up to SU(6) is supported");
  
  double m_real[d][d]; double m_imag[d][d];
  for(int i=0; i<d; i++){
    for(int j=0; j<d; j++){
      m_real[i][j] = KRONECKER(i,j);
      m_imag[i][j] = 0.0;
    }
  }
  
  SU_vector v(d);
  ComponentsFromMatrices(v.components,d,sq_array_2D{d,m_real[0]},sq_array_2D{d,m_imag[0]});
  return(v);
}

SU_vector SU_vector::PosProjector(unsigned int d, unsigned int ii){
  if(d>sqrt(MAXSIZE))
    throw std::runtime_error("SU_vector::PosProjector(unsigned int, unsigned int): Invalid size: only up to SU(6) is supported");
  if(ii>=d)
    throw std::runtime_error("SU_vector::PosProjector(unsigned int, unsigned int): Invalid component: must be smaller than dimension");
  
  double m_real[d][d]; double m_imag[d][d];
  for(int i=0; i<d; i++){
    for(int j=0; j<d; j++){
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
  if(d>sqrt(MAXSIZE))
    throw std::runtime_error("SU_vector::NegProjector(unsigned int, unsigned int): Invalid size: only up to SU(6) is supported");
  if(ii>=d)
    throw std::runtime_error("SU_vector::NegProjector(unsigned int, unsigned int): Invalid component: must be smaller than dimension");
  
  double m_real[d][d]; double m_imag[d][d];
  for(int i=0; i<d; i++){
    for(int j=0; j<d; j++){
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

SU_vector SU_vector::Component(unsigned int d, unsigned int ii){
  if(d>sqrt(MAXSIZE))
    throw std::runtime_error("SU_vector::Component(unsigned int, unsigned int): Invalid size: only up to SU(6) is supported");
  if(ii>=d)
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
  for(int i = 0; i< size; i++)
    components[i] = x;
}


vector<double> SU_vector::GetComponents() const{
  vector<double> x ( dim*dim );
  for ( int i = 0; i < dim*dim ; i ++ )
    x[i] = components[i];
  return x;
}



SU_vector SU_vector::Rescale(double x){
  for(int i=0; i< size; i++){
    components[i] = x*components[i];
  }
  return *this;
};

SU_vector SU_vector::Rotate(unsigned int i, unsigned int j, double th, double del){
  SU_vector suv=*this;
  SU_vector suv_rot(dim);

  assert(i<j && "Components selected for rotation must be in ascending order");
  
  switch (dim){
  case 2:
    if (i == 1 and j == 2){
#include "RotationSU2_12.txt"
    }
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
      }
      break;
    case 2:
      switch (j){
      case 3:
#include "RotationSU3_23.txt"
	break;
      }
    }
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
      }
      break;
    case 2:
      switch (j){
      case 3:
#include "RotationSU4_23.txt"
	break;
      case 4:
#include "RotationSU4_24.txt"
	break;
      }
      break;  
    case 3:
      switch (j){
      case 4:
#include "RotationSU4_34.txt"
	break;
      }
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
      }
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
      }
      break;  
    case 3:
      switch (j){
      case 4:
#include "RotationSU5_34.txt"
	break;
      case 5:
#include "RotationSU5_35.txt"
	break;
      }
      break;
    case 4:
      switch (j){
      case 5:
#include "RotationSU5_45.txt"
	break;
      }
      break;
    }
    break;
  default:
    throw std::runtime_error("SUN_rotation error. \n");
  }
    
  for(int i=0; i < size; i++){
    components[i] = suv_rot.components[i];
  }
    
  return *this;
};

void SU_vector::RotateToB0(const Const& param){
  Const param2;
  param2.th12= -param.th12;
  param2.th13= -param.th13;
  param2.th23= -param.th23;
  param2.th14= -param.th14;
  param2.th24= -param.th24;
  param2.th34= -param.th34;
  param2.th15= -param.th15;
  param2.th25= -param.th25;
  param2.th35= -param.th35;
  param2.th45= -param.th45;
  param2.th16= -param.th16;
  param2.th26= -param.th26;
  param2.th36= -param.th36;
  param2.th46= -param.th46; 
  param2.th56= -param.th56; 
  RotateToB1(param2);
}

void SU_vector::RotateToB1(const Const& param){
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
    th12 = param.th12;
    del12 = 0.0;
#include "RotationSU2.txt"

    break;
  case 3:
    th12 = param.th12;
    del12 = 0.0;
    th13 = param.th13;
    del13 = param.delta1;
    th23 = param.th23;
    del23 = 0.0;
#include "RotationSU3.txt"
    break;
  case 4:
    suv.Rotate(3,4,param.th34,0.0);
    suv.Rotate(2,4,param.th24,0.0);
    suv.Rotate(1,4,param.th14,param.delta2);
    suv.Rotate(2,3,param.th23,0.0);
    suv.Rotate(1,3,param.th13,param.delta1);
    suv.Rotate(1,2,param.th12,0.0);
            
   suv_rot = suv;


    break;
  case 5:

    break;
  default:
    printf("SUN_rotation error. \n");
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


bool SU_vector::operator==(const SU_vector& other){
  if(dim != other.dim) //vectors of different sizes cannot be equal
    return false;
  //if one vector contains data and the other does not they cannot be equal
  if((isinit || isinit_d) != (other.isinit || other.isinit_d))
    return false;
  //if both vectors are empty, consider them equal
  if(!isinit && !isinit_d)
    return true;
  //test whether all components are equal
  for(int i=0; i < SQR(dim); i++){
    if (components[i] != other.components[i])
      return false;
  };
  return true;
};

double SU_vector::operator*(const SU_vector &other){
  if(size!=other.size)
    throw std::runtime_error("Non-matching dimensions in SU_vector inner product");
  return SUTrace(*this,other);
};

detail::MultiplicationProxy SU_vector::operator*(double x) const{
  return(detail::MultiplicationProxy{*this,x});
}

detail::MultiplicationProxy operator*(double x, const SU_vector& v){
  return(detail::MultiplicationProxy{v,x});
}


SU_vector & SU_vector::operator=(const SU_vector& other){
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
};

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

SU_vector & SU_vector::operator+=(const SU_vector &other){
  if(size!=other.size)
    throw std::runtime_error("Non-matching dimensions in SU_vector increment");
  for(int i=0; i < SQR(dim); i++){
    components[i] += other.components[i];
  };
  return *this;
};

SU_vector & SU_vector::operator-=(const SU_vector &other){
  if(size!=other.size)
    throw std::runtime_error("Non-matching dimensions in SU_vector decrement");
  for(int i=0; i < SQR(dim); i++){
    components[i] -= other.components[i];
  };
  return *this;
};

detail::AdditionProxy SU_vector::operator+(const SU_vector& other) const{
  if(size!=other.size)
    throw std::runtime_error("Non-matching dimensions in SU_vector addition");
  return(detail::AdditionProxy{*this,other});
}

detail::SubtractionProxy SU_vector::operator-(const SU_vector& other) const{
  if(size!=other.size)
    throw std::runtime_error("Non-matching dimensions in SU_vector subtraction");
  return(detail::SubtractionProxy{*this,other});
}

detail::EvolutionProxy SU_vector::SUEvolve(const SU_vector& suv1,double t) const{
  return(detail::EvolutionProxy{suv1,*this,t});
}

ostream& operator<<(ostream& os, const SU_vector& V){
  for(int i=0; i< V.size-1; i++)
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
