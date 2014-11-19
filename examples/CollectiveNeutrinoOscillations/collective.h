#ifndef __COLLECTIVE_H
#define __COLLECTIVE_H




#include <iostream>
#include <float.h>
#include <math.h>
#include <complex>
#include <vector>

#include <SQUIDS.h>

class collective: public SQUIDS {
 private:
  //Hamiltonian
  SU_vector H0;

  //
  double mu;
  double n_nu;

 public:
  SU_vector d0;
  std::unique_ptr<SU_vector[]> b0_proj;
  std::unique_ptr<SU_vector[]> b1_proj;
  std::unique_ptr<SU_vector[]> evol_b0_proj;
  std::unique_ptr<SU_vector[]> evol_b1_proj;

  void PreDerive(double t);
  collective(void){};
  collective(double a, double b, double c){init(a,b,c);};
  void init(double, double,double );
  SU_vector H0(double);
  SU_vector HI(int ix,double t);

  void set_evol(void);
};

#endif
