#ifndef __RABI_H
#define __RABI_H




#include <iostream>
#include <float.h>
#include <math.h>
#include <complex>
#include <vector>

#include <SQUIDS.h>

class rabi: public SQUIDS {
 private:
  //Hamiltonian no external field
  SU_vector suH0;
  //dipole
  SU_vector d;

  //Energy difference
  double Delta_E;
  //Lasser Frequency
  double w;
  //Lasser Amplitud
  double A;
  

 public:
  SU_vector d0;
  SU_vector *b0_proj;
  SU_vector *b1_proj;
  SU_vector *evol_b0_proj;
  SU_vector *evol_b1_proj;

  void EvolveProjectors(double t);
  rabi(void){};
  rabi(double a, double b, double c){init(a,b,c);};
  void init(double, double,double );
  SU_vector H0(double);
  SU_vector HI(int ix,double t);

  void set_evol(void);
};

#endif
