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
  //Laser Frequency
  double w;
  //Laser Amplitude
  double A;

 public:
  SU_vector d0;
  std::unique_ptr<SU_vector[]> b0_proj;
  std::unique_ptr<SU_vector[]> b1_proj;
  std::unique_ptr<SU_vector[]> evol_b0_proj;
  std::unique_ptr<SU_vector[]> evol_b1_proj;

  void PreDerive(double t);
  rabi(){};
  rabi(double a, double b, double c){init(a,b,c);};
  void init(double, double,double );
  SU_vector H0(double);
  SU_vector HI(int ix,double t);

  void set_evol();
};

#endif
