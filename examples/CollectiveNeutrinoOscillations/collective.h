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
  SU_vector B;
  SU_vector P;

  //self interaction strengh sqrt(2) G_F n_\nu
  double mu;
  double n_nu;
  double w_min, w_max;
 public:
  SU_vector ex,ey,ez;


  void PreDerive(double t);
  collective(void){};
  collective(double mu, double wmin, double wmax, int Nbins){init(mu, wmin, wmax, Nbins);};
  void init(double mu, double wmin, double wmax, int Nbins);
  SU_vector H0(double);
  SU_vector HI(int ix,double t);

  void set_evol(void);
};

#endif
