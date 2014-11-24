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

 public:
  std::unique_ptr<SU_vector[]> b0_proj;
  std::unique_ptr<SU_vector[]> b1_proj;
  std::unique_ptr<SU_vector[]> evol_b0_proj;
  std::unique_ptr<SU_vector[]> evol_b1_proj;

  void PreDerive(double t);
  collective(void){};
  collective(double mu, double wmin, double wmax, int Nbins){init(mu, wmin, wmax, Nbin);};
  void init(double mu, double wmin, double wmax, int Nbins);
  SU_vector H0(double);
  SU_vector HI(int ix,double t);

  void set_evol(void);
};

#endif
