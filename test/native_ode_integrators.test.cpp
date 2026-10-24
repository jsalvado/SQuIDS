// Every native Runge-Kutta method and step controller of SQuIDS::Evolve against
// an exact solution, with a time-dependent Hamiltonian so that the stage times
// matter. Prints nothing when all checks pass.
#include <cmath>
#include <iostream>
#include <SQuIDS/SQuIDS.h>

using squids::SU_vector;

// H(t) = w(t) P0 in SU(3): diagonal at all times, so the exact solution is the
// rotation by the accumulated phase Phi(t) = int_0^t w.
static double w(double t){ return 1.0 + 0.5*std::sin(t); }
static double Phi(double t){ return t + 0.5*(1.0 - std::cos(t)); }

class Rotor : public squids::SQuIDS {
  SU_vector P0;
public:
  Rotor() : SQuIDS(2, 3, 1, 0, 0.0), P0(SU_vector::Projector(3, 0)) {
    Set_CoherentRhoTerms(true);
    for(unsigned int ix = 0; ix < nx; ix++){
      state[ix].rho[0] = SU_vector::Generator(3, 1)*(1.0 + ix);
      state[ix].rho[0] += SU_vector::Generator(3, 4)*0.7;
      state[ix].rho[0] += SU_vector::Identity(3)*0.3;
    }
  }
  SU_vector HI(unsigned int ix, unsigned int irho, double t) const {
    return P0*((1.0 + 0.5*ix)*w(t));
  }
  // max deviation from the exact solution at the current time
  double Deviation(double T) const {
    double d = 0;
    for(unsigned int ix = 0; ix < nx; ix++){
      SU_vector r0 = SU_vector::Generator(3, 1)*(1.0 + ix);
      r0 += SU_vector::Generator(3, 4)*0.7;
      r0 += SU_vector::Identity(3)*0.3;
      SU_vector exact = r0.Evolve(P0, -(1.0 + 0.5*ix)*Phi(T)); // Evolve(op,t) is the Heisenberg rotation
      for(unsigned int k = 0; k < 9; k++)
        d = std::max(d, std::abs(exact[k] - state[ix].rho[0][k]));
    }
    return d;
  }
};

int main(){
  using squids::ODEMethod;
  using squids::StepController;
  const ODEMethod methods[] = {ODEMethod::RKF45, ODEMethod::BS3, ODEMethod::DP5, ODEMethod::Tsit5,
                               ODEMethod::Vern6, ODEMethod::Vern7, ODEMethod::Vern8, ODEMethod::Vern9,
                               ODEMethod::RKV76, ODEMethod::DOP853};
  const StepController ctls[] = {StepController::I, StepController::PI, StepController::H211b};
  const double T = 25.0;
  for(ODEMethod m : methods){
    for(StepController c : ctls){
      // one call over the whole interval
      Rotor a;
      a.Set_ODEMethod(m); a.Set_StepController(c);
      a.Set_rel_error(1e-11); a.Set_abs_error(1e-11);
      a.Evolve(T);
      if(a.Deviation(T) > 1e-8 || std::abs(a.Get_t() - T) > 1e-14)
        std::cout << "method " << static_cast<int>(m) << " controller " << static_cast<int>(c)
                  << ": deviation " << a.Deviation(T) << " t=" << a.Get_t() << std::endl;
      // the same interval in five pieces (warm-started step size)
      Rotor b;
      b.Set_ODEMethod(m); b.Set_StepController(c);
      b.Set_rel_error(1e-11); b.Set_abs_error(1e-11);
      for(int i = 0; i < 5; i++) b.Evolve(T/5);
      if(b.Deviation(T) > 1e-8)
        std::cout << "method " << static_cast<int>(m) << " controller " << static_cast<int>(c)
                  << " in pieces: deviation " << b.Deviation(T) << std::endl;
      if(a.Get_IntegrationStats().rhs_evaluations == 0)
        std::cout << "method " << static_cast<int>(m) << ": no evaluations counted" << std::endl;
    }
  }
  // fixed-step mode
  Rotor f;
  f.Set_ODEMethod(ODEMethod::Vern9); f.Set_AdaptiveStep(false); f.Set_NumSteps(400);
  f.Evolve(T);
  if(f.Deviation(T) > 1e-9)
    std::cout << "fixed step: deviation " << f.Deviation(T) << std::endl;
  return 0;
}
