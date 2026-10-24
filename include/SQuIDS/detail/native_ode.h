#ifndef SQUIDS_DETAIL_NATIVE_ODE_H
#define SQUIDS_DETAIL_NATIVE_ODE_H

// Native adaptive explicit Runge-Kutta integration for SQuIDS.
//
// This is an alternative to the GSL odeiv2 driver. Compared to driving the same
// methods through GSL it adds: first-same-as-last reuse, a PI or H211b step
// size controller, an initial step size estimate that does not depend on the
// time unit, and warm starting when Evolve is called repeatedly on contiguous
// intervals. The coefficients live in rk_tableaus.h and are generated (and
// checked against all order conditions) by tools/gen_rk_tableaus.py.

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "rk_tableaus.h"

namespace squids{

///\brief Integration method used by SQuIDS::Evolve.
///
/// GSL uses the GSL odeiv2 driver with the stepper set by Set_GSL_step (the
/// historical behavior and the default). All other values select a native
/// embedded explicit Runge-Kutta pair. The number is the order of the
/// propagated solution.
enum class ODEMethod : int {
  GSL = 0,
  RKF45,   ///< Fehlberg 4(5), same coefficients as gsl_odeiv2_step_rkf45
  BS3,     ///< Bogacki-Shampine 3(2), FSAL
  DP5,     ///< Dormand-Prince 5(4), FSAL
  Tsit5,   ///< Tsitouras 5(4), FSAL
  Vern6,   ///< Verner "efficient" 6(5), FSAL
  Vern7,   ///< Verner "efficient" 7(6)
  Vern8,   ///< Verner "efficient" 8(7)
  Vern9,   ///< Verner "efficient" 9(8)
  RKV76,   ///< Verner RKV76.IIa 7(6)
  DOP853   ///< Hairer DOP853, 8th order with combined 5th/3rd order error estimate, FSAL
};

///\brief Step size controller used by the native integrator.
///
/// I is the default: on nuSQuIDS propagation through the Earth it needed less
/// work than PI or H211b at equal accuracy (benchmarks/ode in nuSQuIDS).
enum class StepController : int {
  I,     ///< classic integral controller, h *= 0.9 (1/err)^(1/k)
  PI,    ///< proportional-integral controller (Gustafsson; Soderlind's PI.4.2 gains)
  H211b  ///< Soderlind's H211b digital filter, b = 4
};

///\brief Counters for the work done by SQuIDS::Evolve.
struct IntegrationStats {
  unsigned long rhs_evaluations;
  unsigned long accepted_steps;
  unsigned long rejected_steps;
  IntegrationStats():rhs_evaluations(0),accepted_steps(0),rejected_steps(0){}
};

namespace detail{

struct ExplicitRK {
  const char* name;
  unsigned int stages;
  unsigned int order;      // order of the propagated solution
  unsigned int err_order;  // order of the embedded solution
  bool fsal;               // last stage is f(t+h, y_new)
  const double* c;
  const double* A;         // stages x stages, row major
  const double* b;
  const double* e;         // error = h * sum_j e_j k_j
  const double* e3;        // DOP853 only: second (3rd order) error estimate
};

inline ExplicitRK GetTableau(ODEMethod m){
  using namespace rk;
  switch(m){
    case ODEMethod::RKF45:  return {"rkf45", 6, 5, 4, false, rkf45_c, rkf45_A, rkf45_b, rkf45_e, nullptr};
    case ODEMethod::BS3:    return {"bs3", 4, 3, 2, true, bs3_c, bs3_A, bs3_b, bs3_e, nullptr};
    case ODEMethod::DP5:    return {"dp5", 7, 5, 4, true, dp5_c, dp5_A, dp5_b, dp5_e, nullptr};
    case ODEMethod::Tsit5:  return {"tsit5", 7, 5, 4, true, tsit5_c, tsit5_A, tsit5_b, tsit5_e, nullptr};
    case ODEMethod::Vern6:  return {"vern6", 9, 6, 5, true, vern6_c, vern6_A, vern6_b, vern6_e, nullptr};
    case ODEMethod::Vern7:  return {"vern7", 10, 7, 6, false, vern7_c, vern7_A, vern7_b, vern7_e, nullptr};
    case ODEMethod::Vern8:  return {"vern8", 13, 8, 7, false, vern8_c, vern8_A, vern8_b, vern8_e, nullptr};
    case ODEMethod::Vern9:  return {"vern9", 16, 9, 8, false, vern9_c, vern9_A, vern9_b, vern9_e, nullptr};
    case ODEMethod::RKV76:  return {"rkv76", 10, 7, 6, false, rkv76_c, rkv76_A, rkv76_b, rkv76_e, nullptr};
    case ODEMethod::DOP853: return {"dop853", 13, 8, 7, true, dop853_c, dop853_A, dop853_b, dop853_e, dop853_e3};
    default: throw std::runtime_error("SQuIDS: no native tableau for this ODEMethod");
  }
}

class NativeRK {
public:
  ///\brief Forget the step size carried over from the previous integration.
  void Reset(){ h_next=0; t_next=std::numeric_limits<double>::quiet_NaN(); }

  ///\brief Integrate y from t to t1 in place.
  ///\param f callable f(t, y, dydt) evaluating the right-hand side
  ///\param fixed_steps if nonzero, take this many equal steps without error control
  ///\note t is taken by value on purpose: SQuIDS::Derive overwrites the
  /// object's current time on every right-hand-side evaluation.
  template<typename RHS>
  void Integrate(const ExplicitRK& tab, RHS&& f, double t, double t1, double* y, size_t n,
                 double atol, double rtol, double hmin, double hmax, StepController ctl,
                 unsigned int fixed_steps, IntegrationStats& stats){
    if(t1==t)
      return;
    Prepare(tab, n);
    const double dir = (t1>t) ? 1.0 : -1.0;
    const unsigned int s = tab.stages;
    double* K1 = &K[0];

    f(t, y, K1);
    stats.rhs_evaluations++;

    if(fixed_steps>0){
      const double h = (t1-t)/fixed_steps;
      for(unsigned int i=0; i<fixed_steps; i++){
        Step(tab, f, t, h, y, n, stats);
        std::memcpy(y, &ynew[0], n*sizeof(double));
        t = (i+1==fixed_steps) ? t1 : t+h;
        stats.accepted_steps++;
        if(i+1<fixed_steps)
          AdvanceK1(tab, f, t, y, n, stats);
      }
      return;
    }

    const double abs_hmax = std::min(std::abs(hmax), std::abs(t1-t));
    const double abs_hmin = std::abs(hmin);
    double h;
    if(h_next!=0 && t==t_next && h_next*dir>0)
      h = h_next; // continue an integration that stopped exactly here
    else
      h = dir*InitialStep(tab, f, t, t1, y, n, atol, rtol, stats);
    h = dir*std::max(abs_hmin, std::min(std::abs(h), abs_hmax));

    const double k = tab.err_order+1.0;
    double err_old = 1e-4, h_old = h;
    bool rejected_last = false;
    for(;;){
      bool last = false;
      double h_unclipped = h;
      if(dir*(t+h-t1)>=0 || std::abs(t1-t)<=std::abs(h)*(1.0+1e-12)){
        h = t1-t;
        last = true;
      }

      Step(tab, f, t, h, y, n, stats);
      double err = ErrorNorm(tab, h, y, n, atol, rtol);

      if(err<=1.0){
        stats.accepted_steps++;
        double fac = Controller(ctl, err, err_old, h, h_old, k);
        if(rejected_last)
          fac = std::min(fac, 1.0);
        t = last ? t1 : t+h;
        std::memcpy(y, &ynew[0], n*sizeof(double));
        if(last){
          // the clipped final step says little about the natural step size
          h_next = (std::abs(h_unclipped) > std::abs(h)) ? h_unclipped : h*fac;
          t_next = t1;
          return;
        }
        AdvanceK1(tab, f, t, y, n, stats);
        err_old = std::max(err, 1e-4);
        h_old = h;
        rejected_last = false;
        h = dir*std::max(abs_hmin, std::min(std::abs(h*fac), std::abs(hmax)));
      }
      else{
        stats.rejected_steps++;
        double fac = std::isfinite(err) ? std::max(0.2, 0.9*std::pow(err, -1.0/k)) : 0.1;
        h *= fac;
        rejected_last = true;
        if(std::abs(h)<abs_hmin || t+h==t)
          throw std::runtime_error("SQuIDS::Evolve: native integrator step size underflow at t="
                                   +std::to_string(t));
      }
    }
  }

private:
  std::vector<double> K, ytmp, ynew;
  std::vector<std::vector<std::pair<unsigned int,double> > > rows;
  std::vector<std::pair<unsigned int,double> > bnz, enz, e3nz;
  const char* prepared_for = nullptr;
  double h_next = 0;
  double t_next = std::numeric_limits<double>::quiet_NaN();

  void Prepare(const ExplicitRK& tab, size_t n){
    const unsigned int s = tab.stages;
    if(K.size()!=s*n){
      K.assign(s*n, 0.0);
      ytmp.assign(n, 0.0);
      ynew.assign(n, 0.0);
    }
    if(prepared_for && std::strcmp(prepared_for, tab.name)==0)
      return;
    rows.assign(s, std::vector<std::pair<unsigned int,double> >());
    for(unsigned int i=0; i<s; i++)
      for(unsigned int j=0; j<i; j++)
        if(tab.A[i*s+j]!=0)
          rows[i].push_back(std::make_pair(j, tab.A[i*s+j]));
    bnz.clear(); enz.clear(); e3nz.clear();
    for(unsigned int j=0; j<s; j++){
      if(tab.b[j]!=0) bnz.push_back(std::make_pair(j, tab.b[j]));
      if(tab.e[j]!=0) enz.push_back(std::make_pair(j, tab.e[j]));
      if(tab.e3 && tab.e3[j]!=0) e3nz.push_back(std::make_pair(j, tab.e3[j]));
    }
    prepared_for = tab.name;
    Reset();
  }

  // out = y + h * sum_j a_j K_j over the listed nonzeros
  void Combine(const std::vector<std::pair<unsigned int,double> >& terms, double h,
               const double* y, double* out, size_t n){
    std::memcpy(out, y, n*sizeof(double));
    for(size_t q=0; q<terms.size(); q++){
      const double a = h*terms[q].second;
      const double* Kj = &K[terms[q].first*n];
      for(size_t i=0; i<n; i++)
        out[i] += a*Kj[i];
    }
  }

  // Computes stages 2..s and the new solution in ynew. K1 must hold f(t, y).
  template<typename RHS>
  void Step(const ExplicitRK& tab, RHS&& f, double t, double h, const double* y, size_t n,
            IntegrationStats& stats){
    const unsigned int s = tab.stages;
    for(unsigned int i=1; i<s; i++){
      const bool fsal_stage = tab.fsal && i==s-1;
      double* in = fsal_stage ? &ynew[0] : &ytmp[0];
      Combine(rows[i], h, y, in, n);
      f(t+tab.c[i]*h, in, &K[i*n]);
      stats.rhs_evaluations++;
    }
    if(!tab.fsal)
      Combine(bnz, h, y, &ynew[0], n);
  }

  // Makes K1 = f(t, y) for the step starting at the current (t, y).
  template<typename RHS>
  void AdvanceK1(const ExplicitRK& tab, RHS&& f, double t, const double* y, size_t n,
                 IntegrationStats& stats){
    if(tab.fsal){
      std::memcpy(&K[0], &K[(tab.stages-1)*n], n*sizeof(double));
    }
    else{
      f(t, y, &K[0]);
      stats.rhs_evaluations++;
    }
  }

  // Max norm of the error estimate, scaled by atol + rtol*max(|y_old|, |y_new|).
  double ErrorNorm(const ExplicitRK& tab, double h, const double* y, size_t n,
                   double atol, double rtol){
    double m = 0, m3 = 0;
    for(size_t i=0; i<n; i++){
      double e = 0;
      for(size_t q=0; q<enz.size(); q++)
        e += enz[q].second*K[enz[q].first*n+i];
      const double sc = atol+rtol*std::max(std::abs(y[i]), std::abs(ynew[i]));
      m = std::max(m, std::abs(h*e)/sc);
      if(tab.e3){
        double e3 = 0;
        for(size_t q=0; q<e3nz.size(); q++)
          e3 += e3nz[q].second*K[e3nz[q].first*n+i];
        m3 = std::max(m3, std::abs(h*e3)/sc);
      }
      if(!std::isfinite(m))
        return std::numeric_limits<double>::infinity();
    }
    if(tab.e3){
      // Hairer's DOP853 combination of the 5th and 3rd order estimates
      const double d = m*m+0.01*m3*m3;
      return d>0 ? m*m/std::sqrt(d) : 0.0;
    }
    return m;
  }

  static double Controller(StepController ctl, double err, double err_old, double h,
                           double h_old, double k){
    const double safety = 0.9;
    err = std::max(err, 1e-10);
    double fac;
    switch(ctl){
      case StepController::I:
        fac = safety*std::pow(err, -1.0/k);
        break;
      case StepController::H211b:
        fac = std::pow(safety/err, 0.25/k)*std::pow(safety/err_old, 0.25/k)
              *std::pow(std::abs(h/h_old), -0.25);
        break;
      case StepController::PI:
      default:
        fac = safety*std::pow(err, -0.7/k)*std::pow(err_old, 0.4/k);
        break;
    }
    return std::min(10.0, std::max(0.2, fac));
  }

  // Initial step size from the time scales tau1 = |y|/|y'| and tau2 = sqrt(|y|/|y''|),
  // h = tau * (0.01/|y|)^(1/(p+1)) in the norm scaled by the tolerances. Unlike the
  // textbook formula this is independent of the unit of time, and it still works
  // when f(t0, y0) == 0 (pure flavor states with no interactions start that way).
  template<typename RHS>
  double InitialStep(const ExplicitRK& tab, RHS&& f, double t, double t1, const double* y,
                     size_t n, double atol, double rtol, IntegrationStats& stats){
    const double span = std::abs(t1-t);
    const double dir = (t1>t) ? 1.0 : -1.0;
    const double* f0 = &K[0];
    double d0 = 0, d1 = 0;
    for(size_t i=0; i<n; i++){
      const double sc = atol+rtol*std::abs(y[i]);
      d0 += (y[i]/sc)*(y[i]/sc);
      d1 += (f0[i]/sc)*(f0[i]/sc);
    }
    d0 = std::sqrt(d0/n);
    d1 = std::sqrt(d1/n);
    const double tau1 = (d1>0 && d0>0) ? d0/d1 : std::numeric_limits<double>::infinity();
    double h0 = std::isfinite(tau1) ? 1e-3*tau1 : 1e-6*span;
    h0 = std::min(h0, span);

    double* f1 = &K[n]; // scratch: stage 2 is overwritten by the first real step
    for(size_t i=0; i<n; i++)
      ytmp[i] = y[i]+dir*h0*f0[i];
    f(t+dir*h0, &ytmp[0], f1);
    stats.rhs_evaluations++;
    double d2 = 0;
    for(size_t i=0; i<n; i++){
      const double sc = atol+rtol*std::abs(y[i]);
      d2 += ((f1[i]-f0[i])/sc)*((f1[i]-f0[i])/sc);
    }
    d2 = std::sqrt(d2/n)/h0;
    const double tau2 = (d2>0 && d0>0) ? std::sqrt(d0/d2) : std::numeric_limits<double>::infinity();

    const double tau = std::min(tau1, tau2);
    double h;
    if(!std::isfinite(tau))
      h = 1e-3*span; // f is constant: any step is exact, let the controller grow it
    else
      h = tau*std::pow(0.01/std::max(d0, 1.0), 1.0/(tab.order+1.0));
    return std::min(std::max(h, h0), span);
  }
};

} // namespace detail
} // namespace squids

#endif
