
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
 *      Christopher Weaver (University of Wisconsin Madison)                   *
 *         chris.weaver@icecube.wisc.edu                                       *
 *      Jordi Salvado (University of Wisconsin Madison)                        *
 *         jsalvado@icecube.wisc.edu                                           *
 ******************************************************************************/


#ifndef SQUIDS_CONST_H
#define SQUIDS_CONST_H

#if __cplusplus < 201103L
#error C++11 compiler required. Update your compiler and use the flag -std=c++11
#endif

#include <memory>
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <stdexcept>

#include "version.h"

namespace squids{

///Contains physical and mathematical constants
class Const{
public :
  // mathematical constants //
  double pi;
  double piby2;
  double sqrt2;
  double ln2;
  // astronomical constants //
  double earthradius;
  double sunradius;
  ///// physics constants/////
  double GF;
  double Na;
  double sw_sq;
  double G;
  double alpha;
  double e_charge;
  /////////// units //////////
  // energy
  double TeV;
  double GeV;
  double MeV;
  double keV;
  double eV;
  double Joule;
  // mass
  double kg;
  double gr;
  // time
  double sec;
  double hour;
  double day;
  double year;
  // distance
  double meter;
  double cm;
  double km;
  double fermi;
  double angstrom;
  double AU;
  double ly;
  double parsec;
  // luminocity
  double picobarn;
  double femtobarn;
  // presure
  double Pascal;
  double atm;
  double psi;
  // temperature
  double Kelvin;
  // electromagnetic units
  double C;
  double A;
  double T;
  // angle
  double degree;

  int electron;
  int muon;
  int tau;
  int sterile1;
  int sterile2;
  int sterile3;

  double tau_lifetime;
  double tau_mass;
  double muon_lifetime;
  double muon_mass;
  double electron_mass;

  double proton_mass;
  double neutron_mass;

  ///\brief set the mixing angle between two states of the system
  ///
  ///\param state1 the (zero-based) index of the first state
  ///\param state2 the (zero-based) index of the second state;
  ///              must be larger than state1
  ///\param angle the angle to use
  void SetMixingAngle(unsigned int state1, unsigned int state2, double angle);

  ///\brief get the mixing angle between two states of the system
  ///
  ///\param state1 the (zero-based) index of the first state
  ///\param state2 the (zero-based) index of the second state;
  ///              must be larger than state1
  ///\return the mixing angle between the two states
  double GetMixingAngle(unsigned int state1, unsigned int state2) const;

  ///\brief set the energy splitting between two states of the system
  ///
  ///\param upperState the (zero-based) index of the upper state;
  ///                  the lower state is implicitly the ground state (0)
  ///\param diff the energy difference between the states
  void SetEnergyDifference(unsigned int upperState, double diff);

  ///\brief get the energy splitting between two states of the system
  ///
  ///\param upperState the (zero-based) index of the upper state;
  ///                  the lower state is implicitly the ground state (0)
  ///\return the energy difference between the states
  double GetEnergyDifference(unsigned int upperState) const;

  ///\brief set the complex phase between two states of the system
  ///
  ///\param state1 the (zero-based) index of the first state
  ///\param state2 the (zero-based) index of the second state;
  ///              must be larger than state1
  ///\param phase the phase to use
  void SetPhase(unsigned int state1, unsigned int state2, double phase);

  ///\brief get the complex phase between two states of the system
  ///
  ///\param state1 the (zero-based) index of the first state
  ///\param state2 the (zero-based) index of the second state;
  ///              must be larger than state1
  ///\return the phase between the two states
  double GetPhase(unsigned int state1, unsigned int state2) const;

  Const();
#if !(__PGI && __APPLE__ && __MACH__) //defaulted move constructors frighten pgi
  Const(Const&&)=default;
#endif
  ~Const();
#if !(__PGI && __APPLE__ && __MACH__) //defaulted move assignement frightens pgi
  Const& operator=(Const&&)=default;
#endif

  ///\brief get the complete matrix representing the transformation between bases
  std::unique_ptr<gsl_matrix_complex,void (*)(gsl_matrix_complex*)> GetTransformationMatrix(size_t) const;

private:
  // matrices
  // angles
  std::unique_ptr<gsl_matrix,void (*)(gsl_matrix*)> th;
  // cp-phases
  std::unique_ptr<gsl_matrix,void (*)(gsl_matrix*)> dcp;
  // energy differences
  std::unique_ptr<gsl_matrix,void (*)(gsl_matrix*)> de;
};

} //namespace squids

#endif //ifdef SQUIDS_CONST_H
