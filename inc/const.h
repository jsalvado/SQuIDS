
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


#ifndef __CONST_H
#define __CONST_H

#include <string>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <iostream>
using namespace std;

class Const{
public : 
  // class identifiers
  string name; 
  string linestyle;
  string markerstyle;
  string colorstyle;
  string savefilename;
  // mathematical ants //
  double pi;
  double piby2; 
  double sqrt2;
  double ln2;
  // astronomical ants //
  double earthradius;
  double sunradius;
  ///// physics ants/////
  double GF;
  double Na;
  double sw_sq;
  double G;
  double alpha;
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
  double parsec;
  // luminocity
  double picobarn;
  double femtobarn;
  // presure
  double Pascal;
  double hPascal;
  double atm;
  double psi;
  // temperature
  double Kelvin;
  // angle
  double degree;
  ////// neutrino osc. param. //////////
  // basic
  int numneu; 
  int numneumax;
  int neutype;
  // mixing angles
  double th12;
  double th13;
  double th23;
  double th14;
  double th24;
  double th34;
  double th15;
  double th25;
  double th35;
  double th45;
  double th16;
  double th26;
  double th36;
  double th46;
  double th56;
  // square mass differences
  double dm21sq;
  double dm31sq;
  double dm41sq;
  double dm51sq;
  double dm61sq;
  // cp-phases
  double delta1;
  double delta2;
  double delta3;
  // matrices
  // angles
  gsl_matrix *th;
  // cosine mixing matrix
  gsl_matrix *c;
  // sine mixing matrix
  gsl_matrix *s;
  // cp-phases
  gsl_matrix *dcp;
  // square mass differences
  gsl_matrix *dmsq;
  // mixing matrix
  gsl_matrix_complex *U;
        
  int electron;
  int muon;
  int tau;
  int sterile1;
  int sterile2;
  int sterile3;
        
  double tau_lifetime;
  double tau_mass;
        
  double proton_mass;
  double neutron_mass;

  bool Set(string, double);
  int Refresh(void);
  Const(void);
  ~Const(void);
};

#endif
