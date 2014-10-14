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

#include "const.h"

Const::Const(void){
    /* PHYSICS CONSTANTS
    #===========================================================================
    # NAME
    #===========================================================================
    */
    
    name = "STD";                    // Default values
    linestyle = "solid";             // Default linestyle in plots
    markerstyle = "*";               // Default marker style
    colorstyle = "red";              // Default color style
    savefilename = "output.dat";     // Default color style
    
    /*
    #===============================================================================
    # MATH
    #===============================================================================
    */
    pi=3.14159265358979;	    // Pi
    piby2=1.57079632679490;         // Pi/2
    sqrt2=1.41421356237310;         // Sqrt[2]
    ln2 = log(2.0);                 // log[2]
    
    /*
    #===============================================================================
    # EARTH 
    #===============================================================================
    */
    earthradius = 6371.0;	    // [km] Earth radius
    /*
    #===============================================================================
    # SUN 
    #===============================================================================
    */
    sunradius = 109.0*earthradius;  // [km] Sun radius 
    
    /*
    #===============================================================================
    # # PHYSICAL CONSTANTS
    #===============================================================================
    */
    GF = 1.16639e-23;	            // [eV^-2] Fermi Constant 
    Na = 6.0221415e+23;		    // [mol cm^-3] Avogadro Number
    sw_sq = 0.2312;                 // [dimensionless] sin(th_weinberg) ^2
    G  = 6.67300e-11;               // [m^3 kg^-1 s^-2]
    alpha = 1.0/137.0;              // [dimensionless] fine-structure constant 
    
    /*
    #===============================================================================
    # UNIT CONVERSION FACTORS
    #===============================================================================
    */
    // Energy
    TeV = 1.0e12;                   // [eV/TeV]
    GeV = 1.0e9;                    // [eV/GeV]
    MeV = 1.0e6;                    // [eV/MeV]
    keV = 1.0e3;                    // [eV/keV]
    eV  = 1.0;                      // [eV/eV]
    Joule = 1/1.60225e-19;          // [eV/J]
    // Mass
    kg = 5.62e35;                   // [eV/kg]
    gr = 1e-3*kg;                   // [eV/g] 
    // Time
    sec = 1.523e15;                 // [eV^-1/s]
    hour = 3600.0*sec;              // [eV^-1/h]
    day = 24.0*hour;                // [eV^-1/d]
    year = 365.0*day;               // [eV^-1/yr]
    // Distance
    meter = 5.076e6;                // [eV^-1/m]
    cm = 1.0e-2*meter;              // [eV^-1/cm]
    km = 1.0e3*meter;               // [eV^-1/km]
    fermi = 1.0e-15*meter;          // [eV^-1/fm]
    angstrom = 1.0e-10*meter;       // [eV^-1/A]
    AU = 149.60e9*meter;            // [eV^-1/AU]
    parsec = 3.08568025e16*meter;   // [eV^-1/parsec]
    // luminocity
    picobarn = 1.0e-36*pow(cm,2);       // [eV^-2/pb]
    femtobarn = 1.0e-39*pow(cm,2);      // [eV^-2/fb]
    // Presure
    Pascal = Joule/pow(meter,3);        // [eV^4/Pa]
    hPascal = 100.0*Pascal;         // [eV^4/hPa]
    atm = 101325.0*Pascal;          // [eV^4/atm]
    psi = 6893.0*Pascal;            // [eV^4/psi]
    // Temperature
    Kelvin = 1/1.1604505e4;         // [eV/K]
    // Angle
    degree = pi/180.0;              // [rad/degree]
    
    /*
    #===============================================================================
    # NEUTRINO OSCILLATION PARAMETERS 
    #===============================================================================
    */
    
    numneu = 3;                     // number of neutrinos
    numneumax = 6;                  // maximum neutrino number
    neutype = 0;                    // neutrino or antineutrino
    
    // angles
    th12 = 0.563942;
    th13 = 0.154085;
    th23 = piby2/2.0;
    th14 = 0.0;
    th24 = 0.0;
    th34 = 0.0;
    th15 = 0.0;
    th25 = 0.0;
    th35 = 0.0;
    th45 = 0.0;
    th16 = 0.0;
    th26 = 0.0;
    th36 = 0.0; 
    th46 = 0.0;
    th56 = 0.0;
    
    // square mass differences
    dm21sq = 7.65e-5;
    dm31sq = 2.47e-3;
    dm41sq = 0.0;
    dm51sq = 0.0;
    dm61sq = 0.0;
    
    // cp-phases
    delta1 = 0.0;
    delta2 = 0.0;
    delta3 = 0.0;
    
    // initializing matrices
    
    dmsq = gsl_matrix_alloc(numneumax,1);
    gsl_matrix_set(dmsq,0,0,0.0);
    gsl_matrix_set(dmsq,1,0,dm21sq);
    gsl_matrix_set(dmsq,2,0,dm31sq);
    gsl_matrix_set(dmsq,3,0,dm41sq);
    gsl_matrix_set(dmsq,4,0,dm51sq);
    gsl_matrix_set(dmsq,5,0,dm61sq);
    
    th = gsl_matrix_alloc(numneumax+1,numneumax+1);
    gsl_matrix_set(th,1,2,th12);
    gsl_matrix_set(th,1,2,th12);
    gsl_matrix_set(th,1,3,th13);
    gsl_matrix_set(th,2,3,th23);
    gsl_matrix_set(th,1,4,th14);
    gsl_matrix_set(th,2,4,th24);
    gsl_matrix_set(th,3,4,th34);
    gsl_matrix_set(th,1,5,th15);
    gsl_matrix_set(th,2,5,th25);
    gsl_matrix_set(th,3,5,th35);
    gsl_matrix_set(th,4,5,th45);
    gsl_matrix_set(th,1,6,th16);
    gsl_matrix_set(th,2,6,th26);
    gsl_matrix_set(th,3,6,th36);
    gsl_matrix_set(th,4,6,th46);
    gsl_matrix_set(th,5,6,th56);
    
    c = gsl_matrix_alloc(numneumax+1,numneumax+1);
    gsl_matrix_set(c,1,2,cos(th12));
    gsl_matrix_set(c,1,3,cos(th13));
    gsl_matrix_set(c,1,4,cos(th14));
    gsl_matrix_set(c,2,3,cos(th23));
    gsl_matrix_set(c,2,4,cos(th24));
    gsl_matrix_set(c,3,4,cos(th34));
    gsl_matrix_set(c,1,5,cos(th15));
    gsl_matrix_set(c,2,5,cos(th25));
    gsl_matrix_set(c,3,5,cos(th35));
    gsl_matrix_set(c,4,5,cos(th45));
    gsl_matrix_set(c,1,6,cos(th16));
    gsl_matrix_set(c,2,6,cos(th26));
    gsl_matrix_set(c,3,6,cos(th36));
    gsl_matrix_set(c,4,6,cos(th46));
    gsl_matrix_set(c,5,6,cos(th56));
    
    s = gsl_matrix_alloc(numneumax+1,numneumax+1);
    gsl_matrix_set(s,1,2,sin(th12));
    gsl_matrix_set(s,1,3,sin(th13));
    gsl_matrix_set(s,1,4,sin(th14));
    gsl_matrix_set(s,2,3,sin(th23));
    gsl_matrix_set(s,2,4,sin(th24));
    gsl_matrix_set(s,3,4,sin(th34));
    gsl_matrix_set(s,1,5,sin(th15));
    gsl_matrix_set(s,2,5,sin(th25));
    gsl_matrix_set(s,3,5,sin(th35));
    gsl_matrix_set(s,4,5,sin(th45));
    gsl_matrix_set(s,1,6,sin(th16));
    gsl_matrix_set(s,2,6,sin(th26));
    gsl_matrix_set(s,3,6,sin(th36));
    gsl_matrix_set(s,4,6,sin(th46));
    gsl_matrix_set(s,5,6,sin(th56));      
    
    dcp = gsl_matrix_alloc(numneumax-2+1,1);
    gsl_matrix_set(dcp,0,0,1.0);
    gsl_matrix_set(dcp,1,0,delta1);
    gsl_matrix_set(dcp,2,0,delta2);
    gsl_matrix_set(dcp,3,0,delta3);
    
    electron = 0;
    muon = 1;
    tau = 2;
    sterile1 = 3;
    sterile2 = 4;
    sterile3 = 5;
    
    tau_mass = 1776.82*MeV;
    tau_lifetime = 2.906e-13*sec;
    
    proton_mass = 938.272*MeV;
    neutron_mass = 939.565*MeV;
    
};


Const::~Const(void){
  gsl_matrix_free(dmsq);
  gsl_matrix_free(th);
  gsl_matrix_free(c);
  gsl_matrix_free(s);
  gsl_matrix_free(dcp);

}
int Const::Refresh(void){
    // reinitializing matrices
    
    gsl_matrix_set(dmsq,0,0,0.0);
    gsl_matrix_set(dmsq,1,0,dm21sq);
    gsl_matrix_set(dmsq,2,0,dm31sq);
    gsl_matrix_set(dmsq,3,0,dm41sq);
    gsl_matrix_set(dmsq,4,0,dm51sq);
    gsl_matrix_set(dmsq,5,0,dm61sq);
    
    gsl_matrix_set(th,1,2,th12);
    gsl_matrix_set(th,1,2,th12);
    gsl_matrix_set(th,1,3,th13);
    gsl_matrix_set(th,2,3,th23);
    gsl_matrix_set(th,1,4,th14);
    gsl_matrix_set(th,2,4,th24);
    gsl_matrix_set(th,3,4,th34);
    gsl_matrix_set(th,1,5,th15);
    gsl_matrix_set(th,2,5,th25);
    gsl_matrix_set(th,3,5,th35);
    gsl_matrix_set(th,4,5,th45);
    gsl_matrix_set(th,1,6,th16);
    gsl_matrix_set(th,2,6,th26);
    gsl_matrix_set(th,3,6,th36);
    gsl_matrix_set(th,4,6,th46);
    gsl_matrix_set(th,5,6,th56);
    
    gsl_matrix_set(c,1,2,cos(th12));
    gsl_matrix_set(c,1,3,cos(th13));
    gsl_matrix_set(c,1,4,cos(th14));
    gsl_matrix_set(c,2,3,cos(th23));
    gsl_matrix_set(c,2,4,cos(th24));
    gsl_matrix_set(c,3,4,cos(th34));
    gsl_matrix_set(c,1,5,cos(th15));
    gsl_matrix_set(c,2,5,cos(th25));
    gsl_matrix_set(c,3,5,cos(th35));
    gsl_matrix_set(c,4,5,cos(th45));
    gsl_matrix_set(c,1,6,cos(th16));
    gsl_matrix_set(c,2,6,cos(th26));
    gsl_matrix_set(c,3,6,cos(th36));
    gsl_matrix_set(c,4,6,cos(th46));
    gsl_matrix_set(c,5,6,cos(th56));
    
    gsl_matrix_set(s,1,2,sin(th12));
    gsl_matrix_set(s,1,3,sin(th13));
    gsl_matrix_set(s,1,4,sin(th14));
    gsl_matrix_set(s,2,3,sin(th23));
    gsl_matrix_set(s,2,4,sin(th24));
    gsl_matrix_set(s,3,4,sin(th34));
    gsl_matrix_set(s,1,5,sin(th15));
    gsl_matrix_set(s,2,5,sin(th25));
    gsl_matrix_set(s,3,5,sin(th35));
    gsl_matrix_set(s,4,5,sin(th45));
    gsl_matrix_set(s,1,6,sin(th16));
    gsl_matrix_set(s,2,6,sin(th26));
    gsl_matrix_set(s,3,6,sin(th36));
    gsl_matrix_set(s,4,6,sin(th46));
    gsl_matrix_set(s,5,6,sin(th56));      
    

    gsl_matrix_set(dcp,0,0,1.0);
    gsl_matrix_set(dcp,1,0,delta1);
    gsl_matrix_set(dcp,2,0,delta2);
    gsl_matrix_set(dcp,3,0,delta3);    
    
    return 0;
}

void Const::Set_th12(double val){
  th12 = val;
  Refresh();
}
void Const::Set_th13(double val){
  th13 = val;
  Refresh();
}
void Const::Set_th23(double val){
  th23 = val;
  Refresh();
}
void Const::Set_th14(double val){
  th14 = val;
  Refresh();
}
void Const::Set_th24(double val){
  th24 = val;
  Refresh();
}
void Const::Set_th34(double val){
  th34 = val;
  Refresh();
}
void Const::Set_th15(double val){
  th15 = val;
  Refresh();
}
void Const::Set_th25(double val){
  th25 = val;
  Refresh();
}
void Const::Set_th35(double val){
  th35 = val;
  Refresh();
}
void Const::Set_th45(double val){
  th45 = val;
  Refresh();
}
void Const::Set_th16(double val){
  th16 = val;
  Refresh();
}
void Const::Set_th26(double val){
  th26 = val;
  Refresh();
}
void Const::Set_th36(double val){
  th36 = val;
  Refresh();
}
void Const::Set_th46(double val){
  th46 = val;
  Refresh();
}
void Const::Set_th56(double val){
  th56 = val;
  Refresh();
}
void Const::Set_dm21sq(double val){
  dm21sq = val;
  Refresh();
}
void Const::Set_dm31sq(double val){
  dm31sq = val;
  Refresh();
}
void Const::Set_dm41sq(double val){
  dm41sq = val;
  Refresh();
}
void Const::Set_dm51sq(double val){
  dm51sq = val;
  Refresh();
}
void Const::Set_dm61sq(double val){
  dm61sq = val;
  Refresh();
}
void Const::Set_delta1(double val){
  delta1 = val;
  Refresh();
}
void Const::Set_delta2(double val){
  delta2 = val;
  Refresh();
}
void Const::Set_delta3(double val){
  delta3 = val;
  Refresh();
}
 


