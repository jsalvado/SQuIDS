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

#include <cmath>
#include "SQuIDS/const.h"
#include "SQuIDS/SU_inc/dimension.h"

Const::Const(){
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
    // luminosity
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
    
    // initializing matrices
    
    dmsq = gsl_matrix_alloc(SQUIDS_MAX_HILBERT_DIM-1,1);
    for(unsigned int i=0; i<SQUIDS_MAX_HILBERT_DIM-1; i++)
        gsl_matrix_set(dmsq,i,0,0.0);
    
    th = gsl_matrix_alloc(SQUIDS_MAX_HILBERT_DIM,SQUIDS_MAX_HILBERT_DIM);
    for(unsigned int i=0; i<SQUIDS_MAX_HILBERT_DIM; i++){
        for(unsigned int j=0; j<SQUIDS_MAX_HILBERT_DIM; j++)
            gsl_matrix_set(th,i,j,0.0);
    }
    
    dcp = gsl_matrix_alloc(SQUIDS_MAX_HILBERT_DIM,SQUIDS_MAX_HILBERT_DIM);
    for(unsigned int i=0; i<SQUIDS_MAX_HILBERT_DIM; i++){
        for(unsigned int j=0; j<SQUIDS_MAX_HILBERT_DIM; j++)
            gsl_matrix_set(dcp,i,j,0.0);
    }
    
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

Const::~Const(){
    gsl_matrix_free(dmsq);
    gsl_matrix_free(th);
    gsl_matrix_free(dcp);
}

void Const::SetMixingAngle(unsigned int state1, unsigned int state2, double angle){
    if(state2<=state1)
        throw std::runtime_error("Const::SetMixingAngle: state indices should be ordered and unequal"
                                 " (Got "+std::to_string(state1)+" and "+std::to_string(state2)+")");
    if(state1>=SQUIDS_MAX_HILBERT_DIM-1)
        throw std::runtime_error("Const::SetMixingAngle: First state index must be less than " SQUIDS_MAX_HILBERT_DIM_STR-1);
    if(state2>=SQUIDS_MAX_HILBERT_DIM)
        throw std::runtime_error("Const::SetMixingAngle: Second mass state index must be less than " SQUIDS_MAX_HILBERT_DIM_STR);
    
    gsl_matrix_set(th,state1,state2,angle);
}

double Const::GetMixingAngle(unsigned int state1, unsigned int state2) const{
    if(state2<=state1)
        throw std::runtime_error("Const::GetMixingAngle: state indices should be ordered and unequal"
                                 " (Got "+std::to_string(state1)+" and "+std::to_string(state2)+")");
    if(state1>=SQUIDS_MAX_HILBERT_DIM-1)
        throw std::runtime_error("Const::GetMixingAngle: First state index must be less than " SQUIDS_MAX_HILBERT_DIM_STR-1);
    if(state2>=SQUIDS_MAX_HILBERT_DIM)
        throw std::runtime_error("Const::GetMixingAngle: Second mass state index must be less than " SQUIDS_MAX_HILBERT_DIM_STR);
    
    return(gsl_matrix_get(th,state1,state2));
}

void Const::SetSquaredEnergyDifference(unsigned int upperState, double sqdiff){
    if(upperState==0)
        throw std::runtime_error("Const::SetSquaredEnergyDifference: Upper state index must be greater than 0");
    if(upperState>=SQUIDS_MAX_HILBERT_DIM)
        throw std::runtime_error("Const::SetSquaredEnergyDifference: Upper state index must be less than " SQUIDS_MAX_HILBERT_DIM_STR);
    
    gsl_matrix_set(dmsq,upperState-1,0,sqdiff);
}

double Const::GetSquaredEnergyDifference(unsigned int upperState) const{
    if(upperState==0)
        throw std::runtime_error("Const::GetSquaredEnergyDifference: Upper state index must be greater than 0");
    if(upperState>=SQUIDS_MAX_HILBERT_DIM)
        throw std::runtime_error("Const::GetSquaredEnergyDifference: Upper state index must be less than " SQUIDS_MAX_HILBERT_DIM_STR);
    
    return(gsl_matrix_get(dmsq,upperState-1,0));
}

void Const::SetPhase(unsigned int state1, unsigned int state2, double phase){
    if(state2<=state1)
        throw std::runtime_error("Const::SetPhase: State indices should be ordered and unequal"
                                 " (Got "+std::to_string(state1)+" and "+std::to_string(state2)+")");
    if(state2>=SQUIDS_MAX_HILBERT_DIM)
        throw std::runtime_error("Const::SetPhase: Upper state index must be less than " SQUIDS_MAX_HILBERT_DIM_STR);
    
    gsl_matrix_set(dcp,state1,state2,phase);
}

double Const::GetPhase(unsigned int state1, unsigned int state2) const{
    if(state2<=state1)
        throw std::runtime_error("Const::GetPhase: State indices should be ordered and unequal"
                                 " (Got "+std::to_string(state1)+" and "+std::to_string(state2)+")");
    if(state2>=SQUIDS_MAX_HILBERT_DIM)
        throw std::runtime_error("Const::GetPhase: Upper state index must be less than " SQUIDS_MAX_HILBERT_DIM_STR);
    
    return(gsl_matrix_get(dcp,state1,state2));
}
