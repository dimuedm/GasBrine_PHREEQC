/*
For the use of phreeqc by ReSoC
Created by Li Jun, Dimue company 
Date: 2017-11-12
*/

#ifndef PHREEQC_RESOC_H
#define PHREEQC_RESOC_H


///*#define USE_LONG_DOUBLE*/
//#ifdef LDBLE_LONG_DOUBLE
//#define LDBLE long double
////#define SCANFORMAT "%Lf"
//#else
//#define LDBLE double
////#define SCANFORMAT "%lf"
////#endif
//#endif /* _INC_PHRQTYPE_H */

#include "EOS/EnumDef.h"
#include "EOS/EOSPhaseProp.h"
#include "EOS/GasComponentDataBase.h"

struct gasPhase_inputs
{
    long double volume; // Unit: L = 1.e-3 m3
    bool gasFugacityProvided = false;
    std::map<std::string, double> moleNumber;
    std::map<std::string, double> fugacity;
};

struct conditionChange
{
    double temperature;
    double pressure;
    std::map<std::string, double> masterChange;
    bool gasPhaseIn;
    double mass_of_H2O;
    gasPhase_inputs gasPhaseInfo;
    std::map<std::string, double> mineralsUpdate;
    float deltaTime;
    std::map<std::string, double> eqMineralChange;
    std::map<std::string, double> kineticMineralChange;
};


struct gasPhase_results
{
    double volume=0;
    double total_moles=0;
    std::map<std::string, double> moleNumber;
};

struct speciesComp
{
    string name;
    double molality_species;
    double molality_master;
    double activityCoeff;
};

struct waterSolution_results
{
    double mass_of_H2O;
    double activityOfH2O;
    double volume;
    double density;
//    LDBLE viscosity;
    double dDensity_dp;
//    LDBLE dViscosity_dp;
    double salinity; //Definiation: sum of (charge_of_cation*molality_of_cation)
    double ph;
    double pe;
    double specificConductance;
//    std::map<std::string, master*> masterSpecies;
    std::map<std::string, double> masterSpecies_molality;
    std::map<std::string, double> masterSpecies_moleNumber;
    std::map<std::string, double> species_molality;
    
};

struct minerals
{
    std::string formula;
    double initial_mole;
    double final_mole;
    double saturationIndex;
};

struct kineticPhases
{
    std::string formula;
    double initial_mole;
    double final_mole;
  //  double saturationIndex;
};







#endif
