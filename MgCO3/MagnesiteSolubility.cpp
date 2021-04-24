#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"
#include "../EOS/EOSPhaseProp.h"
#include "../logError.h"

#include "MagnesiteSolubility.h"


MagnesiteSolubility::MagnesiteSolubility()
{
}

MagnesiteSolubility::~MagnesiteSolubility()
{
}

void MagnesiteSolubility::initialize()
{
    string inputFileName = "MagSolubility.phr";

    std::istream *db_cookie = NULL;
    std::istream *input_cookie = NULL;
    int errors;

    //    std::string inputFileName = "O2_solubility_water_iniSolution_purePhase_NoH2O.dat";
    std::string outFileName;
    std::string databaseFileName;

    phreeqc.ReSoCMode = 0;
    phreeqc.ReSoCMode_DEBUG = 1;

    errors = phreeqc.process_file_names(&db_cookie, &input_cookie, &inputFileName, &outFileName, &databaseFileName, TRUE);

    errors = phreeqc.do_initialize();

    phreeqc.Get_phrq_io()->push_istream(db_cookie);
    errors = phreeqc.read_database();
    phreeqc.Get_phrq_io()->clear_istream();
    //    std::cout << endl;
    phreeqc.Get_phrq_io()->push_istream(input_cookie);
    phreeqc.run_simulations();

    phreeqc.gasComponentDatabase.initialize();
}

void MagnesiteSolubility::calcMagnesiteSolubility(double temperatureK, double pressureBar, double mNaCl, double CO2moleNumber, double &magnesiiteMolality, double &Mg2Molality, double &pH)
{
    LDBLE atm_1bar = 0.9869233;
    phreeqc.gasComponentDatabase.initialize();
    conditionChange  update_timeStep;
    update_timeStep.pressure = pressureBar;
    update_timeStep.pressure *= atm_1bar;
    update_timeStep.temperature = temperatureK - 273.15;
    update_timeStep.masterChange["Na"] = mNaCl;
    update_timeStep.masterChange["Cl"] = mNaCl;
    update_timeStep.mass_of_H2O = 1.0;

    double volumeL = 60;

    update_timeStep.gasPhaseIn = TRUE;
    update_timeStep.gasPhaseInfo.volume = volumeL;
    update_timeStep.gasPhaseInfo.moleNumber["CO2"] = CO2moleNumber;
    update_timeStep.gasPhaseInfo.moleNumber["N2"] = 1.0 - CO2moleNumber;
    update_timeStep.mineralsUpdate["Magnesite"] = 1.0;


    phreeqc.run_simulation_timeStep(update_timeStep);

    gasPhase_results gas;
    waterSolution_results water;
    std::map<std::string, minerals> mineral;
    phreeqc.results_for_ReSoC(gas, water, mineral);

    map<string, double>::iterator it = water.masterSpecies_molality.find("Mg");
    if (it != water.masterSpecies_molality.end())
    {
        magnesiiteMolality = it->second;
    }

    it = water.species_molality.find("Mg+2");
    if (it != water.species_molality.end())
    {
        Mg2Molality = it->second;
    }

    pH = water.ph;
    //    update_timeStep.masterChange = water.masterSpecies_molality;
}

void MagnesiteSolubility::calcMagnesiteSolubilityNoCO2(double temperatureK, double pressureBar, double mNaCl, double &magMolality)
{
    LDBLE atm_1bar = 0.9869233;
    phreeqc.gasComponentDatabase.initialize();
    conditionChange  update_timeStep;
    update_timeStep.pressure = pressureBar;
    update_timeStep.pressure *= atm_1bar;
    update_timeStep.temperature = temperatureK - 273.15;
    update_timeStep.masterChange["Na"] = mNaCl;
    update_timeStep.masterChange["Cl"] = mNaCl;
    update_timeStep.mass_of_H2O = 1.0;

    double volumeL = 60;

    update_timeStep.gasPhaseIn = FALSE;
    /*update_timeStep.gasPhaseInfo.volume = volumeL;
    update_timeStep.gasPhaseInfo.moleNumber["CO2"] = 1.0;
    update_timeStep.mineralsUpdate["Magnesite"] = 1.0;*/

    update_timeStep.mineralsUpdate["Magnesite"] = 1.0;
    phreeqc.run_simulation_timeStep(update_timeStep);

    gasPhase_results gas;
    waterSolution_results water;
    std::map<std::string, minerals> mineral;
    phreeqc.results_for_ReSoC(gas, water, mineral);

    map<string, double>::iterator it = water.masterSpecies_molality.find("Mg");
    if (it != water.masterSpecies_molality.end())
    {
        magMolality = it->second;
    }


    //    update_timeStep.masterChange = water.masterSpecies_molality;
}
