#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"
#include "CH4SolubilityInBrine.h"

CH4SolubilityInBrine::CH4SolubilityInBrine()
{
}

CH4SolubilityInBrine::~CH4SolubilityInBrine()
{
}

double CH4SolubilityInBrine::calCH4SolubilityInBrine(double temperatureK, double pressureBar, double mNaCl, double volumeL)
{
    Phreeqc phreeqc;
    std::istream *db_cookie = NULL;
    std::istream *input_cookie = NULL;
    int errors;

    std::string inputFileName = "CH4_solubility_in_brine.phr";
    std::string outFileName;
    std::string databaseFileName;

    phreeqc.ReSoCMode = 0;
    phreeqc.ReSoCMode_DEBUG = 1;

    errors = phreeqc.process_file_names(&db_cookie, &input_cookie, &inputFileName, &outFileName, &databaseFileName, TRUE);

    errors = phreeqc.do_initialize();

    phreeqc.Get_phrq_io()->push_istream(db_cookie);
    errors = phreeqc.read_database();
    phreeqc.Get_phrq_io()->clear_istream();
    std::cout << endl;
    phreeqc.Get_phrq_io()->push_istream(input_cookie);
    phreeqc.run_simulations();

    LDBLE atm_1bar = 0.9869233;
    phreeqc.gasComponentDatabase.initialize();
    conditionChange  update_timeStep;
    update_timeStep.pressure = pressureBar;
    update_timeStep.pressure *= atm_1bar;
    update_timeStep.temperature = temperatureK - 273.15;
    update_timeStep.masterChange["Na"] = mNaCl;
    update_timeStep.masterChange["Cl"] = mNaCl;
    update_timeStep.mass_of_H2O = 1.0;

    update_timeStep.gasPhaseIn = TRUE;
    update_timeStep.gasPhaseInfo.volume = volumeL;

    //DEFINE: we use molecular formula. NO "(g)" any more.
    update_timeStep.gasPhaseInfo.moleNumber["CH4"] = 1;

    phreeqc.run_simulation_timeStep(update_timeStep);

    gasPhase_results gas;
    waterSolution_results water;
    std::map<std::string, minerals> mineral;
    phreeqc.results_for_ReSoC(gas, water, mineral);

    double CH4Solubility = 0.0;
    std::map<std::string, double>::iterator it = water.masterSpecies_molality.find("C(-4)");
    if (it != water.masterSpecies_molality.end())
    {
        CH4Solubility = it->second;
    }

    cout << endl;
    cout << endl;
    cout << "CH4 solubility in molality: " << CH4Solubility << endl;
    cout << endl;
    cout << endl;

    return CH4Solubility;
}

int CH4SolubilityInBrine::calCH4SolubilityInBrine()
{
    return 0;
}

double CH4SolubilityInBrine::calCH4SolubilityInBrineNewMethod(double temperatureK, double pressureBar, double mNaCl)
{
    std::string inputFileName = "CH4_solubility_in_brine.phr";
    gasWater.initializePhreeqc(inputFileName);

    map<string, double> gasMole;
    gasMole["CH4"] = 10.0;
    map<string, double> speciesMolality;
    speciesMolality["Na"] = mNaCl;
    speciesMolality["Cl"] = mNaCl;
    gasWater.calPhaseEquilibria(gasMole, speciesMolality, temperatureK, pressureBar);
    double CO2solubility = gasWater.xAqueous[0] / gasWater.xAqueous[1] * 55.508;
    yH2O = gasWater.yGas[1];
    return CO2solubility;
}