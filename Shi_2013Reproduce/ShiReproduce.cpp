#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"

#include "ShiReproduce.h"

using namespace std;

ShiReproduce::ShiReproduce()
{
}

ShiReproduce::~ShiReproduce()
{
}

void ShiReproduce::initialize()
{
    string inputFileName = "ShiReproduce.phr";

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


void ShiReproduce::CO2saturation(double temperatureK, double pressureBar, waterSolution_results &water, gasPhase_results &gas, double &yH2O)
{
    map<string, double> masterSpecies;
    map<string, double> gasMoleNumber;
    map<string, double> mineralComp;
    map<string, minerals> finalMineral;

    gasMoleNumber["CO2"] = 10.0;
    masterSpecies["Na"] = 0.1;
    masterSpecies["Cl"] = 0.1;

    /*double temperatureK = 273.15 + 25.0;
    double pressureBar = 1.0;

    waterSolution_results water;
    gasPhase_results gas;
    double yH2O;*/

    gasWaterEquilibria.calPhaseEquilibria(&phreeqc, gasMoleNumber, masterSpecies, mineralComp, temperatureK, pressureBar, gas, water, finalMineral, yH2O);

    map<string, double> ::iterator it = water.masterSpecies_molality.begin();
    for (; it != water.masterSpecies_molality.end(); it++)
    {
        cout << it->first << "  " << it->second << endl;
    }
    cout << "yH2O = " << yH2O << endl;
}

void ShiReproduce::calciteDisolution(double temperatureK, double pressureBar, gasPhase_results &gas, waterSolution_results &water, map<string, minerals> &mineral)
{
    /*map<string, double> ::iterator it = water.masterSpecies_molality.begin();
    for (; it != water.masterSpecies_molality.end(); it++)
    {

    }*/
    LDBLE atm_1bar = 0.9869233;
    conditionChange updateCondition;
    updateCondition.gasPhaseIn = false;
    updateCondition.masterChange = water.masterSpecies_molality;
    updateCondition.pressure = pressureBar*atm_1bar;
    updateCondition.temperature = temperatureK - 273.15;
    updateCondition.mass_of_H2O = water.mass_of_H2O;

    updateCondition.mineralsUpdate["Calcite"] = 1.0;
    phreeqc.run_simulation_timeStep(updateCondition);
    phreeqc.results_for_ReSoC(gas, water, mineral);

    cout << endl;
    map<string, double> ::iterator it = water.masterSpecies_molality.begin();
    for (; it != water.masterSpecies_molality.end(); it++)
    {
        cout << it->first << "  " << it->second << endl;
    }
}