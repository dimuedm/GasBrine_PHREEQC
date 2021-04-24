#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"
#include "ToolGasWaterEquilibria.h"
#include "../logError.h"

#include "WaterMineralEquilibria.h"

WaterMineralEquilibria::WaterMineralEquilibria()
{
}

WaterMineralEquilibria::~WaterMineralEquilibria()
{
}

void WaterMineralEquilibria::initialize(Phreeqc *phreeqc, string inputFileName)
{
    std::istream *db_cookie = NULL;
    std::istream *input_cookie = NULL;
    int errors;

    std::string outFileName;
    std::string databaseFileName;

    phreeqc->ReSoCMode = 1;
    phreeqc->ReSoCMode_DEBUG = 0;

    errors = phreeqc->process_file_names(&db_cookie, &input_cookie, &inputFileName, &outFileName, &databaseFileName, TRUE);

    errors = phreeqc->do_initialize();

    phreeqc->Get_phrq_io()->push_istream(db_cookie);
    errors = phreeqc->read_database();
    phreeqc->Get_phrq_io()->clear_istream();

    phreeqc->Get_phrq_io()->push_istream(input_cookie);

    phreeqc->run_simulations();

    phreeqc->gasComponentDatabase.initialize();
}

void WaterMineralEquilibria::calPhaseEquilibria(Phreeqc *phreeqc, double temperatureK, double pressureBar, map<string, double> &mastrSpecies, map<string, double> &mineralComp, waterSolution_results &finalWater, map<string, minerals> &finalMineral)
{
    const double moleWeightH2O = 18.105; //  g/mole

    LDBLE atm_1bar = 0.9869233;
    phreeqc->gasComponentDatabase.initialize();
    conditionChange  update_timeStep;
    update_timeStep.pressure = pressureBar;
    update_timeStep.pressure *= atm_1bar;
    update_timeStep.temperature = temperatureK - 273.15;
    update_timeStep.mass_of_H2O = 1.0;
    update_timeStep.masterChange = mastrSpecies;

    update_timeStep.gasPhaseIn = FALSE;
    update_timeStep.mineralsUpdate = mineralComp;
    phreeqc->run_simulation_timeStep(update_timeStep);

    gasPhase_results gas;
    waterSolution_results water;
    std::map<std::string, minerals> mineral;
    phreeqc->results_for_ReSoC(gas, water, mineral);
    finalWater = water;
    finalMineral = mineral;
}