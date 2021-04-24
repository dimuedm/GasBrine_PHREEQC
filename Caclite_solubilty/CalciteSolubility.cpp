
#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"
#include "../EOS/EOSPhaseProp.h"
#include "../logError.h"

#include "CalciteSolubility.h"

CalciteSolubility::CalciteSolubility()
{
}

CalciteSolubility::~CalciteSolubility()
{
}

void  CalciteSolubility::initialize()
{
    string inputFileName = "ShiReproduce.phr";

    gasWaterEquili.initialize(&phreeqc, inputFileName);
}

void CalciteSolubility::calCalciteSolubility(double temperatureK, double pressureBar, map<string, double> gasFeed, map<string, double> masterSpecies, map<string, double> mineralComp, waterSolution_results &water)
{
    gasWaterEquili.calPhaseEquilibria(&phreeqc, gasFeed, masterSpecies, mineralComp, temperatureK, pressureBar, gas_final, water_final, mineral_final, yH2O);
    water = water_final;
}

void CalciteSolubility::calCalciteSolubility(double temperatureK, double pressureBar, map<string, double> masterSpecies, map<string, double> mineralComp, waterSolution_results &water)
{
    waterMinEquili.calPhaseEquilibria(&phreeqc, temperatureK, pressureBar, masterSpecies, mineralComp, water_final, mineral_final);
    water = water_final;
}