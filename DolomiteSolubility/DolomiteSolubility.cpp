
#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"
#include "../EOS/EOSPhaseProp.h"
#include "../logError.h"
#include "DolomiteSolubility.h"

DolomiteSolubility::DolomiteSolubility()
{
}

DolomiteSolubility::~DolomiteSolubility()
{
}

void DolomiteSolubility::initialize()
{
    string inputFileName = "ShiReproduce.phr";
    waterDolomite.initialize(&phreeqc, inputFileName);
}

void DolomiteSolubility::initializeCO2()
{
    string inputFileName = "ShiReproduce.phr";
    gasWaterEquili.initialize(&phreeqc, inputFileName);
}

void DolomiteSolubility::calcDolomiteSolubilityExpReproduce(double temperatureK, double pressureBar, double mNaCl, waterSolution_results &water)
{
    map<string, double> masterSpecies;
    masterSpecies["Na"] = 0.0975;
    masterSpecies["Cl"] = mNaCl;
    masterSpecies["H(1)"] = 2.5e-3;
    map<string, double> mineralComp;
    mineralComp["Dolomite"] = 0.96;
    mineralComp["Calcite"] = 0.0;
    mineralComp["Magnesite"] = 0.0;
    waterDolomite.calPhaseEquilibria(&phreeqc, temperatureK, pressureBar, masterSpecies, mineralComp, finalWater, finalMineral);
    water = finalWater;
}

void DolomiteSolubility::calcDolomiteSolubility(double temperatureK, double pressureBar, double mNaCl, waterSolution_results &water)
{
    map<string, double> masterSpecies;
    masterSpecies["Na"] = mNaCl;
    masterSpecies["Cl"] = mNaCl;
//    masterSpecies["H(1)"] = 2.5e-3;
    map<string, double> mineralComp;
    mineralComp["Dolomite"] = 1.0;
    mineralComp["Calcite"] = 0.0;
    mineralComp["Magnesite"] = 0.0;
    waterDolomite.calPhaseEquilibria(&phreeqc, temperatureK, pressureBar, masterSpecies, mineralComp, finalWater, finalMineral);
    water = finalWater;
}

void DolomiteSolubility::calcDolomiteSolubilityWithCO2(double temperatureK, double pressureBar, double mNaCl, map<string, double> &gasMoleNumber, waterSolution_results &water)
{
    map<string, double> masterSpecies;
    masterSpecies["Na"] = mNaCl;
    masterSpecies["Cl"] = mNaCl;
    //    masterSpecies["H(1)"] = 2.5e-3;
    map<string, double> mineralComp;
    mineralComp["Dolomite"] = 1.0;
    mineralComp["Calcite"] = 0.0;
    mineralComp["Magnesite"] = 0.0;
//    waterDolomite.calPhaseEquilibria(&phreeqc, temperatureK, pressureBar, masterSpecies, mineralComp, finalWater, finalMineral);
    gasWaterEquili.calPhaseEquilibria(&phreeqc, gasMoleNumber, masterSpecies, mineralComp, temperatureK, pressureBar, finalGas, finalWater, finalMineral, yH2O);
    water = finalWater;
}
