#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"
#include "CO2N2withPureWater.h"
#include "../EOS/EOSPhaseProp.h"
#include "../logError.h"

CO2N2withPureWater::CO2N2withPureWater()
{
}

CO2N2withPureWater::~CO2N2withPureWater()
{
}

void CO2N2withPureWater::calCO2N2H2Oeqilibria(double temperatureK, double pressureBar, map<string, double> &gasFeed, map<string, double> &masterSpecies)
{
 //   map<string, double> masterSpecies;
//    map<string, double> gasMoleNumber;

    double volume = 60.0; //m3

    int gasFeedSize = gasFeed.size();
    if (gasFeedSize <= 0)
    {
        logError logerror;
        logerror.LOGERROR("CO2N2withPureWater: Invalid input for gasFeed!");
        exit(-1);
    }
    vector<string> formula(gasFeedSize);
    vector<double> zFeed_gas(gasFeedSize);
    map<string, double>::iterator it = gasFeed.begin();

    double totalGas = 0.0;
    for (; it != gasFeed.end(); it++)
    {
        if (it->second < 0)
        {
            logError logerror;
            logerror.LOGERROR("CO2N2withPureWater: Invalid input for gasFeed!");
            exit(-1);
        }

        totalGas += it->second;
    }

    if (totalGas <= 0)
    {
        logError logerror;
        logerror.LOGERROR("CO2N2withPureWater: Invalid input for gasFeed!");
        exit(-1);
    }

    it = gasFeed.begin();
    int iGas = 0;
    for (; it != gasFeed.end(); it++)
    {
        it->second /= totalGas;
        formula[iGas] = it->first;
        zFeed_gas[iGas] = it->second;
        iGas++;
    }

    EOSPhaseProp phaseProp;
    GasComponentDataBase gasCompDataBase;
    gasCompDataBase.initialize();
    phaseProp.initialize(formula, zFeed_gas, &gasCompDataBase, GAS, temperatureK, pressureBar, PR78);
    phaseProp.calDensity(false);
    double density = phaseProp.densityMole;
    double gasMoleNumberTotal = volume * density;

    map<string, double> gasMoleNumbers;
    for (int i = 0; i < gasFeedSize; i++)
    {
        gasMoleNumbers[formula[i]] = zFeed_gas[i] * gasMoleNumberTotal;
    }
    gasWaterEquili.calPhaseEquilibria(gasMoleNumbers, masterSpecies, temperatureK, pressureBar);

    /*it = gasWaterEquili.finalWater.masterSpecies_molality.begin();
    for (; it != gasWaterEquili.finalWater.masterSpecies_molality.end(); it++)
    {

    }*/

    double N2Malality = gasWaterEquili.finalWater.masterSpecies_molality.find("N(0)")->second;
    N2Malality += gasWaterEquili.finalWater.masterSpecies_molality.find("N(3)")->second;
    N2Malality += gasWaterEquili.finalWater.masterSpecies_molality.find("N(-3)")->second;
    N2Malality += gasWaterEquili.finalWater.masterSpecies_molality.find("N(5)")->second;
    N2Malality /= 2.0;

    double CO2Molality = gasWaterEquili.finalWater.masterSpecies_molality.find("C(4)")->second;
    CO2Molality += gasWaterEquili.finalWater.masterSpecies_molality.find("C(-4)")->second;

    double H2OMolality = gasWaterEquili.finalWater.mass_of_H2O*1000.0 / 18.015;

    waterPhaseMoleFraction["N2"]  = N2Malality / (N2Malality + CO2Molality + H2OMolality);
    waterPhaseMoleFraction["CO2"] = CO2Molality / (N2Malality + CO2Molality + H2OMolality);
    waterPhaseMoleFraction["H2O"] = H2OMolality / (N2Malality + CO2Molality + H2OMolality);
}

void CO2N2withPureWater::initialize()
{
    string inputFile = "CO2N2withWater.phr";
    gasWaterEquili.initializePhreeqc(inputFile);
}