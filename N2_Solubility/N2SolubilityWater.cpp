#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"

#include "N2SolubilityWater.h"


N2SolubilityWater::N2SolubilityWater()
{
}

N2SolubilityWater::~N2SolubilityWater()
{
}


double N2SolubilityWater::calN2Solubility(double temperatureK, double pressureBar)
{
//    GasWaterEquilibria gasWaterEquili;
    
    map<string, double> masterSpecies;
    map<string, double> gasMoleNumber;

    gasMoleNumber["N2"] = 1.0;
    
    gasWaterEquili.calPhaseEquilibria(gasMoleNumber, masterSpecies, temperatureK, pressureBar);
    
    map<string, double>::iterator it = gasWaterEquili.finalWater.masterSpecies_molality.begin();
    double nitrogenMolality = 0.0;
    for (; it != gasWaterEquili.finalWater.masterSpecies_molality.end(); it++)
    {
        string element = it->first;
        string nitrogen = "N(";
        if (element.find(nitrogen)!=string::npos)
        {
            nitrogenMolality += it->second;
        }
    }
    m_yH2O = gasWaterEquili.yH2O_final;
    return nitrogenMolality / 2.0;
}

void N2SolubilityWater::calN2Solubility()
{

}

void N2SolubilityWater::initialize()
{
    string inputFileName = "N2_solubility_Pitzer.phr";
    gasWaterEquili.initializePhreeqc(inputFileName);
}