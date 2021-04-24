#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"
#include "N2SolubilityInBrine.h"

N2SolubilityInBrine::N2SolubilityInBrine()
{
}

N2SolubilityInBrine::~N2SolubilityInBrine()
{
}

void N2SolubilityInBrine::initialize()
{
    string inputFileName = "N2_solubility_Pitzer.phr";
    gasWaterEquili.initializePhreeqc(inputFileName);
}

double N2SolubilityInBrine::calN2Solubility(double temperatureK, double pressureBar, double mNaCl)
{
    //    GasWaterEquilibria gasWaterEquili;

    map<string, double> masterSpecies;
    map<string, double> gasMoleNumber;

    gasMoleNumber["N2"] = 1.0;
    masterSpecies["Na"] = mNaCl;
    masterSpecies["Cl"] = mNaCl;

    gasWaterEquili.calPhaseEquilibria(gasMoleNumber, masterSpecies, temperatureK, pressureBar);

    map<string, double>::iterator it = gasWaterEquili.finalWater.masterSpecies_molality.begin();
    double nitrogenMolality = 0.0;
    for (; it != gasWaterEquili.finalWater.masterSpecies_molality.end(); it++)
    {
        string element = it->first;
        string nitrogen = "N(";
        if (element.find(nitrogen) != string::npos)
        {
            nitrogenMolality += it->second;
        }
    }
    m_yH2O = gasWaterEquili.yH2O_final;
    return nitrogenMolality / 2.0;
}

double N2SolubilityInBrine::calN2SolubilityInBrineNewMethod(double temperatureK, double pressureBar, double mNaCl)
{
    std::string inputFileName = "N2_solubility_Pitzer.phr";
    gasWater.initializePhreeqc(inputFileName);

    map<string, double> gasMole;
    gasMole["N2"] = 10.0;
    map<string, double> speciesMolality;
    speciesMolality["Na"] = mNaCl;
    speciesMolality["Cl"] = mNaCl;
    gasWater.calPhaseEquilibria(gasMole, speciesMolality, temperatureK, pressureBar);
    double N2solubility = gasWater.xAqueous[0] / gasWater.xAqueous[1] * 55.508;
    yH2O = gasWater.yGas[1];
    return N2solubility;
}