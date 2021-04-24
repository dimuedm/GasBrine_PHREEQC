

#include "H2S_solubility.h"

H2S_solubility::H2S_solubility()
{
}

H2S_solubility::~H2S_solubility()
{
}

double H2S_solubility::calH2SSolubilityInBrineNewMethod(double temperatureK, double pressureBar, double mNaCl)
{
    std::string inputFileName = "H2S_solubility_Pitzer.phr";
    gasWater.initializePhreeqc(inputFileName);

    map<string, double> gasMole;
    gasMole["H2Sg"] = 10.0;
    map<string, double> speciesMolality;
    speciesMolality["Na"] = mNaCl;
    speciesMolality["Cl"] = mNaCl;
    gasWater.calPhaseEquilibria(gasMole, speciesMolality, temperatureK, pressureBar);
    double CO2solubility = gasWater.xAqueous[0] / gasWater.xAqueous[1] * 55.508;
    yH2O = gasWater.yGas[1];
    return CO2solubility;
}