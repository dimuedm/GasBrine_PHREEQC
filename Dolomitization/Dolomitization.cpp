
#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"
#include "../EOS/EOSPhaseProp.h"
#include "../logError.h"
#include "Dolomitization.h"

Dolomitization::Dolomitization()
{
}

Dolomitization::~Dolomitization()
{
}

void Dolomitization::initialize()
{
    string inputFile = "Dolomization.phr";
    waterMineral.initialize(&phreeqc, inputFile);
}

void Dolomitization::initializeCO2()
{
    string inputFile = "Dolomization.phr";
    gasWaterMineral.initialize(&phreeqc, inputFile);

}

void Dolomitization::dolomitization(double temperatureK, double pressureBar, map<string, double> &mastrSpecies, map<string, double> &mineralComp)
{
    
    waterMineral.calPhaseEquilibria(&phreeqc, temperatureK, pressureBar, mastrSpecies, mineralComp, waterFinal, mineralFinal);
    //waterFinal = water;
}

void Dolomitization::dolomitizationCO2(double temperatureK, double pressureBar, map<string, double> &gasFeed, map<string, double> &mstrSpecies, map<string, double> &mineralComp)
{
    gasWaterMineral.calPhaseEquilibria(&phreeqc, gasFeed, mstrSpecies, mineralComp, temperatureK, pressureBar, gasFinal, waterFinal, mineralFinal, yH2O);
}