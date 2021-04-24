/*
    Try to simulate dolomitization process with Mg contained fluid and calcite.
    Created by Jun Li, Wuhan Dimue company
    Created: 2018-5-22  
*/

#ifndef DOLOMITIZATION_H_
#define DOLOMITIZATION_H_

#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"
#include "../MultiPhaseEquilibria/ToolGasWaterEquilibria.h"
#include "../MultiPhaseEquilibria/WaterMineralEquilibria.h"


class Dolomitization
{
public:
    Dolomitization();
    ~Dolomitization();

    void initialize();
    void initializeCO2();
    void dolomitization(double temperatureK, double pressureBar, map<string, double> &mastrSpecies, map<string, double> &mineralComp);
    void dolomitizationCO2(double temperatureK, double pressureBar,map<string, double> &gasFeed, map<string, double> &mstrSpecies, map<string,double> &mineralComp);
    waterSolution_results waterFinal;
    gasPhase_results gasFinal;
    map<string, minerals> mineralFinal;
    double yH2O;

private:
    Phreeqc phreeqc;
    WaterMineralEquilibria waterMineral;
    ToolGasWaterEquilibria gasWaterMineral;
};




#endif