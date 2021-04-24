/*

Calculate calcite solubility, 
Created by Jun Li, Dimue technology company.
Created: 2018-5-22

*/

#ifndef CALCITESOLUBILITY_H_
#define CALCITESOLUBILITY_H_
#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"
#include "../MultiPhaseEquilibria/ToolGasWaterEquilibria.h"
#include "../MultiPhaseEquilibria/WaterMineralEquilibria.h"

class CalciteSolubility
{
public:
    CalciteSolubility();
    ~CalciteSolubility();

    void initialize();
    void calCalciteSolubility(double temperatureK, double pressureBar, map<string, double> gasFeed, map<string, double> masterSpecies, map<string, double> mineralComp, waterSolution_results &water);
    void calCalciteSolubility(double temperatureK, double pressureBar, map<string, double> masterSpecies, map<string, double> mineralComp, waterSolution_results &water);
    gasPhase_results gas_final;
    waterSolution_results water_final;
    map<string, minerals> mineral_final;
    double yH2O;

private:
    ToolGasWaterEquilibria gasWaterEquili;
    WaterMineralEquilibria waterMinEquili;
    Phreeqc phreeqc;
};




#endif