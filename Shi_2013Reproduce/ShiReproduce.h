/*

Try to reproduce experimental results from Shi et al., 2013
Created by Jun Li, Wuhan Dimue Company
Created: 2018-4-28

*/

#ifndef SHIREPRODUCE_H_
#define SHIREPRODUCE_H_
#include <vector>
#include "../Phreeqc_ReSoC.h"
#include "../Phreeqc.h"
#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"
#include "../MultiPhaseEquilibria/ToolGasWaterEquilibria.h"
class ShiReproduce
{
public:
    ShiReproduce();
    ~ShiReproduce();
    void initialize();
    void CO2saturation(double temperatureK, double pressureBar, waterSolution_results &water, gasPhase_results &gas, double &yH2O);
    void calciteDisolution(double temperatureK, double pressueBar, gasPhase_results &gas, waterSolution_results &water, map<string, minerals> &mineral);

private:
    Phreeqc phreeqc;
    ToolGasWaterEquilibria gasWaterEquilibria;
};

#endif