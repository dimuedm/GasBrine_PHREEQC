/*

    This class calculates gas-water equilibrium as seperate tool.
    Created by Jun Li, Wuhan Dimue company.
    Created: 2018-4-28

*/

#ifndef TOOLGASWATEREQUILIBRIA_H_
#define TOOLGASWATEREQUILIBRIA_H_

#include <vector>
#include "../Phreeqc_ReSoC.h"
#include "../Phreeqc.h"

class ToolGasWaterEquilibria
{
public:
    ToolGasWaterEquilibria();
    ~ToolGasWaterEquilibria();

    void initialize(Phreeqc *a_phreeqc, string inputFileName);
    void calPhaseEquilibria(Phreeqc *a_phreeqc, map<string, double> &gasMoleNumber, map<string, double> &masterSpecies, map<string, double> &mineralComp, double temperatureK, double pressureBar, gasPhase_results &finalGas, waterSolution_results &finalWater, map<string, minerals> &finalMineral, double &yH2O);//, map<string, double> minerals);


private:
    //Phreeqc *m_phreeqc_ = NULL;
    double KEquiliH2O(double temperatureK, double pressureBar);
};


#endif