/*

    Try to reproduce experimental results of dolomite solubility from Pascale Be´ne´zeth, Ulf-Niklas Berninger, Nicolas Bovet c, Jacques Schott a,
Eric H. Oelkers, Geochimica et Cosmochimica Acta 224 (2018) 262–275.
    Created by Jun Li, Wuhan Dimue company
    Created: 2018-5-23

*/

#ifndef DOLOMITESOLUBILITY_H_
#define DOLOMITESOLUBILITY_H_

#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"
#include "../MultiPhaseEquilibria/ToolGasWaterEquilibria.h"
#include "../MultiPhaseEquilibria/WaterMineralEquilibria.h"

class DolomiteSolubility
{
public:
    DolomiteSolubility();
    ~DolomiteSolubility();

    void initialize();
    void initializeCO2();
    void calcDolomiteSolubility(double temperatureK, double pressureBar, double mNaCl, waterSolution_results &water);
    void calcDolomiteSolubilityExpReproduce(double temperatureK, double pressureBar, double mNaCl, waterSolution_results &water);
    void calcDolomiteSolubilityWithCO2(double temperatureK, double pressureBar, double mNaCl, map<string, double> &gasMoleNumber, waterSolution_results &water);
    waterSolution_results finalWater;
    gasPhase_results finalGas;
    double yH2O;
    map<string, minerals> finalMineral;

private:
    WaterMineralEquilibria waterDolomite;
    ToolGasWaterEquilibria gasWaterEquili;
    Phreeqc phreeqc;
};





#endif