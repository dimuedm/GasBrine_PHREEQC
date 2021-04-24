// CO2 solubility in NaCl solubitons, by phreeqc EQUILIBRIUM_PHASES method
// Created by Jun Li, Wuhan Dimue company
// Created: 2018-4-10

#ifndef CO2SOLUINBRINE_H_
#define CO2SOLUINBRINE_H_

#include "../MultiPhaseEquilibria/GasIniWaterEquilibria.h"
#include "../MultiPhaseEquilibria//GasWaterEquilibria.h"

using namespace std;

class CO2SoluInBrine
{
public:
    CO2SoluInBrine();
    ~CO2SoluInBrine();

public:
    double calCO2Solubility(double temperatureK, double pressureBar, double mNaCl, double volumeL);
    int calCO2Solubility();

    double calCO2SolubilityNewMethod(double temperatureK , double pressureBar, double mNaCl);
    double yH2O;
    GasWaterEquilibria * getEvaporite();

public:
    int CO2evaporationCal(double temperatureK, double pressureBar, double mNaCl);

private:
    GasIniWaterEquilibria gasWater;

private:
    GasWaterEquilibria evaporite;
};



#endif