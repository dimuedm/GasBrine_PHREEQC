// Calculate O2 solubility in brine with phreeqc EQUILIBRIUM_PHASES method
// Created by Jun Li, Wuhan Dimue company
// Created: 2018-4-11

#ifndef O2SOLUBILITYINBRINE_H_
#define O2SOLUBILITYINBRINE_H_
#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"

class O2SolubilityInBrine
{
public:
    O2SolubilityInBrine();
    ~O2SolubilityInBrine();

public:
    double calO2SolubilityInBrine(double temperatureK, double pressureBar, double mNaCl, double volumeL);
    int calO2SolubilityInBrine();
    void initialize();
    double calO2Solubility_GasWaterEquili(double temperatureK, double pressureBar, double mNaCl);

private:
    GasWaterEquilibria gasWaterEquili;
};




#endif