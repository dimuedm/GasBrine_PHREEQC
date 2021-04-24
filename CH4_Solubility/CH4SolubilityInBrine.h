/*

CH4 solubility in brine with modified pitzer database
Created by Jun Li, Dimue Technology Company
Created: 2018-4-11

*/

#ifndef CH4SOLUBILITYINBRINE_H_
#define CH4SOLUBILITYINBRINE_H_

#include "..\MultiPhaseEquilibria\GasIniWaterEquilibria.h"

class CH4SolubilityInBrine
{
public:
    CH4SolubilityInBrine();
    ~CH4SolubilityInBrine();

public:
    double calCH4SolubilityInBrine(double temperatureK, double pressureBar, double mNaCl, double volumeL);
    int calCH4SolubilityInBrine();

public:
    double calCH4SolubilityInBrineNewMethod(double temperatureK, double pressureBar, double mNaCl);
    double yH2O;

private:
    GasIniWaterEquilibria gasWater;
};





#endif