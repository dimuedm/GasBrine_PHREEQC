/*

N2 solubility test by phreeqc code with EQUILIBRIUM_PHASES method
Gas-water equilibrium is calculated iteratively.
Created by Jun Li, Wuhan Dimue company
Created: 2018-4-26

*/

#ifndef N2SOLUBILITYINBRINE_H_
#define N2SOLUBILITYINBRINE_H_

#include "..\MultiPhaseEquilibria\GasWaterEquilibria.h"
#include "..\MultiPhaseEquilibria\GasIniWaterEquilibria.h"

class N2SolubilityInBrine
{
public:
    N2SolubilityInBrine();
    ~N2SolubilityInBrine();

    void initialize();
    double calN2Solubility(double temperatureK, double pressureBar, double mNaCl);
    void calN2Solubility();

    double m_yH2O;
public:
    double calN2SolubilityInBrineNewMethod(double temperatureK, double pressureBar, double mNaCl);
    double yH2O;
private:
    GasWaterEquilibria gasWaterEquili;
private:
    GasIniWaterEquilibria gasWater;
};


#endif