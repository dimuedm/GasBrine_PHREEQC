/*

 N2 solubility test by phreeqc code with EQUILIBRIUM_PHASES method
 Gas-water equilibrium is calculated iteratively.
 Created by Jun Li, Wuhan Dimue company
 Created: 2018-4-25

*/

#ifndef N2SOLUBILITYWATER_H_
#define N2SOLUBILITYWATER_H_

#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"

class N2SolubilityWater
{
public:
    N2SolubilityWater();
    ~N2SolubilityWater();

    void initialize();
    double calN2Solubility(double temperatureK, double pressureBar);
    void calN2Solubility();

    double m_yH2O;
private:
   // GasWaterEquilibria gasWaterEquili;
    GasWaterEquilibria gasWaterEquili;
    
};


#endif