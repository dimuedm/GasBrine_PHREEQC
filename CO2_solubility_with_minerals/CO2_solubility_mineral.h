/*

   Calculte CO2 solubility in water with saturated minerals: to evaluate the influcence of carbonates
   Created by Jun Li, Wuhan Dimue company
   Created: 2018-6-24

*/

#ifndef CO2_SOLUBILITY_MINREAL_
#define CO2_SOLUBILITY_MINREAL_

#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"
#include "../MultiPhaseEquilibria/ToolGasWaterEquilibria.h"
#include "../MultiPhaseEquilibria/WaterMineralEquilibria.h"

class CO2SolubilityWithMineral
{
public:
    CO2SolubilityWithMineral();
    ~CO2SolubilityWithMineral();

    double CO2Solu_with_minerals(double temperatureK, double pressureBar, map<string, double> &gasFeed, map<string, double> &mstrSpecies, map<string, double> &mineralComp);

private:

};



#endif