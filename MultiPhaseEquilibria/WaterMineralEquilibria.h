/*

   Calculate water-mineral equilibria (without gas present);
   Created by Jun Li, Wuhan Dimue company
   Created: 2018-5-23

*/

#ifndef WATERMINERALEQUILIBRIA_H_
#define WATERMINERALEQUILIBRIA_H_
#include <vector>
#include "../Phreeqc_ReSoC.h"
#include "../Phreeqc.h"

class WaterMineralEquilibria
{
public:
    WaterMineralEquilibria();
    ~WaterMineralEquilibria();
    void initialize(Phreeqc *phreeqc, string inputFileName);
    void calPhaseEquilibria(Phreeqc *phreeqc, double temperatureK, double pressureBar, map<string, double> &mastrSpecies, map<string, double> &mineralComp, waterSolution_results &finalWater, map<string, minerals> &finalMineral);

private:

};




#endif