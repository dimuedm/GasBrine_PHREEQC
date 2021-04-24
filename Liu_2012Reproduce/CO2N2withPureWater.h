/*

    This class tries to reproduce the experimental results from Liu et al., 2012, dx.doi.org/10.1021/je3000958 | J. Chem. Eng. Data 2012, 57, 1928−1932
    Created by Jun Li, Wuhan Dimue company. 
    Created: 2018-4-30.

*/

#ifndef CO2N2WITHPUREWATER_H_
#define CO2N2WITHPUREWATER_H_
#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"

class CO2N2withPureWater
{
public:
    CO2N2withPureWater();
    ~CO2N2withPureWater();

public:
    void calCO2N2H2Oeqilibria(double temperatureK, double pressureBar, map<string, double> &gasFeed, map<string, double> &masterSpecies);
    map<string, double> waterPhaseMoleFraction;
    map<string, double> gasMoleFraction;
    void initialize();

private:
    GasWaterEquilibria gasWaterEquili;
};




#endif