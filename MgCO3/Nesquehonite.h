/*
This class calculates nesquehonite solubility at various temperature and pressure.
Created by Jun Li, Dimue, Wuhan.
Created: 2019-1-16
*/
#ifndef NESQUEHONITE_H_
#define NESQUEHONITE_H_

#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"
class Nesquehonite
{
public:
    Nesquehonite();
    ~Nesquehonite();

public:
    void initialize();
    void calcNesquehoniteSolubility(double temperature, double pressure, double mNaCl, double CO2moleNumber, double &nesMolality,  double &Mg2Molality, double &pH);
    void calcNesquehoniteSolubilityNoCO2(double temperature, double pressure, double mNaCl, double &magMolality);
//    void updateLogK(string mineralName, double logK);
    double searchLogK(double a_temperature, double a_pressure, double a_mNaCl, double a_moleFracCO2gas, double a_nesSolubility_measure, double logK_up, double logK_down);
private:
    GasWaterEquilibria gasWaterEquili;
    Phreeqc phreeqc;
};



#endif