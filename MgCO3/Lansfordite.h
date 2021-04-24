/*
This class calculates nesquehonite solubility at various temperature and pressure.
Created by Jun Li, Dimue, Wuhan.
Created: 2019-1-31
*/

#ifndef LANSFORDITE_H_
#define LANSFORDITE_H_
#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"

class Lansfordite
{
public:
    Lansfordite();
    ~Lansfordite();

public:
    void initialize();
    void calcLansforditeSolubility(double temperature, double pressure, double mNaCl, double CO2moleNumber, double &lansMolality, double &Mg2Molality, double &pH);
    void calcLansforditeSolubilityNoCO2(double temperature, double pressure, double mNaCl, double &magMolality);
    //    void updateLogK(string mineralName, double logK);
    double searchLogK(double a_temperature, double a_pressure, double a_mNaCl, double a_moleFracCO2gas, double a_lansSolubility_measure, double logK_up, double logK_down);
private:
    GasWaterEquilibria gasWaterEquili;
    Phreeqc phreeqc;
};


#endif