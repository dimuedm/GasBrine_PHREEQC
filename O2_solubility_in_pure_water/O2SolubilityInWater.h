// Calculate O2 solubility in pure water with phreeqc EQUILIBRIUM_PHASES method
// Created by Jun Li, Wuhan Dimue Company
// Created: 2018-4-10

#ifndef O2SOLUBILITYINWATER_H_
#define O2SOLUBILITYINWATER_H_

class O2SolubilityInWater
{
public:
    O2SolubilityInWater();
    ~O2SolubilityInWater();

public:
    double calO2Solubility(double temperatureK, double pressureBar, double volumeL);
    int calO2Solubility();

private:

};



#endif