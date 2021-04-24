// CO2 solubility test by phreeqc code with EQUILIBRIUM_PHASES method
// Created by Jun Li, Wuhan Dimue company
// Created: 2018-4-9

#ifndef CO2SOLUINWATER_H_
#define CO2SOLUINWATER_H_

class CO2SoluInWater
{
public:
    CO2SoluInWater();
    ~CO2SoluInWater();

public:
    double calCO2SolubilityInPureWater(double temperatureK, double pressureBar, double volumeL);
    int calCO2SolubilityInPureWater(std::string experimentDataFileName);
private:

};




#endif