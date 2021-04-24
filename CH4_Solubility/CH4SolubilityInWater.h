/*

Calculate CH4 solubility in pure water by phreeqc.dat
Created by Jun Li, Wuhan Dimue Technology
Created: 2018-4-11

*/

#ifndef CH4SOLUBILITYINWATER_H_
#define CH4SOLUBILITYINWATER_H_

class CH4SolubilityInwater
{
public:
    CH4SolubilityInwater();
    ~CH4SolubilityInwater();

public:
    double calCH4SolubilityInWater(double temperatureK, double pressureBar, double volumeL);
    int calCH4SolubilityInWater();

private:

};





#endif