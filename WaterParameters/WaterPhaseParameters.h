// Gas and water parameters which will be used in water propery and gas-water mutual solubility moduel
// Created by Jun Li, Dimue Company, Wuhan
// Created date: 2017-5-15
// Modifed by: 2017-6-28
// Modified: 2018-1-10, add incompressible water model;

#ifndef WATERPARAMTERS_WATERPHASEPARAMETERS_H_
#define WATERPARAMTERS_WATERPHASEPARAMETERS_H_
//#include "PVTenvironment.h"
#include <vector>
using namespace std;

class WaterPhaseParameters
{
public:
    WaterPhaseParameters();
    ~WaterPhaseParameters();
    //void initialize(PVTenvironment *environment, double temperaure, double pressure);

    void getEquiliConst(double temperature, double pressure);
    void getActivityCoefficient(double molalityNaCl, double temperature, double pressure);
    void getDactivityCoefficient_dP(double molalityNaCl, double temperature, double pressure);

    double calActivityCoefficient_lamda(int NO_forWaterSolu, double temperature, double pressure);
    double calActivityCoefficient_zeta(int NO_forWaterSolu, double temperature, double pressure);
    
    double specificVolume(int NO_forWaterSolu, double temperature, double pressure, double molalityNaCl);
    vector<double> activityCoeff;
    vector<double> equilibriumConst;

    double calEquiliConst(int NO_forWaterSolu, double temperature, double pressure);
private:
//    PVTenvironment *environment;
    int iComponentNumber;
//    vector<Component>* components;

    double calWaterEquiliConst(double temperature, double pressure);
    double calCO2EquiliConst(double temperature, double pressure);
    double calCH4EquiliConst(double temperature, double pressure);
    double calN2EquiliConst(double temperature, double pressure);
    double calH2SEquiliConst(double temperature, double pressure);
    double calO2EquiliConst(double temperature, double pressure);

    double calCO2SpecificVolume(double temperature, double molalityNaCl);
    double calCH4SpecificVolume(double temperature, double pressure, double molalityNaCl);
    double calH2SSpecificVolume(double temperature, double pressure, double molalityNaCl);
    double calN2SpecificVolume(double temperature, double pressure, double molalityNaCl);

    const double R = 8.31446;
};



#endif // !WATERPHASEPARAMETERS_H
