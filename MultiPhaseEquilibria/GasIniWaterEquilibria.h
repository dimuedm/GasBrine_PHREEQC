#ifndef GASINIWATEREQUILIBRIA_H_
#define GASINIWATEREQUILIBRIA_H_

#include <vector>
#include "../Phreeqc_ReSoC.h"
#include "../Phreeqc.h"
#include "../WaterParameters/WaterPhaseParameters.h"
using namespace std;

class GasIniWaterEquilibria
{
public:
    GasIniWaterEquilibria();
    ~GasIniWaterEquilibria();

public:
    void                    initializePhreeqc(string inputFileName);
    void                    updateReaction(conditionChange &newCondition);
    double                  calPhaseEquilibria(map<string, double> &gasMoleNumber, map<string, double> &masterSpecies, double temperatureK, double pressureBar);
    double                  calPhaseEquilibria_RealH2Omass(map<string, double>& gasMoleNumber, double &a_massOfH2OInWater, map<string, double>& masterSpecies, double temperatureK, double pressureBar, double &totalMoleNumber);
    gasPhase_results        finalGas;
    waterSolution_results   finalWater;
    map<string, minerals>   finalMineral;
    map<string, kineticPhases> finalKineticMineral;
    double                  yH2O_final;
    vector<double>          xAqueous;
    vector<double>          yGas;
    //map<string, double>     solubilityInMolality;
    //map<string, double>     yGasMap;
    map<string, int>        gasSequence;
    void                    generateGasCompName(map<string,double> gasMoleNumber);

    string                  getElemName(string a_gasCompName);
    double                  getElemCoeff(string a_gasCompName);
    int                     getGasCompIndex(string gasName);
    vector<string>          gasCompName;

private:
    Phreeqc                 phreeqc;
    WaterPhaseParameters    waterParam;
    
    int                     iterations = 0;
    vector<string>          elemName;
    vector<double>          elemCoeff;

    double                  KEquiliH2O(double temperatureK, double pressureBar);
    double                  waterSaturationPressure(double temperatureK);
    double                  calEquiliConst(int NO_forWaterSolu, double temperatureK, double pressureBar);
    double                  calWaterEquiliConst(double temperatureK, double pressureBar);
    double                  calCH4EquiliConst(double temperatureK, double pressureBar);
    double                  calCO2EquiliConst(double temperatureK, double pressureBar);
    double                  calH2SEquiliConst(double temperatureK, double pressureBar);
    double                  calN2EquiliConst(double temperatureK, double pressureBar);
    double                  calO2EquiliConst(double temperatureK, double pressureBar);

};




#endif