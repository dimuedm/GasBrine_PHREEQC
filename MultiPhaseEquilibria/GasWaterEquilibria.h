/*

This class handles fully gas-water equilbria with phreeqc functionality
For two component systems, the calculatin routine is:
(1)
(2)
(3)

For multi-component (more than 2, including H2O) system, the calculation routine is:
(1)
(2)
(3)

Created by Jun Li, Wuhan Dimue technology company.
Created: 2018-4-13


*/

#ifndef GASWATEREQUILIBRIA_H_
#define GASWATEREQUILIBRIA_H_

#include <vector>
#include "../Phreeqc_ReSoC.h"
#include "../Phreeqc.h"

using namespace std;


class GasWaterEquilibria
{
public:
    GasWaterEquilibria();
    ~GasWaterEquilibria();

public:
    void initializePhreeqc(string inputFileName);
    void calPhaseEquilibria(map<string, double> &gasMoleNumber, map<string, double> &masterSpecies, double temperatureK, double pressureBar);//, map<string, double> minerals);
    void calPhaseEquilibria2(map<string, double> &gasMoleNumber, map<string, double> &masterSpecies, double temperatureK, double pressureBar);
    void print_gasWaterEquilibria_screen();
    map<int, waterSolution_results> map_water;
    map<int, gasPhase_results> map_gas;
    map<int, map<string,minerals>> map_mineral;

    gasPhase_results finalGas;
    waterSolution_results finalWater;
    double yH2O_final;
private:
    Phreeqc                 phreeqc;
    vector <double>         KEquiliConst;
    vector<string>          gasCompositionName;
    int iterations = 0;

    double KEquiliH2O(double temperatureK, double pressureBar);

    double calGasWaterPhasePartition(vector<double> zFeed, double temperatureK, double pressureBar);
    double cal2ComponentGasWaterPartition(vector<double> zFeed, double temperatureK, double pressureBar);
    double calMultiComponentGasWaterPartition(vector<double> zFeed, double temperatureK, double pressureBar);

    int calKEquilConst(map<string,double> &moleFraction, double temperatureK, double pressureBar);

    bool calcGasMole(double &beta_GasMole, vector<double> &zFeed, vector<double> &Kvalue);
    double calRReqn(double &beta_gasMole, vector<double> &zFeed, vector<double> & Kvalue);
    double calRReqnPrim(double &beta, vector<double>&zFeed, vector<double>&Kvalue);
};






#endif