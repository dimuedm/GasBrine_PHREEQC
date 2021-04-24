//Caclualte phase property.
//Created by Jun Li, Dimue company, Wuhan.
//Created: 2017-3-24
//Last modified: 2017-6-28

#ifndef EOSPHASEPROP_H
#define EOSPHASEPROP_H
#include <iostream>
#include <vector>
#include "Component.h"
#include "GasComponentDataBase.h"
//#include "PVTenvironment.h"
using namespace std;


class EOSPhaseProp
{
public:
    EOSPhaseProp();
    ~EOSPhaseProp();

    void initialize(vector<string> componentFormula, vector<double> zFeed, GasComponentDataBase *a_componentDatabase_, PhaseType phase, double temperature, double pressure, EOSChoice a_eosType);

    void calcFugacity(bool calculateDireviatives);
    void calcFugacity();
    void calDensity(bool derivative);

    void calViscosity(double densityMole);

    vector<double> zFeed;
    double A, B, a, b, c, dA_dp, dB_dp;
    double dA_dT, dB_dT, da_dT;
    vector<double> df_dp;
    vector<double> fugacity;
    vector<double> logPhi;
    vector<double> d2a_dxidT;
    vector< vector<double> > df_dx;
    vector<vector<double> > dlogPhi_dx;
    vector<double> da_dx, db_dx, dc_dx, dA_dx, dB_dx;
    vector<double> dlogPhi_dT;
    vector<double> df_dT;
    void calDT();
    double Z;
    double densityMole, densitySI;
    double dDensity_dP, dDensity_dT;
    double viscosity;
    double dViscosity_dP, dViscosity_dT;
    vector<double> dViscosity_dx;
    vector<double> dDensity_dx;
    double u, w;
    double moleVolume;
 
private:
    double      pressure;
    double      temperature;
    double      moleWeight;
    PhaseType   Phase;

//    vector<double> xComp;
    vector<Component>       components;
    map<string, Component>  *componentProperty;
//    GasComponentDataBase    *componentDatabase;

    //PVTenvironment *environment;
    int ComponentNumber;
    void calMoleWeight();
    void cal_abc();
    void calZ();
    EOSChoice eosType = PR78;

    double dZ_dp, dZ_dT;
    double wB, wBB;
    vector<double> dZ_dx;
    vector < vector<double> > kij;

    inline double FindMax(double x1, double x2, double x3);
    inline double FindMin(double x1, double x2, double x3);

    inline void calPseudoCriticalProperties();
    double Pc_pseudo, Tc_pseudo, Vc_pseudo, densityPseudoCritical;
    void calComponentViscosity();
    vector<double> viscosityComponent;
};


#endif

