//Global enum or structures.
//Created by Jun Li, Dimue company, Wuhan.
//Created: 2017-3-24
//Modified: 2017-7-18, BO(black-oil) is added in EOSChoice. struct SBOmodelInputData is added!

#ifndef ENUMDEF_H
#define ENUMDEF_H
#include <cmath>
#include <vector>
#include <string>
using namespace std;

const double R = 8.31446; //     J/mole/K
const double PI = 4.0*atan(1.0);
const double bigValue = 1.0e10;

//enum gasComponents
//{
//    H2O,
//    CO2,
//    CH4,
//    N2,
//    O2,
//    H2S,
//};

#define H2Og 0
#define CO2g 1
#define CH4g 2
#define N2g  3
#define O2g  4
#define H2Sg 5

enum EOSChoice
{
    PR76,
    PR78,
    RK,
    SRK,
    BO,
};

enum CalculationMethod
{
    Newton,
    SSI,
};

enum PhaseType
{
    GAS,
    OIL,
    WATER,
    TWOPHASE,
};

enum SystemPhase
{
    GASsys,
    OILsys,
    WATERsys,
    GASWATER,
    GASOIL,
    WATEROIL,
    GASWATEROIL,
};

struct SPhaseProperty
{
    PhaseType phase;
    double densitySI;
    double densityMole;
    double viscosity;
    vector<double> xMoleFraction;
    vector<double> dDensitySI_dm;
    double dDensitySI_dp;
    vector<double> dDensityMole_dm;
    double dDensityMole_dp;
    double dViscosityDp;
    vector<double> dViscosityDm;
    double saturation;
//    double zFactor;    //Actually, we don't need zFactor...
};

struct SBOPhaseProperty
{
    PhaseType phase;
    double density;
    double viscosity;
    vector<double> massFraction;
    double dDensityDp;
    double dViscosityDp;
    vector<double> dMassFractionDp;
    double saturation;
    double dSaturationDp;
    vector<double> dSaturationDm;
};

struct SInputDataBlock
{
    string dataBlockName;
    vector <string> lineData;
};

struct SzCompositionVsHight
{
    double hight;
    int componentNumber;
    vector<double> zComposition;
};

struct SsaltMolalityVsHight
{
    double hight;
    int saltNumber;
    vector<double> saltMolality;
};

struct ScomponentProperty
{
    int No;
    string name;
    double Mw;
    double Tcrit;
    double Pcrit;
    double zcrit;
    double acentricFactor;
    double shift;
};

struct SBOmodelInputData
{
    PhaseType phaseName;
    double density_standard;
 //   vector <vector <double> > table;
    vector<double> pressure;
    vector<double> B;
    vector<double> viscosity;
    vector<double> R;
};

struct phasePropertyForIni_compositional
{
    PhaseType phase;
    double densityMass;
    double densityMole;
    double dDensityMassDp;
    double dDensityMoleDp;
    vector<double> dDensityMassDm;
    vector<double> dDensityMoleDm;
	double viscosity;
    double dVisDp;
    vector<double> dVisDm; //to mole fraction for compositional model;
};

struct phasePropertyForIni_BO
{
    PhaseType phase;
    double densityMass;

    double dDensityMassDp;

    double viscosity;
    double dVisDp;
};

#endif

