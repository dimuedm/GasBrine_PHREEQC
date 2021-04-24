//Provide Component property
//Created by Jun Li, Dimue company, Wuhan.
//Created: 2017-3-24
//Last modified: 2017-6-28

#ifndef COMPONENT_H
#define COMPONENT_H
#include "EnumDef.h"
#include <string>
using namespace std;

class Component
{
public:
    
    Component();
    ~Component();
    string Name;
    
    double Pc;
    double Tc;
    double Vc;
    double MoleWeight;
    double Omega;
    double Zcrit;
    double s_cb;
    int number;
    int NO_forWaterSolu;
    double aCubic, bCubic, cCubic;
    void calcEOSProperty(double Temperature, EOSChoice eos);
    double daCubic_dT, d2aCubic_dT2, aCubic_sqrt;

private:
    

};



#endif