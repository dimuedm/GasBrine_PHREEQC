#include<iostream>
#include "Component.h"
using namespace std;

Component::Component()
{
}

Component::~Component()
{
}

// Peng-Robinson, 1976;
void Component::calcEOSProperty(double Temperature, EOSChoice eosType)
{
    double fw;

    switch (eosType)
    {
    case PR76:

        fw = 0.37464 + 1.54226*Omega - 0.26992*Omega*Omega;

        aCubic = 0.4572355*R*R*Tc*Tc / Pc;
        aCubic *= pow(1.0 + fw*(1.0 - sqrt(Temperature / Tc)), 2);
        bCubic = 0.0779607*R*Tc / Pc;
        daCubic_dT = -fw*(0.4572355*R*R*Tc*Tc / Pc)*(1 + fw*(1 - sqrt(Temperature / Tc))) / sqrt(Tc*Temperature);
        d2aCubic_dT2 = 0.0;
        
        break;
    case PR78:
        if (Omega < 0.5)
        {
            fw = 0.37464 + 1.54226*Omega - 0.26992*Omega*Omega;
        }
        else
        {
            fw = 0.3796 + 1.485*Omega - 0.1644*Omega*Omega + 0.01667*Omega*Omega*Omega;
        }
        aCubic = 0.4572355*R*R*Tc*Tc / Pc;
        aCubic *= pow(1.0 + fw*(1.0 - sqrt(Temperature / Tc)), 2);
        bCubic = 0.0779607*R*Tc / Pc;
        daCubic_dT = -fw*(0.4572355*R*R*Tc*Tc / Pc)*(1 + fw*(1 - sqrt(Temperature / Tc))) / sqrt(Tc*Temperature);
        d2aCubic_dT2 = 0.0;
        break;
    case RK:
        aCubic = 0.42748*R*R*pow(Tc, 2.5) / Pc / sqrt(Temperature);
        bCubic = 0.08664*R*Tc / Pc;
        daCubic_dT = 0.0;
        d2aCubic_dT2 = 0.0;
        break;
    case SRK:
        fw = 0.48 + 1.574*Omega - 0.176*Omega*Omega;
        aCubic = 0.42748*R*R*Tc*Tc / Pc;
        aCubic *= pow(1 + fw*(1 - sqrt(Temperature / Tc)), 2);
        bCubic = 0.08664*R*Tc / Pc;
        daCubic_dT = 0.0;
        d2aCubic_dT2 = 0.0;
        break;
    default:
        break;
    }
    aCubic_sqrt = sqrt(aCubic);
    cCubic = s_cb*bCubic;
}

