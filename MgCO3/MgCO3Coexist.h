#pragma once

#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"

class MgCO3Coexit
{
public:
    MgCO3Coexit();
    ~MgCO3Coexit();

public:
    void initialize();
    void calcMgCO3Equilibria(double temperature, double pressure, double mNaCl, double CO2moleNumber, double &nesMolality, double &Mg2Molality, double &pH, vector<bool> &mineralIndex);
    void calcMgCO3EquilibriaNoCO2(double temperature, double pressure, double mNaCl, double &magMolality);

    waterSolution_results m_water;
    map<string, minerals> m_mineral;

private:
    GasWaterEquilibria gasWaterEquili;
    Phreeqc phreeqc;

};

