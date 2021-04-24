#pragma once
#include "..\MultiPhaseEquilibria\GasWaterEquilibria.h"
#include "..\MultiPhaseEquilibria\GasIniWaterEquilibria.h"

class H2S_solubility
{
public:
    H2S_solubility();
    ~H2S_solubility();

    void initialize();
    double calH2SSolubility(double temperatureK, double pressureBar, double mNaCl);
    double m_yH2O;
public:
    double calH2SSolubilityInBrineNewMethod(double temperatureK, double pressureBar, double mNaCl);
    double yH2O;

private:
    GasWaterEquilibria gasWaterEquili;
private:
    GasIniWaterEquilibria gasWater;


};

