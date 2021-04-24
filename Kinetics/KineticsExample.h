#pragma once
#include "../MultiPhaseEquilibria/WaterMineralEquilibria.h"

class KineticsExample
{
public:
    KineticsExample();
    ~KineticsExample();

    void initialize();
    
    waterSolution_results water;
    void m_getResults();
    void m_stepRun(conditionChange &condition);
private:
    WaterMineralEquilibria m_equilibrium_;
    Phreeqc m_phreqc_;
};