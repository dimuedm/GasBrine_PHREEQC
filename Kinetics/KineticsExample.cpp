#include "KineticsExample.h"

KineticsExample::KineticsExample()
{
}

KineticsExample::~KineticsExample()
{
}

void KineticsExample::initialize()
{
  //  string l_inputFile = "ex6_test1.phr";
   // string l_inputFile = "ex6C_test1.phr";
   // string l_inputFile = "Quartz_eq.phr";
    string l_inputFile = "GeoChemExample15.dat";
    m_equilibrium_.initialize(&m_phreqc_, l_inputFile);
}

void KineticsExample::m_getResults()
{
    gasPhase_results l_gas;
    map<string, minerals> l_mineral;
    m_phreqc_.results_for_ReSoC(l_gas, water, l_mineral);
    
}

void KineticsExample::m_stepRun(conditionChange &condition)
{
    m_phreqc_.run_simulation_timeStep(condition);
}