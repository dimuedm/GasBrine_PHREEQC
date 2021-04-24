#include "MgCO3Coexist.h"
MgCO3Coexit::MgCO3Coexit()
{
}

MgCO3Coexit::~MgCO3Coexit()
{
}

void MgCO3Coexit::initialize()
{
    string inputFileName = "MgCO3.phr";

    std::istream *db_cookie = NULL;
    std::istream *input_cookie = NULL;
    int errors;

    //    std::string inputFileName = "O2_solubility_water_iniSolution_purePhase_NoH2O.dat";
    std::string outFileName;
    std::string databaseFileName;

    phreeqc.ReSoCMode = 1;
    phreeqc.ReSoCMode_DEBUG = 0;

    errors = phreeqc.process_file_names(&db_cookie, &input_cookie, &inputFileName, &outFileName, &databaseFileName, TRUE);

    errors = phreeqc.do_initialize();

    phreeqc.Get_phrq_io()->push_istream(db_cookie);
    errors = phreeqc.read_database();
    phreeqc.Get_phrq_io()->clear_istream();
    //    std::cout << endl;
    phreeqc.Get_phrq_io()->push_istream(input_cookie);
    phreeqc.run_simulations();

    phreeqc.gasComponentDatabase.initialize();
}

void MgCO3Coexit::calcMgCO3Equilibria(double temperatureK, double pressureBar, double mNaCl, double CO2moleNumber, double &MgCO3Mineral, double &Mg2Molality, double &pH, vector<bool> &mineralIndex)
{
    LDBLE atm_1bar = 0.9869233;
    phreeqc.gasComponentDatabase.initialize();
    conditionChange  update_timeStep;
    update_timeStep.pressure = pressureBar;
    update_timeStep.pressure *= atm_1bar;
    update_timeStep.temperature = temperatureK - 273.15;
    update_timeStep.masterChange["Na"] = mNaCl;
    update_timeStep.masterChange["Cl"] = mNaCl;
    update_timeStep.mass_of_H2O = 1.0;

    double volumeL = 60;

    update_timeStep.gasPhaseIn = TRUE;
    update_timeStep.gasPhaseInfo.volume = volumeL;
    update_timeStep.gasPhaseInfo.moleNumber["CO2"] = CO2moleNumber;
    update_timeStep.gasPhaseInfo.moleNumber["N2"] = 1.0 - CO2moleNumber;
    update_timeStep.mineralsUpdate["Magnesite"] = 10.0;


    phreeqc.run_simulation_timeStep(update_timeStep);

    gasPhase_results gas;
  //  waterSolution_results water;
  //  std::map<std::string, minerals> mineral;
    phreeqc.results_for_ReSoC(gas, m_water, m_mineral);

    map<string, double>::iterator it = m_water.masterSpecies_molality.find("Mg");
    if (it != m_water.masterSpecies_molality.end())
    {
        MgCO3Mineral = it->second;
    }

    it = m_water.species_molality.find("Mg+2");
    if (it != m_water.species_molality.end())
    {
        Mg2Molality = it->second;
    }

    pH = m_water.ph;
    
    std::map<std::string, minerals>::iterator itMineral;
    mineralIndex.resize(3, false);

    itMineral = m_mineral.find("Magnesite");

    mineralIndex[0] = false;
    if (itMineral->second.final_mole > 0  && fabs(itMineral->second.saturationIndex)<1.e-6)
    {
        mineralIndex[0] = true;
    }

    itMineral = m_mineral.find("Nesquehonite");
    mineralIndex[1] = false;
    if (itMineral->second.final_mole > 0 && fabs(itMineral->second.saturationIndex)<1.e-6)
    {
        mineralIndex[1] = true;
    }

    itMineral = m_mineral.find("Lansfordite");
    mineralIndex[2] = false;
    if (itMineral->second.final_mole > 0 && fabs(itMineral->second.saturationIndex)<1.e-6)
    {
        mineralIndex[2] = true;
    }
}

void MgCO3Coexit::calcMgCO3EquilibriaNoCO2(double temperatureK, double pressureBar, double mNaCl, double &MgCO3Molality)
{
    LDBLE atm_1bar = 0.9869233;
    phreeqc.gasComponentDatabase.initialize();
    conditionChange  update_timeStep;
    update_timeStep.pressure = pressureBar;
    update_timeStep.pressure *= atm_1bar;
    update_timeStep.temperature = temperatureK - 273.15;
    update_timeStep.masterChange["Na"] = mNaCl;
    update_timeStep.masterChange["Cl"] = mNaCl;
    update_timeStep.mass_of_H2O = 1.0;

    double volumeL = 60;

    update_timeStep.gasPhaseIn = FALSE;
    /*update_timeStep.gasPhaseInfo.volume = volumeL;
    update_timeStep.gasPhaseInfo.moleNumber["CO2"] = 1.0;
    update_timeStep.mineralsUpdate["Magnesite"] = 1.0;*/

    update_timeStep.mineralsUpdate["Magnesite"] = 1.0;
    phreeqc.run_simulation_timeStep(update_timeStep);

    gasPhase_results gas;
    //waterSolution_results water;
    //std::map<std::string, minerals> mineral;
    phreeqc.results_for_ReSoC(gas, m_water, m_mineral);

    map<string, double>::iterator it = m_water.masterSpecies_molality.find("Mg");
    if (it != m_water.masterSpecies_molality.end())
    {
        MgCO3Molality = it->second;
    }
}