
#include "../global_structures.h"
#include "../EOS/EOSPhaseProp.h"
#include "../logError.h"
#include "Lansfordite.h"
using namespace std;

Lansfordite::Lansfordite()
{
}

Lansfordite::~Lansfordite()
{
}

void Lansfordite::initialize()
{
    string inputFileName = "Lansfordite.phr";

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

//void Nesquehonite::updateLogK(string a_mineralName, double logK)
//{
//    phreeqc.updateLogK();
//}

double Lansfordite::searchLogK(double a_temperature, double a_pressure, double a_mNaCl, double a_moleFracCO2gas, double a_lansSolubility_measure, double logK_up, double logK_down)
{
    double logK = 0.5*(logK_down + logK_up);
    double pH = 7.0;
    double lansMolality, Mg2Molality;

    if (a_lansSolubility_measure <= 0)
    {
        logError logerror;
        logerror.LOGERROR("Error input of lansfordite solubilityŁĄ");
        return -999;
    }

    int iterations = 0;
    while (true)
    {
        iterations++;
        if (iterations > 100)
        {
            logError logerror;
            logerror.LOGERROR("Cannot converge!");
            return -999;
        }
        phreeqc.updateLogK("Lansfordite", logK);
        calcLansforditeSolubility(a_temperature, a_pressure, a_mNaCl, a_moleFracCO2gas, lansMolality, Mg2Molality, pH);
        //     cout << "nesMolality solubility = " << nesMolality << endl;
        if (abs((lansMolality - a_lansSolubility_measure) / a_lansSolubility_measure) < 1.0e-5)
        {
            return logK;
        }
        if (lansMolality > a_lansSolubility_measure)
        {
            logK_up = logK;

        }
        else
        {
            logK_down = logK;
        }
        logK = (logK_down + logK_up) / 2.0;
    }

}

void Lansfordite::calcLansforditeSolubility(double temperatureK, double pressureBar, double mNaCl, double CO2moleNumber, double &lansMolality, double &Mg2Molality, double &pH)
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
    update_timeStep.mineralsUpdate["Lansfordite"] = 10.0;


    phreeqc.run_simulation_timeStep(update_timeStep);

    gasPhase_results gas;
    waterSolution_results water;
    std::map<std::string, minerals> mineral;
    phreeqc.results_for_ReSoC(gas, water, mineral);

    map<string, double>::iterator it = water.masterSpecies_molality.find("Mg");
    if (it != water.masterSpecies_molality.end())
    {
        lansMolality = it->second;
    }

    it = water.species_molality.find("Mg+2");
    if (it != water.species_molality.end())
    {
        Mg2Molality = it->second;
    }

    pH = water.ph;
}

void Lansfordite::calcLansforditeSolubilityNoCO2(double temperatureK, double pressureBar, double mNaCl, double &lansMolality)
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

    update_timeStep.mineralsUpdate["Lansfordite"] = 1.0;
    phreeqc.run_simulation_timeStep(update_timeStep);

    gasPhase_results gas;
    waterSolution_results water;
    std::map<std::string, minerals> mineral;
    phreeqc.results_for_ReSoC(gas, water, mineral);

    map<string, double>::iterator it = water.masterSpecies_molality.find("Mg");
    if (it != water.masterSpecies_molality.end())
    {
        lansMolality = it->second;
    }

}