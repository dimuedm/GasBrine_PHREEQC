#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"
#include "O2SolubilityInWater.h"

O2SolubilityInWater::O2SolubilityInWater()
{
}

O2SolubilityInWater::~O2SolubilityInWater()
{
}

double O2SolubilityInWater::calO2Solubility(double temperatureK, double pressureBar, double volumeL)
{
    Phreeqc phreeqc;
    std::istream *db_cookie = NULL;
    std::istream *input_cookie = NULL;
    int errors;

    std::string inputFileName = "O2_solubility_water_iniSolution_purePhase_NoH2O.dat";
    std::string outFileName;
    std::string databaseFileName;

    phreeqc.ReSoCMode = 0;
    phreeqc.ReSoCMode_DEBUG = 1;

    errors = phreeqc.process_file_names(&db_cookie, &input_cookie, &inputFileName, &outFileName, &databaseFileName, TRUE);

    errors = phreeqc.do_initialize();

    phreeqc.Get_phrq_io()->push_istream(db_cookie);
    errors = phreeqc.read_database();
    phreeqc.Get_phrq_io()->clear_istream();
    std::cout << endl;
    phreeqc.Get_phrq_io()->push_istream(input_cookie);
    phreeqc.run_simulations();

    LDBLE atm_1bar = 0.9869233;
    phreeqc.gasComponentDatabase.initialize();
    conditionChange  update_timeStep;
    update_timeStep.pressure = pressureBar;
    update_timeStep.pressure *= atm_1bar;
    update_timeStep.temperature = temperatureK - 273.15;
    update_timeStep.mass_of_H2O = 1.0;

    update_timeStep.gasPhaseIn = TRUE;
    update_timeStep.gasPhaseInfo.volume = volumeL;

    //DEFINE: we use molecular formula. NO "(g)" any more.
    update_timeStep.gasPhaseInfo.moleNumber["O2"] = 1;
    //    update_timeStep.gasPhaseInfo.moleFraction["H2O"] = 0.00;

    phreeqc.run_simulation_timeStep(update_timeStep);

    gasPhase_results gas;
    waterSolution_results water;
    std::map<std::string, minerals> mineral;
    phreeqc.results_for_ReSoC(gas, water, mineral);

    double O2Solubility = 0.0;
    std::map<std::string, double>::iterator it = water.masterSpecies_molality.find("O(0)");
    if (it != water.masterSpecies_molality.end())
    {
        O2Solubility = it->second*0.5;
    }

    cout << endl;
    cout << endl;
    cout << "O2 solubility in molality: " << O2Solubility << endl;
    cout << endl;
    cout << endl;

    return O2Solubility;
}

int O2SolubilityInWater::calO2Solubility()
{

    Phreeqc phreeqc;
    std::istream *db_cookie = NULL;
    std::istream *input_cookie = NULL;
    int errors;

    std::string inputFileName = "O2_solubility_water_iniSolution_purePhase_NoH2O.dat";
    std::string outFileName;
    std::string databaseFileName;

    phreeqc.ReSoCMode = 1;
    phreeqc.ReSoCMode_DEBUG = 0;

    errors = phreeqc.process_file_names(&db_cookie, &input_cookie, &inputFileName, &outFileName, &databaseFileName, TRUE);

    errors = phreeqc.do_initialize();

    phreeqc.Get_phrq_io()->push_istream(db_cookie);
    errors = phreeqc.read_database();
    phreeqc.Get_phrq_io()->clear_istream();
    std::cout << endl;
    phreeqc.Get_phrq_io()->push_istream(input_cookie);
    phreeqc.run_simulations();

    string inPutFile = "O2_solubility_in_water_data.txt";
    //   string inPutFile = experimentalDataFile;
    string outPutFile = "O2_in_Water_by_phreeqc.txt";
    ifstream inFile(inPutFile.c_str());
    ofstream oFile(outPutFile.c_str());
 //   int iLiterature;
    double temperature, pressure, mO2;
    if (!inFile.eof())
    {
        inFile.ignore(2000, '\n');
    }

    LDBLE atm_1bar = 0.9869233;
    int i = 0;
    phreeqc.gasComponentDatabase.initialize();

    gasPhase_results gas;
    waterSolution_results water;
    water.mass_of_H2O = 1.0;

    while (inFile >> pressure >> temperature >> mO2)
    {
        i++;
        conditionChange update_timeStep;
        update_timeStep.pressure = pressure*atm_1bar;
        update_timeStep.temperature = temperature - 273.15;
        update_timeStep.mass_of_H2O = 1.0;

        update_timeStep.gasPhaseIn = TRUE;
        update_timeStep.gasPhaseInfo.volume = 20.0;

        //DEFINE: we use molecular formula. NO "(g)" any more.
        update_timeStep.gasPhaseInfo.moleNumber["O2"] = 1;
        //        update_timeStep.gasPhaseInfo.moleFraction["H2O"] = 0.00;
        if (i > 1)
        {
            master *carbon = phreeqc.master_bsearch("O(0)");
            if (carbon != NULL)
            {
                update_timeStep.masterChange["O2"] = 0.0;
            }
        }

        phreeqc.run_simulation_timeStep(update_timeStep);

        gasPhase_results gas;
        waterSolution_results water;
        std::map<std::string, minerals> mineral;
        phreeqc.results_for_ReSoC(gas, water, mineral);

        double O2Solubility = 0.0;
        std::map<std::string, double>::iterator it = water.masterSpecies_molality.find("O(0)");
        if (it != water.masterSpecies_molality.end())
        {
            O2Solubility = it->second*0.5;
        }

        /*map<string, LDBLE>::iterator it;
        it = water.masterSpecies_molality.find("C(4)");*/
        LDBLE diff = (O2Solubility - mO2) / mO2;
        oFile << i << " " << temperature << "  " << pressure << "  " << mO2 << "  " << O2Solubility << "  " << diff << endl;
        cout << i << " "  << temperature << "  " << pressure << "  " << mO2 << "  " << O2Solubility << "  " << diff << endl;
    }

    return OK;
}