#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"
#include "ToolGasWaterEquilibria.h"
#include "../logError.h"

using namespace std;

ToolGasWaterEquilibria::ToolGasWaterEquilibria()
{
}

ToolGasWaterEquilibria::~ToolGasWaterEquilibria()
{
}

void ToolGasWaterEquilibria::initialize(Phreeqc *a_phreeqc, string inputFileName)
{
    Phreeqc *phreeqc = a_phreeqc;
    std::istream *db_cookie = NULL;
    std::istream *input_cookie = NULL;
    int errors;

    //    std::string inputFileName = "O2_solubility_water_iniSolution_purePhase_NoH2O.dat";
    std::string outFileName;
    std::string databaseFileName;

    phreeqc->ReSoCMode = 1;
    phreeqc->ReSoCMode_DEBUG = 0;

    errors = phreeqc->process_file_names(&db_cookie, &input_cookie, &inputFileName, &outFileName, &databaseFileName, TRUE);

    errors = phreeqc->do_initialize();

    phreeqc->Get_phrq_io()->push_istream(db_cookie);
    errors = phreeqc->read_database();
    phreeqc->Get_phrq_io()->clear_istream();
    //    std::cout << endl;
    phreeqc->Get_phrq_io()->push_istream(input_cookie);
    phreeqc->run_simulations();

    phreeqc->gasComponentDatabase.initialize();
}

void ToolGasWaterEquilibria::calPhaseEquilibria(Phreeqc *a_phreeqc, map<string, double> &gasPhaseMoleNumber, map<string, double> &masterSpecies, map<string, double> &mineralComp, double temperatureK, double pressureBar, gasPhase_results &finalGas, waterSolution_results &finalWater, map<string, minerals> &finalMineral, double &yH2O_final)
{
    Phreeqc *phreeqc = a_phreeqc;
    vector<string> gasCompositionName;

    int sizeOfGasPhaseComponents = (int)gasPhaseMoleNumber.size();
    vector<double> zFeed(sizeOfGasPhaseComponents);
    gasCompositionName.resize(sizeOfGasPhaseComponents);
    const double moleWeightH2O = 18.105; //  g/mole
    double totalMole = 0.0;
    bool H2OFound = false;
    string nameH2O = "";
    double moleH2O = 0.0;
    double KH2O = KEquiliH2O(temperatureK, pressureBar);
    map<string, double>::iterator it = gasPhaseMoleNumber.begin();
    int iComponent = 0;
    for (; it != gasPhaseMoleNumber.end(); it++)
    {
        if (it->second<0)
        {
            logError logerror;
            logerror.LOGERROR("GasWaterEquilibria: Negative mole fraction!");
            exit(-1);
        }
        if ((it->first == "H2O") || (it->first == "h2o") || (it->first == "H2o") || (it->first == "h2O"))
        {
            if (H2OFound)
            {
                logError logerror;
                logerror.LOGERROR("GasWaterEquilibria: More than 1 H2O component information are found!");
                exit(-1);
            }
            nameH2O = it->first;
            moleH2O = it->second;
            H2OFound = true;
        }

        totalMole += it->second;
    }


    it = gasPhaseMoleNumber.begin();
    for (; it != gasPhaseMoleNumber.end(); it++)
    {
//        it->second /= totalMole;
        if (it->first != nameH2O)
        {
            zFeed[iComponent] = it->second/totalMole;
            gasCompositionName[iComponent] = it->first;
            iComponent++;
        }
        else
        {
            zFeed[sizeOfGasPhaseComponents - 1] = it->second/totalMole;
            gasCompositionName[sizeOfGasPhaseComponents - 1] = nameH2O;
        }
    }

    if (!H2OFound)
    {
        zFeed.push_back(0.0);
        nameH2O = "H2O";
        gasCompositionName.push_back("H2O");
        gasPhaseMoleNumber[nameH2O] = 0.0;
        sizeOfGasPhaseComponents++;
        moleH2O = 0.0;
    }

    LDBLE atm_1bar = 0.9869233;
    phreeqc->gasComponentDatabase.initialize();
    conditionChange  update_timeStep;
    update_timeStep.pressure = pressureBar;
    update_timeStep.pressure *= atm_1bar;
    update_timeStep.temperature = temperatureK - 273.15;
    update_timeStep.mass_of_H2O = 1.0;

    update_timeStep.gasPhaseIn = TRUE;
    map <string, double> gasPhaseMoleNumber_update;
 //   update_timeStep.gasPhaseInfo.volume = 60.0;

    gasPhaseMoleNumber_update = gasPhaseMoleNumber;
    gasPhaseMoleNumber_update.erase(nameH2O);
    double moleNumberOfH2OinGas_old = moleH2O;
    update_timeStep.masterChange = masterSpecies;

    int iterations = 0;
    while (true)
    {
        iterations++;

        if (iterations >= 100)
        {
            logError logerror;
            logerror.LOGERROR("GasWaterEquilibria: Too many iterations!");
            exit(-1);
        }

        EOSPhaseProp eosProp;

        eosProp.initialize(gasCompositionName, zFeed, &phreeqc->gasComponentDatabase, GAS, temperatureK, pressureBar, PR78);
        eosProp.calcFugacity();

        int sizeOfgasComposition_H2OFree = sizeOfGasPhaseComponents - 1;

        for (int i = 0; i < sizeOfgasComposition_H2OFree; i++)
        {
            map<string, double> ::iterator it = update_timeStep.gasPhaseInfo.fugacity.find(gasCompositionName[i]);
            if (it == update_timeStep.gasPhaseInfo.fugacity.end())
            {
                update_timeStep.gasPhaseInfo.fugacity[gasCompositionName[i]] = eosProp.fugacity[i];
            }
            else
            {
                it->second = eosProp.fugacity[i];
            }

            update_timeStep.gasPhaseInfo.gasFugacityProvided = true;

            it = update_timeStep.gasPhaseInfo.moleNumber.find(gasCompositionName[i]);
            if (it == update_timeStep.gasPhaseInfo.moleNumber.end())
            {
                update_timeStep.gasPhaseInfo.moleNumber[gasCompositionName[i]] = gasPhaseMoleNumber_update[gasCompositionName[i]];
            }
        }

        update_timeStep.mineralsUpdate = mineralComp;

        phreeqc->run_simulation_timeStep(update_timeStep);
        gasPhase_results gas;
        waterSolution_results water;
        std::map<std::string, minerals> mineral;
        phreeqc->results_for_ReSoC(gas, water, mineral);

        map<string, minerals> ::iterator it1 = mineral.begin();
        map<string, double> ::iterator it2;
        for (; it1 != mineral.end(); it1++)
        {
            it2 = mineralComp.find(it1->first);
            if (it2 != mineralComp.end())
            {
                it2->second = it1->second.final_mole;

            }
            else
            {
                mineralComp[it1->first] = it1->second.final_mole;
            }
        }

        double yH2O = KH2O*water.activityOfH2O / pressureBar / exp(eosProp.logPhi[sizeOfgasComposition_H2OFree]);// eosProp.fugacity[sizeOfgasComposition_H2OFree];

        double maxDiffGasChange = 0.0;
        map<string, double>::iterator it = gasPhaseMoleNumber_update.begin();
        double totalGasMole = 0;
        for (; it != gasPhaseMoleNumber_update.end(); it++)
        {
            if (it->second > 0)
            {
                //                totalGasMole += it->second;
                if (gas.moleNumber.find(it->first) != gas.moleNumber.end())
                {
                    double updatedGasMole = (gas.moleNumber.find(it->first)->second);
                    totalGasMole += updatedGasMole;
                    double diff = fabs((updatedGasMole - (it->second)) / it->second);
                    if (diff > maxDiffGasChange)
                    {
                        maxDiffGasChange = diff;
                    }
                }
                else
                {
                    logError logerror;
                    logerror.LOGERROR("GasWaterEquilibria:: CANNOT find gas " + it->first + "! Please check!");
                    exit(-1);
                }
            }
        }

        if (maxDiffGasChange < 1.0e-5)
        {
            finalGas = gas;
            finalWater = water;
            finalMineral = mineral;
            yH2O_final = yH2O;
            return;
        }

        double moleNumberOfH2OinGas = totalGasMole / (1.0 - yH2O)*yH2O;
        double moleNumberOfH2OtoGas = moleNumberOfH2OinGas - moleNumberOfH2OinGas_old;
        moleNumberOfH2OinGas_old = moleNumberOfH2OinGas;
        double mass_of_H2O_old = water.mass_of_H2O;
        water.mass_of_H2O -= moleNumberOfH2OtoGas*moleWeightH2O / 1000.0; // KG;
        update_timeStep.mass_of_H2O = water.mass_of_H2O;

        gasPhaseMoleNumber_update = gas.moleNumber;
        update_timeStep.gasPhaseInfo.moleNumber = gasPhaseMoleNumber_update;

        for (int i = 0; i < sizeOfgasComposition_H2OFree; i++)
        {
            map<string, double> ::iterator iter = gasPhaseMoleNumber_update.find(gasCompositionName[i]);
            if (iter == gasPhaseMoleNumber_update.end())
            {
                logError logerror;
                logerror.LOGERROR("GasWaterEquilibria: CANNOT find component name from gasPhaseMoleNumber_update!");
                exit(-1);
            }
            else
            {
                zFeed[i] = iter->second / totalGasMole*(1.0 - yH2O);
            }
        }
        zFeed[sizeOfgasComposition_H2OFree] = yH2O;

        it = water.masterSpecies_moleNumber.begin();
        for (; it != water.masterSpecies_moleNumber.end(); it++)
        {
            map<string, double> ::iterator it2 = update_timeStep.masterChange.find(it->first);
            if (it2 != update_timeStep.masterChange.end())
            {
                it2->second = it->second / update_timeStep.mass_of_H2O;
            }
            else
            {
                update_timeStep.masterChange[it->first] = it->second / update_timeStep.mass_of_H2O;
            }
        }
    }
}

double ToolGasWaterEquilibria::KEquiliH2O(double temperature, double pressure)
{
    vector<double> coeff(7);
    if (temperature > 373.15)
    {
        coeff[0] = -0.902831272090587;
        coeff[1] = 0.0364929382599573;
        coeff[2] = 0.000436100194951838;
        coeff[3] = -3.10936036838763e-6;
        coeff[4] = 4.59205301435314e-9;
        coeff[5] = 16.2996873190268;
        coeff[6] = 0.0281119409320635;
    }
    else
    {
        coeff[0] = 9.31063597147498;
        coeff[1] = -0.189286700479115;
        coeff[2] = 0.00130713565196933;
        coeff[3] = -3.80022376294946e-6;
        coeff[4] = 4.00913697169677e-9;
        coeff[5] = 22.7692468626879;
        coeff[6] = -0.0112913301884868;
    }
    double equiliConstValue = (coeff[0] + coeff[1] * (temperature)+coeff[2] * temperature*temperature + coeff[3] * temperature*temperature*temperature + coeff[4] * pow(temperature, 4))*exp((pressure - 1.0)*(coeff[5] + coeff[6] * temperature) / temperature / R*0.1);
    return equiliConstValue;
}
