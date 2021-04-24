#include <iostream>
#include <map>

#include "../Phreeqc.h"
#include "../global_structures.h"
#include "GasWaterEquilibria.h"
#include "../logError.h"

GasWaterEquilibria::GasWaterEquilibria()
{
}

GasWaterEquilibria::~GasWaterEquilibria()
{
}

void GasWaterEquilibria::initializePhreeqc(string inputFileName)
{
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

void GasWaterEquilibria::calPhaseEquilibria(map<string, double> &gasPhaseMoleNumber, map<string, double> &masterSpecies, double temperatureK, double pressureBar)//, map<string, double> masterSpecies, map<string, double> minerals)
/*
    Here, gasPhaseMoleNumber is scaled each gas component mole numbers if H2O mass in water phase is 1 KG initially. 
*/
{
    int sizeOfGasPhaseComponents = (int)gasPhaseMoleNumber.size();
    vector<double> zFeed(sizeOfGasPhaseComponents);
    gasCompositionName.resize(sizeOfGasPhaseComponents);
    const double moleWeightH2O = 18.105; //  g/mole
    double totalMole = 0.0;
    bool H2OFound = false;
    string nameH2O="";
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
            zFeed[sizeOfGasPhaseComponents-1] = it->second/totalMole;
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
    phreeqc.gasComponentDatabase.initialize();
    conditionChange  update_timeStep;
    update_timeStep.pressure = pressureBar;
    update_timeStep.pressure *= atm_1bar;
    update_timeStep.temperature = temperatureK - 273.15;
    update_timeStep.mass_of_H2O = 1.0;

    update_timeStep.gasPhaseIn = TRUE;
    map <string, double> gasPhaseMoleNumber_update;
    gasPhaseMoleNumber_update = gasPhaseMoleNumber;
    gasPhaseMoleNumber_update.erase(nameH2O);
    double moleNumberOfH2OinGas_old = moleH2O;
    update_timeStep.masterChange = masterSpecies;

    iterations = 0;
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

        eosProp.initialize(gasCompositionName, zFeed, &phreeqc.gasComponentDatabase, GAS, temperatureK, pressureBar, PR78);
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

        phreeqc.run_simulation_timeStep(update_timeStep);
        gasPhase_results gas;
        waterSolution_results water;
        std::map<std::string, minerals> mineral;
        phreeqc.results_for_ReSoC(gas, water, mineral);
        double yH2O = KH2O*water.activityOfH2O / pressureBar / exp(eosProp.logPhi[sizeOfgasComposition_H2OFree]);// eosProp.fugacity[sizeOfgasComposition_H2OFree];
        
        double maxDiffGasChange = 0.0;
        map<string, double>::iterator it = gasPhaseMoleNumber_update.begin();
        double totalGasMole = 0;
        for (; it != gasPhaseMoleNumber_update.end(); it++)
        {
            if (it->second > 0)
            {
//                totalGasMole += it->second;
                if (gas.moleNumber.find(it->first ) != gas.moleNumber.end())
                {
                    double updatedGasMole = (gas.moleNumber.find(it->first )->second);
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
                update_timeStep.masterChange[it->first] = it->second/update_timeStep.mass_of_H2O;
            }
        }
    }
}

void GasWaterEquilibria::calPhaseEquilibria2(map<string, double> &gasPhaseMoleNumber, map<string, double> &masterSpecies, double temperatureK, double pressureBar)//, map<string, double> masterSpecies, map<string, double> minerals)
/*
Here, gasPhaseMoleNumber is scaled each gas component mole numbers if H2O mass in water phase is 1 KG initially.
*/
{
    int sizeOfGasPhaseComponents = (int) gasPhaseMoleNumber.size();
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
            zFeed[iComponent] = it->second / totalMole;
            gasCompositionName[iComponent] = it->first;
            iComponent++;
        }
        else
        {
            zFeed[sizeOfGasPhaseComponents - 1] = it->second / totalMole;
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
    phreeqc.gasComponentDatabase.initialize();
    conditionChange  update_timeStep;
    update_timeStep.pressure = pressureBar;
    update_timeStep.pressure *= atm_1bar;
    update_timeStep.temperature = temperatureK - 273.15;
    update_timeStep.mass_of_H2O = 1.0;

    update_timeStep.gasPhaseIn = TRUE;
    map <string, double> gasPhaseMoleNumber_update;
    gasPhaseMoleNumber_update = gasPhaseMoleNumber;
    gasPhaseMoleNumber_update.erase(nameH2O);
    double moleNumberOfH2OinGas_old = moleH2O;
    update_timeStep.masterChange = masterSpecies;

    iterations = 0;
    while (iterations<50)
    {
        iterations++;

        if (iterations >= 100)
        {
            logError logerror;
            logerror.LOGERROR("GasWaterEquilibria: Too many iterations!");
            exit(-1);
        }

        EOSPhaseProp eosProp;

        eosProp.initialize(gasCompositionName, zFeed, &phreeqc.gasComponentDatabase, GAS, temperatureK, pressureBar, PR78);
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

        phreeqc.run_simulation_timeStep(update_timeStep);
        gasPhase_results gas;
        waterSolution_results water;
        std::map<std::string, minerals> mineral;
        phreeqc.results_for_ReSoC(gas, water, mineral);
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

        map_water[iterations] = water;
        map_gas[iterations] = gas;
        map_mineral[iterations] = mineral;

        if (water.mass_of_H2O < 5.0e-6)
        {
            return;
        }

       /* if (maxDiffGasChange < 1.0e-5)
        {
            finalGas = gas;
            finalWater = water;
            yH2O_final = yH2O;
            return;
        }*/

        double moleNumberOfH2OtoGas = water.mass_of_H2O*0.3*1000.0 / moleWeightH2O;
        double moleNumberOfH2OinGas = moleNumberOfH2OinGas_old + moleNumberOfH2OtoGas;
        yH2O = moleNumberOfH2OinGas / (totalGasMole + moleNumberOfH2OinGas);
        double mass_of_H2O_old = water.mass_of_H2O;
        water.mass_of_H2O *= 0.7;
        update_timeStep.mass_of_H2O = water.mass_of_H2O;

        /*double moleNumberOfH2OinGas = totalGasMole / (1.0 - yH2O)*yH2O;
        double moleNumberOfH2OtoGas = moleNumberOfH2OinGas - moleNumberOfH2OinGas_old;
        moleNumberOfH2OinGas_old = moleNumberOfH2OinGas;*/

        //double mass_of_H2O_old = water.mass_of_H2O;
        //water.mass_of_H2O -= moleNumberOfH2OtoGas*moleWeightH2O / 1000.0; // KG;
        //update_timeStep.mass_of_H2O = water.mass_of_H2O;

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

        map<string, minerals> ::iterator it3 = mineral.begin();
        for (; it3 != mineral.end(); it3++)
        {
            update_timeStep.mineralsUpdate[it3->first] = it3->second.final_mole;
        }
    }
}


void GasWaterEquilibria::print_gasWaterEquilibria_screen()
{
    cout << "Iteration Numbers: " << iterations << endl;
    cout << endl;

    cout << "-------------- Water information --------------" << endl;
    cout << "Water mass: " << finalWater.mass_of_H2O << endl;
    cout << "Water density: " << finalWater.density << endl;
    cout << "Water activity: " << finalWater.activityOfH2O << endl;
    map<string, double> ::iterator it = finalWater.masterSpecies_molality.begin();
    for (; it != finalWater.masterSpecies_molality.end(); it++)
    {
        cout << "Master species " << it->first << ": " << it->second << endl;
    }
    cout << endl;

    cout << "-------------- Gas information ---------------" << endl;
    cout <<"Total mole number of gas: "<< finalGas.total_moles << endl;
 //   cout << endl;
    int gasSize = (int) finalGas.moleNumber.size();
    map<string, double>::iterator it2 = finalGas.moleNumber.begin();
    for (; it2 != finalGas.moleNumber.end(); it2++)
    {
        cout << "Gas component " << it2->first << " " << it2->second << endl;
    }
    cout << endl;
}

int GasWaterEquilibria::calKEquilConst(map<string,double> &moleFraction, double temperatureK, double pressureBar)
{
    LDBLE atm_1bar = 0.9869233;
    phreeqc.gasComponentDatabase.initialize();
    conditionChange  update_timeStep;
    update_timeStep.pressure = pressureBar;
    update_timeStep.pressure *= atm_1bar;
    update_timeStep.temperature = temperatureK - 273.15;
    update_timeStep.mass_of_H2O = 1.0;

    update_timeStep.gasPhaseIn = TRUE;
 //   update_timeStep.gasPhaseInfo.

    map<string, double> gasComposition_H2OFree;

    map<string, double> ::iterator it = moleFraction.begin();
    string nameH2O;
    double moleTotal = 0;
    double moleH2O = 0.0;
    bool H2OFound = false;
    double yH2Oini = 0.0;
    double KH2O = KEquiliH2O(temperatureK, pressureBar);
    for (; it != moleFraction.end(); it++)
    {
 //       moleTotal += it->second;

        if (it->second <= 0)
        {
            moleFraction.erase(it);
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
        }
        else
        {
            moleTotal += it->second;
        }
    }

    if (moleTotal <= 0)
    {
        logError logerror;
        logerror.LOGERROR("GasWaterEquilibria: Negative total mole gas fraction!");
        exit(-1);
    }

    it = moleFraction.begin();
    const double totalGasMole_H2OFree = 15.0; //mole
    for (; it != moleFraction.end(); it++)
    {
        if (it->first != nameH2O)
        {
            gasComposition_H2OFree[it->first] = (it->second / moleTotal) *totalGasMole_H2OFree;
        }
    }

    update_timeStep.gasPhaseInfo.moleNumber = gasComposition_H2OFree;

    int sizeOfgasComposition_H2OFree = (int) gasComposition_H2OFree.size();
    if (sizeOfgasComposition_H2OFree <= 0)
    {
        logError logerror;
        logerror.LOGERROR("GasWaterEquilibria: Error size of gas composition");
        exit(-1);
    }

    gasCompositionName.resize(sizeOfgasComposition_H2OFree+1);
    vector<double>    gasMoleFraction(sizeOfgasComposition_H2OFree+1);

    int i = 0;
    it = gasComposition_H2OFree.begin();
    for (; it != gasComposition_H2OFree.end(); it++)
    {
        gasCompositionName[i] = it->first;
        gasMoleFraction[i] = it->second / (totalGasMole_H2OFree)*(1.0-yH2Oini);
        i++;
    }

    gasCompositionName[sizeOfgasComposition_H2OFree] = "H2O";
    gasMoleFraction[sizeOfgasComposition_H2OFree] = yH2Oini;

    vector<double> KEquili_old(sizeOfgasComposition_H2OFree + 1);
    for (int i = 0; i < sizeOfgasComposition_H2OFree + 1; i++)
    {
        KEquili_old[i] = 1.0;
    }
    int iterations = 0;
    while (true)
    {
        EOSPhaseProp eosProp;

        eosProp.initialize(gasCompositionName, gasMoleFraction, &phreeqc.gasComponentDatabase, GAS, update_timeStep.temperature, update_timeStep.pressure, PR78);
        eosProp.calcFugacity();
        eosProp.calDensity(false);

        for (int i = 0; i < sizeOfgasComposition_H2OFree; i++)
        {
            update_timeStep.gasPhaseInfo.fugacity[gasCompositionName[i]] = eosProp.fugacity[i];
        }

        phreeqc.run_simulation_timeStep(update_timeStep);
        gasPhase_results gas;
        waterSolution_results water;
        std::map<std::string, minerals> mineral;
        phreeqc.results_for_ReSoC(gas, water, mineral);
        map<string, double> gasMoleInWater;
        double moleNumberTotalOfWater = 0;
        double yH2O = KH2O*water.activityOfH2O/pressureBar/eosProp.fugacity[sizeOfgasComposition_H2OFree];

        if (gas.total_moles > 0)
        {
            map<string, double>::iterator it = gas.moleNumber.begin();
            for (; it != gas.moleNumber.end(); it++)
            {
                double gasMoleTemp = (gasComposition_H2OFree.find(it->first)->second) - (it->second);
                gasMoleInWater[it->first] = gasMoleTemp;
                moleNumberTotalOfWater += gasMoleTemp;
            }
        }

        double moleNumberOfH2O = water.mass_of_H2O*1000.0 / 18.015;
        moleNumberTotalOfWater += moleNumberOfH2O;

        for (int i = 0; i < sizeOfgasComposition_H2OFree; i++)
        {
            double Ktemp = gas.moleNumber.find(gasCompositionName[i])->second / gas.total_moles*(1.0 - yH2O);
            gasMoleFraction[i] = Ktemp;
            Ktemp /= (gasMoleInWater.find(gasCompositionName[i])->second) / moleNumberTotalOfWater;           
            KEquiliConst[i] = Ktemp;
        }
        gasMoleFraction[sizeOfgasComposition_H2OFree] = yH2O;
        KEquiliConst[sizeOfgasComposition_H2OFree] = yH2O / moleNumberOfH2O*moleNumberTotalOfWater;

        double maxKeqiliDiff = fabs((KEquiliConst[sizeOfgasComposition_H2OFree] - KEquili_old[sizeOfgasComposition_H2OFree]) / KEquili_old[sizeOfgasComposition_H2OFree]);
        for (int i = 0; i < sizeOfgasComposition_H2OFree; i++)
        {
            double tempKequiDiff = fabs((KEquiliConst[sizeOfgasComposition_H2OFree] - KEquili_old[i]) / KEquili_old[i]);
            if (maxKeqiliDiff < tempKequiDiff)
            {
                maxKeqiliDiff = tempKequiDiff;
            }
        }

        if (maxKeqiliDiff < 1.0e-5)
        {
            return iterations;
        }

        KEquili_old = KEquiliConst;
        iterations++;

        update_timeStep.gasPhaseInfo.moleNumber = gas.moleNumber;
        update_timeStep.masterChange = water.masterSpecies_molality;
    }

}

double GasWaterEquilibria::KEquiliH2O(double temperature, double pressure)
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

bool GasWaterEquilibria::calcGasMole(double &beta_GasMole, vector<double> &zFeed, vector<double> &Kvalue)
{
    double RReqn0, RReqn1;
    double beta0 = 0.0, beta1 = 1.0;
    int i;

    if (zFeed.size() != Kvalue.size())
    {
        cerr << "Error size of composition or K value!" << endl;
    }

    RReqn0 = calRReqn(beta0, zFeed, Kvalue);
    RReqn1 = calRReqn(beta1, zFeed, Kvalue);

    if (RReqn0 <= 0.0)
    {
        beta_GasMole = -1.0;
        return false;
    }

    if (RReqn1 >= 0.0)
    {
        beta_GasMole = 2.0;
        return false;
    }

    //Find tighter boundary..
    double betaMin = 0.0, betaMax = 1.0;
    double LowerBound, HigherBound;
    int ComponentNumber = (int) zFeed.size();

    for (i = 0; i < ComponentNumber; i++)
    {
        if (Kvalue[i] > 1.0)
        {
            LowerBound = (Kvalue[i] * zFeed[i] - 1.0) / (Kvalue[i] - 1.0);
            if (betaMin < LowerBound) betaMin = LowerBound;
        }

        if (Kvalue[i] < 1.0)
        {
            HigherBound = (1.0 - zFeed[i]) / (1.0 - Kvalue[i]);
            if (betaMax > HigherBound) betaMax = HigherBound;
        }
    }

    beta_GasMole = 0.5*(betaMax + betaMin);

    //Find beta_GasMole by Newton iteration!
    const double eps_beta = 1.0e-6;
    const int MaxIteration_beta = 10000;
    int iterNum = 0;
    double RReqn_Current, RReqnPrim_Current;
    double beta_GasMole_new;
    while (true)
    {
        if (iterNum > MaxIteration_beta)
        {
            //          if (fabs(RReqn_Current) < 1.e-4)
            //         {
            return true;
            //         }
            cerr << "Too many iterations for gas mole fraction calculation!!!" << endl;
            exit(-1);
        }

        RReqn_Current = calRReqn(beta_GasMole, zFeed, Kvalue);
        //        cout << RReqn_Current << endl;

        if (fabs(RReqn_Current) <= eps_beta)
        {
            return true;
        }

        if (RReqn_Current > 0)
        {
            betaMin = beta_GasMole;
        }
        else
        {
            betaMax = beta_GasMole;
        }

        iterNum += 1;

        RReqnPrim_Current = calRReqnPrim(beta_GasMole, zFeed, Kvalue);

        //      cout << "RReqn_Current = " << RReqn_Current << endl;

        if (RReqnPrim_Current == 0.0)
        {
            cerr << "RReqnPrim is 0!!!" << endl;
            exit(1);
        }
        beta_GasMole_new = beta_GasMole + (-RReqn_Current / RReqnPrim_Current);
        //    cout << "beta_Gasmole = " << beta_GasMole_new << endl;
        if (beta_GasMole_new<betaMax && beta_GasMole_new>betaMin)
        {
            beta_GasMole = beta_GasMole_new;
        }
        else
        {
            beta_GasMole = 0.5*(betaMin + betaMax);
        }
    }
}

double GasWaterEquilibria::calRReqn(double &beta, vector<double>&zFeed, vector<double> &Kvalue)
{
    double RReqn = 0.0;
    if (Kvalue.size() != zFeed.size())
    {
        cerr << "Error size of composition or K value!" << endl;
    }

    int i;
    int ComponentNumber =(int) zFeed.size();
    for (i = 0; i < ComponentNumber; i++)
    {
        RReqn += zFeed[i] * (Kvalue[i] - 1.0) / (1.0 - beta + beta*Kvalue[i]);
    }

    return RReqn;
}

double GasWaterEquilibria::calRReqnPrim(double &beta, vector<double>&zFeed, vector<double>&Kvalue)
{
    double RReqnPrim = 0.0;
    if (Kvalue.size() != zFeed.size())
    {
        cerr << "Error size of composition or K value!" << endl;
    }

    int i;
    int ComponentNumber = (int)zFeed.size();
    for (i = 0; i < ComponentNumber; i++)
    {
        RReqnPrim -= zFeed[i] * (Kvalue[i] - 1.0)*(Kvalue[i] - 1.0) / (1.0 - beta + beta*Kvalue[i]);
    }

    return RReqnPrim;
}