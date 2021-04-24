
#include <iostream>
#include <map>
#define testRun 0

#if testRun
#include "Phreeqc.h"
#include "global_structures.h"
#include "CO2_solubility_in_pure_water/CO2SoluInWater.h"
#include "CO2_solubility_in_NaCl_solution/CO2SoluInBrine.h"
#include "O2_solubility_in_pure_water/O2SolubilityInWater.h"
#include "O2_solubility_in_brine/O2SolubilityInBrine.h"
#include "CH4_Solubility/CH4SolubilityInWater.h"
#include "CH4_Solubility/CH4SolubilityInBrine.h"
#include "MultiPhaseEquilibria/GasWaterEquilibria.h"
#include "N2_Solubility/N2SolubilityWater.h"
#include "N2_Solubility/N2SolubilityInBrine.h"
#include "Shi_2013Reproduce/ShiReproduce.h"
#include "Liu_2012Reproduce/CO2N2withPureWater.h"
#include "Ellis1959_1963/EllisReproduce.h"
#include "MgCO3/MagnesiteSolubility.h"
#include "MgCO3/Nesquehonite.h"
#include "MgCO3/Lansfordite.h"
#include "MgCO3/MgCO3Coexist.h"
#include "Caclite_solubilty/CalciteSolubility.h"
#include "MultiPhaseEquilibria/WaterMineralEquilibria.h"
#include "DolomiteSolubility/DolomiteSolubility.h"
#include "Dolomitization/Dolomitization.h"
#include "Kinetics/KineticsExample.h"
#include "Kinetics/QuestBatch.h"
#include "H2S_solubility/H2S_solubility.h"
#endif

#include "General/GeoChemCalc.h"
#include "InputManager.h"

using namespace std;

void main(int argc, char** argv)
{
#if !testRun
    string inputFileName;
    if (argc > 1)
    {
        inputFileName = argv[1];
    }
    else
    {
        inputFileName = "GeoChemInput.dat";
    }
    InputManager inputs;

    if (!inputs.readInputFile(inputFileName))
    {
        cerr << "ERROR: Reading input file error! Please check the inputs!" << endl;
        exit(-1);
    }
    
    GeoChemCalc geoChemSim;

    geoChemSim.initilize(inputs.getInputDataBlock());
  //  geoChemSim.calculate(inputFileName);
    geoChemSim.calculateWithIniWater(inputFileName);

#else
    CO2SoluInWater CO2solubility;
    double temperatureK = 353.15;
    double pressureBar = 100;
    double volumeL = 20.0;
//    CO2solubility.calCO2SolubilityInPureWater(temperatureK, pressureBar, volumeL);
//    CO2solubility.calCO2SolubilityInPureWater("CO2_in_Water_exp.txt");
    CO2SoluInBrine CO2solubilityBrine;
    double mNaCl = 4.0;
//    CO2solubilityBrine.calCO2Solubility(temperatureK, pressureBar, mNaCl, volumeL);
//    CO2solubilityBrine.calCO2Solubility();

    O2SolubilityInWater O2Solubility; 
//    O2Solubility.calO2Solubility(temperatureK, pressureBar, volumeL);
//    O2Solubility.calO2Solubility();
    O2SolubilityInBrine  O2solubilityBrine;
 //   O2solubilityBrine.calO2SolubilityInBrine(temperatureK, pressureBar, mNaCl, volumeL);
//    O2solubilityBrine.calO2SolubilityInBrine();

    CH4SolubilityInwater CH4solubility;
 //   CH4solubility.calCH4SolubilityInWater(temperatureK, pressureBar, volumeL);

    CH4SolubilityInBrine CH4solubilityBrne;

    temperatureK = 398.15;
    pressureBar = 100.0;
    mNaCl = 4.6;

    int iChoice = 35;
    if (iChoice == 1)
    {
        N2SolubilityWater N2SolubilityInWater;
        N2SolubilityInWater.initialize();
        cout << "T(K)    " << "P(bar)    " << "N2 solubility" << endl;
        while (true)
        {
            if (pressureBar > 660)
            {
                return;
            }

            double N2Solubility = N2SolubilityInWater.calN2Solubility(temperatureK, pressureBar);
            cout << temperatureK << "  " << pressureBar << "  " << N2Solubility << endl;

            pressureBar += 30.0;
        }
    }
    else if (iChoice == 2)
    {
        N2SolubilityWater N2SolubilityInWater;
        N2SolubilityInWater.initialize();
        double N2Solubility = N2SolubilityInWater.calN2Solubility(temperatureK, pressureBar);
        cout << "N2 solubility in molality (mole/KgW): " << N2Solubility << " !" << endl;
        cout << "H2O solubility in N2 rich phase: " << N2SolubilityInWater.m_yH2O << " !" << endl;
    }
    else if (iChoice == 3)
    {
        N2SolubilityInBrine N2solubility;
        N2solubility.initialize(); 
        cout << "T(K)    " << "P(bar)    " <<"mNaCl   "<< "N2 solubility" << endl;
        while (true)
        {
            if (pressureBar > 660)
            {
                return;
            }

            double N2Solu = N2solubility.calN2Solubility(temperatureK, pressureBar, mNaCl);
            cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << N2Solu << endl;

            pressureBar += 30.0;
        }
    }
    else if (iChoice == 4)
    {
        O2SolubilityInBrine O2solubility;
        O2solubility.initialize();
        mNaCl = 0.0;
        cout << "T(K)    " << "P(bar)    " << "mNaCl   " << "O2 solubility" << endl;
        temperatureK = 533;
        pressureBar = 50;
        while (true)
        {
            if (pressureBar > 210)
            {
                return;
            }

            double O2solu = O2solubility.calO2Solubility_GasWaterEquili(temperatureK, pressureBar, mNaCl);
            cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << O2solu << endl;
            pressureBar += 5.0;
        }
    }
    else if (iChoice == 5)
    {
        O2SolubilityInBrine O2solubility;
        O2solubility.initialize();
        mNaCl = 0.0;
        cout << "T(K)    " << "P(bar)    " << "mNaCl   " << "O2 solubility" << endl;
        temperatureK = 273;
        pressureBar = 1.0;
        while (true)
        {
            if (mNaCl > 2.6)
            {
                return;
            }

            double O2solu = O2solubility.calO2Solubility_GasWaterEquili(temperatureK, pressureBar, mNaCl);
            cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << O2solu << endl;
            mNaCl += 0.2;
        }
    }
    else if (iChoice == 6)
    {
        ShiReproduce shiReproduce;
        waterSolution_results water;
        gasPhase_results gas;
        map<string, minerals> mineral;
        double temperatureK = 273.15+42.0;
        double pressureBar = 41.37;

        double yH2O;
        shiReproduce.initialize();
        shiReproduce.CO2saturation(temperatureK, pressureBar, water, gas, yH2O);
        /*temperatureK = 273.15;
        pressureBar = 34.5;
        shiReproduce.calciteDisolution(temperatureK, pressureBar, gas, water, mineral);*/

    }
    else if (iChoice == 7)
    {
        CO2N2withPureWater CO2N2H2O;
        map<string, double> gasFeed;
        map<string, double> masterSpecies;
        gasFeed["N2"] = 0.123;
        gasFeed["CO2"] = 0.877;
        temperatureK = 308.15;
        pressureBar = 80.0;
        CO2N2H2O.initialize();
        CO2N2H2O.calCO2N2H2Oeqilibria(temperatureK, pressureBar, gasFeed, masterSpecies);
        cout << "xN2 : " << CO2N2H2O.waterPhaseMoleFraction.find("N2")->second << endl;
        cout << "xCO2 : " << CO2N2H2O.waterPhaseMoleFraction.find("CO2")->second << endl;
        cout << "xH2O : " << CO2N2H2O.waterPhaseMoleFraction.find("H2O")->second << endl;
    }
    else if (iChoice == 8)
    {
        CO2N2withPureWater CO2N2H2O;
        CO2N2H2O.initialize();

        ifstream inFile("Liu_2012-Data.txt");
        double gasN2, gasCO2;
        double waterN2, waterCO2, waterH2O;
        temperatureK = 308.15;
        pressureBar = 80.0;
        cout << endl;
        while (inFile >> gasN2 >> gasCO2 >> pressureBar >> temperatureK >> waterN2 >> waterCO2 >> waterH2O)
        {
            map<string, double> gasFeed;
            map<string, double> masterSpecies;
            gasFeed["N2"] = gasN2;
            gasFeed["CO2"] = gasCO2;
            CO2N2H2O.calCO2N2H2Oeqilibria(temperatureK, pressureBar, gasFeed, masterSpecies);
            cout << CO2N2H2O.waterPhaseMoleFraction.find("N2")->second << " " << CO2N2H2O.waterPhaseMoleFraction.find("CO2")->second << " " << CO2N2H2O.waterPhaseMoleFraction.find("H2O")->second << endl;
        }
        cout << endl;
    }
    else if (iChoice == 9)
    {
        CO2N2withPureWater CO2N2H2O;
        CO2N2H2O.initialize();

        ifstream inFile("Liu_2012_with_salts.txt");
//        double waterN2, waterCO2, waterH2O;
        double gasN2, gasCO2;
        map<string, double> masterSpecies;
        /*masterSpecies["Na"] = 0.059979;
        masterSpecies["K"] = 0.047098;
        masterSpecies["Ca"] = 0.031611;
        masterSpecies["Cl"] = 0.170298;*/

        /*masterSpecies["Na"] = 0.05698;
        masterSpecies["K"] = 0.044743;
        masterSpecies["Ca"] = 0.03003;
        masterSpecies["Cl"] = 0.161783;*/

        masterSpecies["Na"] = 0.170298;

        masterSpecies["Cl"] = 0.170298;

        temperatureK = 308.15;
        pressureBar = 80.0;
        cout << endl;
        while (inFile >> gasN2 >> gasCO2)
        {
            map<string, double> gasFeed;
            gasFeed["N2"] = gasN2;
            gasFeed["CO2"] = gasCO2;
            CO2N2H2O.calCO2N2H2Oeqilibria(temperatureK, pressureBar, gasFeed, masterSpecies);
            cout << CO2N2H2O.waterPhaseMoleFraction.find("N2")->second << " " << CO2N2H2O.waterPhaseMoleFraction.find("CO2")->second << " " << CO2N2H2O.waterPhaseMoleFraction.find("H2O")->second << endl;
        }
        cout << endl;
    }
    else if (iChoice==10)
    {
        EllisReproduce ellis;
        ellis.initialize();
        double pressureBar = 12;
        double temperatureK = 404;
        double mNaCl = 0.0;
        double calciteMolality = 0.0;
        cout << endl;
        while (true)
        {
            if (mNaCl > 1.2)
            {
                cout << endl;
                return;
            }

            ellis.calcCalciteSolubility(temperatureK, pressureBar, mNaCl, calciteMolality);
            cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << calciteMolality << endl;

            mNaCl += 0.1;
        }
    }
    else if (iChoice == 11)
    {
        EllisReproduce ellis;
        ellis.initialize();
        double pressureBar = 1;
        double temperatureK = 298.15;
        double mNaCl = 0.0;
        double calciteMolality = 0.0;
        cout << endl;
        while (true)
        {
            if (temperatureK > 530)
            {
                cout << endl;
                return;
            }

            ellis.calcCalciteSolubility(temperatureK, pressureBar, mNaCl, calciteMolality);
            cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << calciteMolality << endl;

            temperatureK += 10;
        }
    }
    else if (iChoice == 12)
    {
        MagnesiteSolubility mgsSolubility;
        mgsSolubility.initialize();
        double pressureBar = 5;
        double temperatureK = 150+273.15;
        double mNaCl = 0.1;
        double calciteMolality = 0.0;
        double Mg2Molality = 0.0;
        double pH = 7.0;
        cout << endl;
        while (true)
        {
            if (pressureBar > 40)
            {
                cout << endl;
                return;
            }

            mgsSolubility.calcMagnesiteSolubility(temperatureK, pressureBar, mNaCl, 1.0, calciteMolality, Mg2Molality, pH);
            cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << calciteMolality << "  " << Mg2Molality << endl;

            pressureBar += 2.5;
        }
    }
    else if (iChoice == 13)
    {
        MagnesiteSolubility mgsSolubility;
        mgsSolubility.initialize();
        double pressureBar = 4;
        double temperatureK = 50 + 273.15;
        double mNaCl = 0.1;
        double calciteMolality = 0.0;
        cout << endl;
        
        mgsSolubility.calcMagnesiteSolubilityNoCO2(temperatureK, pressureBar, mNaCl, calciteMolality);
        cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << calciteMolality << endl;
    }
    else if (iChoice == 14)
    {
        MagnesiteSolubility mgsSolubility;
        mgsSolubility.initialize();
        double pressureBar = 4;
        double temperatureK = 50 + 273.15;
        double mNaCl = 0.1;
        double calciteMolality = 0.0;
        double Mg2Molality = 0.0;
        double pH = 7.0;
        cout << endl;

        mgsSolubility.calcMagnesiteSolubility(temperatureK, pressureBar, mNaCl, 1.0, calciteMolality, Mg2Molality, pH);
        cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << calciteMolality << "  " << Mg2Molality << endl;
    }
    else if (iChoice == 15)
    {
        CalciteSolubility calSolubility;
        calSolubility.initialize();
        map<string, double> gasFeed;
        map<string, double> masterSpecies;
        map<string, double> minralComp;
        waterSolution_results water;

        double temperatureK = 70+273.15;
        double pressureBar = 1000.0;
        double mNaCl = 0.0;
        
        double iniGas = 10.0;
//        gasFeed["CO2"] = iniGas*1.0;
//        gasFeed["N2"] = iniGas*0.1;

        minralComp["Magnesite"] = 1.0;

        cout << endl;

        while (true)
        {
            if (mNaCl > 6)
            {
                cout << endl;
                return;
            }

            masterSpecies["Na"] = mNaCl;
            masterSpecies["Cl"] = mNaCl;

            int sizeOfGas = (int) gasFeed.size();
            if (sizeOfGas < 1)
            {
                calSolubility.calCalciteSolubility(temperatureK, pressureBar, masterSpecies, minralComp, water);
 //               continue;
            }
            else
            {
              calSolubility.calCalciteSolubility(temperatureK, pressureBar, gasFeed, masterSpecies, minralComp, water);
            }

            map<string, double> ::iterator it = water.masterSpecies_molality.find("Mg");
            double calciteSolubility = 0.0; 
            if (it != water.masterSpecies_molality.end())
            {
                calciteSolubility = it->second;
            }
            else
            {
                cerr << "CANNOT find calcium species!" << endl;
                return;
            }

            double CO2solu = 0.0;
            it = water.masterSpecies_molality.find("C(4)");
            if (it != water.masterSpecies_molality.end())
            {
                CO2solu = it->second - calciteSolubility;
            }

            cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << water.ph << "  " << CO2solu << "  " << calciteSolubility << endl;

            mNaCl += 0.2;
        }
    }
    else if (iChoice == 16)
    {
        double temperatureK = 247.5 + 273.15;
        double pressureBar = 1.0*PASCAL_PER_ATM / 1.e5;
        double mNaCl = 0.1;
        
        DolomiteSolubility dolomite;
        waterSolution_results water;
        dolomite.initialize();
        dolomite.calcDolomiteSolubilityExpReproduce(temperatureK, pressureBar, mNaCl, water);
        double logMg2 = 0.0;
        map<string, double>::iterator it = water.species_molality.find("Mg+2");
        if (it != water.species_molality.end())
        {
            logMg2 = it->second;
            logMg2 = log10(logMg2);
        }
        else
        {
            cerr << "CANNOT find Mg+2" << endl;
            return;
        }

        double logCa2 = 0.0;
        it = water.species_molality.find("Ca+2");
        if (it != water.species_molality.end())
        {
            logCa2 = it->second;
            logCa2 = log10(logCa2);
        }
        else
        {
            cerr << "CANNOT find Ca+2" << endl;
            return;
        }

        double logH1 = 0.0;
        it = water.species_molality.find("H+");
        if (it != water.species_molality.end())
        {
            logH1 = it->second;
            logH1 = log10(logH1);
        }
        else
        {
            cerr << "CANNOT find H+" << endl;
            return;
        }

        cout << temperatureK << " " << pressureBar << " " << -logH1 << " " << logCa2 << " " << logMg2 << endl;
        cout << endl;
    }
    else if (iChoice == 17)
    {
        double temperatureK = 70 + 273.15;
        double pressureBar = 1000.0;
        double mNaCl = 0.0;

        DolomiteSolubility dolomiteWithCO2;
        waterSolution_results water;
        dolomiteWithCO2.initializeCO2();
        map<string, double> gasMole;
        gasMole["CO2"] = 10.0;

        while (true)
        {
            if (mNaCl > 6)
            {
                cout << endl;
                return;
            }
            dolomiteWithCO2.calcDolomiteSolubilityWithCO2(temperatureK, pressureBar, mNaCl, gasMole, water);
            double Mg2 = 0.0;
            double Ca2 = 0.0;

            map<string, double> ::iterator it = water.masterSpecies_molality.find("Mg");
            if (it != water.masterSpecies_molality.end())
            {
                Mg2 = it->second;
            }
            else
            {
                cerr << "CANNOT find Mg" << endl;
                return;
            }

            it = water.masterSpecies_molality.find("Ca");
            if (it != water.masterSpecies_molality.end())
            {
                Ca2 = it->second;
 //          logCa2 = log10(logCa2);
            }
            else
            {
                cerr << "CANNOT find Ca+2" << endl;
                return;
            }

            cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << Ca2 << "  " << Mg2 << endl;
            mNaCl += 0.2;
        }

    }
    else if (iChoice == 18)
    {
        double temperatureK = 298.15;
        double pressureBar = 20.0;
        map<string, double> mstrSpecies;
        mstrSpecies["Mg"] = 0.1;
        mstrSpecies["Ca"] = 0.1;
        mstrSpecies["C(4)"] = 0.2;

        while (true)
        {
            if (pressureBar > 1000)
            {
                cout << endl;
                return;
            }
            map<string, double> mineralComp;
            mineralComp["Calcite"] = 0.0;
            Dolomitization dolomi;
            dolomi.initialize();
            dolomi.dolomitization(temperatureK, pressureBar, mstrSpecies, mineralComp);
            map<string, minerals>::iterator it = dolomi.mineralFinal.find("Dolomite");
            double dolomiteComp = 0.0;
            double calciteComp = 0.0;
            double mgsComp = 0.0;
            double nesComp = 0.0;
            if (it != dolomi.mineralFinal.end())
            {
                dolomiteComp = it->second.final_mole;
            }
            it = dolomi.mineralFinal.find("Calcite");
            if (it != dolomi.mineralFinal.end())
            {
                calciteComp = it->second.final_mole;
            }
            it = dolomi.mineralFinal.find("Magnesite");
            if (it != dolomi.mineralFinal.end())
            {
                mgsComp = it->second.final_mole;
            }
            it = dolomi.mineralFinal.find("Nesquehonite");
            {
                nesComp = it->second.final_mole;
            }
            cout << temperatureK << "  " << pressureBar << "  " << dolomiteComp << "  " << calciteComp << "  " << mgsComp << "  " << nesComp << endl;
            pressureBar += 20.0;
        }
    }
    else if (iChoice==19)
    {
        double temperatureK = 298.15;
        double pressureBar = 20.0;
        map<string, double> mstrSpecies;
        mstrSpecies["Mg"] = 0.1;
        mstrSpecies["Na"] = 0.2;
        mstrSpecies["Cl"] = 0.4;

        while (true)
        {
            if (pressureBar > 1000)
            {
                cout << endl;
                return;
            }
            map<string, double> mineralComp;
            mineralComp["Calcite"] = 1.0;
            Dolomitization dolomi;

            map<string, double> gasFeed;
            gasFeed["CO2"] = 10.0;

            dolomi.initializeCO2();
            dolomi.dolomitizationCO2(temperatureK, pressureBar, gasFeed, mstrSpecies, mineralComp);
            map<string, minerals>::iterator it = dolomi.mineralFinal.find("Dolomite");
            double dolomiteComp = 0.0;
            double calciteComp = 0.0;
            double mgsComp = 0.0;
            double nesComp = 0.0;
            if (it != dolomi.mineralFinal.end())
            {
                dolomiteComp = it->second.final_mole;
            }
            it = dolomi.mineralFinal.find("Calcite");
            if (it != dolomi.mineralFinal.end())
            {
                calciteComp = it->second.final_mole;
            }
            it = dolomi.mineralFinal.find("Magnesite");
            if (it != dolomi.mineralFinal.end())
            {
                mgsComp = it->second.final_mole;
            }
            it = dolomi.mineralFinal.find("Nesquehonite");
            {
                nesComp = it->second.final_mole;
            }
            cout << temperatureK << "  " << pressureBar << "  " << dolomiteComp << "  " << calciteComp <<"  "<<mgsComp<< "  " << nesComp << endl;
            pressureBar += 20.0;
        }
    }
    else if (iChoice == 20)
    {
        Phreeqc phrqc;
        string inputFile = "CO2SolubilityWithMineral.phr";
        ToolGasWaterEquilibria CO2WaterMineralEquili;
        CO2WaterMineralEquili.initialize(&phrqc, inputFile);

        map<string, double> gasMoleNumber;
        map<string, double> masterSpecies;
        map<string, double> mineralComp;

        gasMoleNumber["N2"] = 2.0;
        gasMoleNumber["CO2"] = 8.0;
        mineralComp["Calcite"] = 1.0;

        double temperatureK = 373.15;
        double pressureBar = 300.0;
        double mNaCl = 0.0;
        gasPhase_results finalGas;
        waterSolution_results finalWater;
        map<string, minerals> finalMineral;
        double yH2O;

        while (true)
        {
            if (mNaCl > 6)
            {
                cout << endl;
                return;
            }
            mineralComp["Calcite"] = 1.0;
            masterSpecies["Na"] = mNaCl;
            masterSpecies["Cl"] = mNaCl;
            CO2WaterMineralEquili.calPhaseEquilibria(&phrqc, gasMoleNumber, masterSpecies, mineralComp, temperatureK, pressureBar, finalGas, finalWater, finalMineral, yH2O);
            double CO2molality = finalWater.masterSpecies_molality.find("C(4)")->second;
            double mineralMolality = finalWater.masterSpecies_molality.find("Ca")->second;
            double O2molality = finalWater.masterSpecies_molality.find("N(0)")->second / 2.0;
            double densityWater = finalWater.density;
            double ph = finalWater.ph;
            cout << pressureBar << "  " << mNaCl << "  " <<ph<<"  " << CO2molality<<"  "<<O2molality <<"  "<<densityWater << "  " << mineralMolality << endl;
 //           cout << pressureBar << "  " << CO2molality << "  " << mineralMolality << endl;
            mNaCl += 0.2;
        }
    }
    else if (iChoice == 21)
    {
        Phreeqc phrqc;
        CO2SoluInBrine co2Water;
        double temperatureK = 328.15;
        double pressureBar = 100.0;
        mNaCl = 3.0;
        double CO2solu=   co2Water.calCO2SolubilityNewMethod(temperatureK, pressureBar, mNaCl);
        cout << "CO2 solubility:" << CO2solu << endl;
        cout << "H2O solubility in gas phase:" << co2Water.yH2O << endl;
    }
    else if (iChoice == 22)
    {
        ifstream expData("CO2_in_Brine_exp.txt");
        ofstream calData("CO2_brine_calc_phreeqcNew_2018_12.txt");

        int iLiterature;
        double pressure, temperature, mNaCl, mCO2;


        if (!expData.eof())
        {
            expData.ignore(2000, '\n');
        }

        int i = 0;
        while (expData >> iLiterature >> temperature >> pressure >> mNaCl >> mCO2)
        {
            i++;
            CO2SoluInBrine co2Water;
            double CO2solu = co2Water.calCO2SolubilityNewMethod(temperature, pressure, mNaCl);
            calData << i << "  " << iLiterature << "  " << temperature << "  " << pressure << "  " << mNaCl << "  " << mCO2 << "  " << CO2solu << endl;
            cout << i << "  " << iLiterature << "  " << temperature << "  " << pressure << "  " << mNaCl << "  " << mCO2 << "  " << CO2solu << endl;
        }
    }
    else if (iChoice == 23)
    {
        ifstream expData("CO2_in_Water_exp.txt");
        ofstream calData("CO2_water_calc_phreeqcNew_2018_12.txt");

        int iLiterature;
        double pressure, temperature,  mCO2;


        if (!expData.eof())
        {
            expData.ignore(2000, '\n');
        }

        int i = 0;
        while (expData >> iLiterature >> temperature >> pressure >> mCO2)
        {
            i++;
            CO2SoluInBrine co2Water;
            double CO2solu = co2Water.calCO2SolubilityNewMethod(temperature, pressure, 0.0);
            calData << i << "  " << iLiterature << "  " << temperature << "  " << pressure << "  " << mCO2 << "  " << CO2solu << endl;
            cout << i << "  " << iLiterature << "  " << temperature << "  " << pressure << "  "  << mCO2 << "  " << CO2solu << endl;
        }
    }
    else if (iChoice == 24)
    {
        ifstream expData("CH4_in_Water_exp.txt");
        ofstream calData("CH4_water_calc_phrqcNew_2018_12.txt");

        int iLiterature;
        double pressure, temperature, mCH4;
        if (!expData.eof())
        {
            expData.ignore(2000, '\n');
        }

        int i = 0;
        while (expData>>iLiterature>>temperature>>pressure>>mNaCl>>mCH4)
        {
            i++;
            CH4SolubilityInBrine ch4Water;
            double CH4solu = ch4Water.calCH4SolubilityInBrineNewMethod(temperature, pressure, mNaCl);
            calData << i << " " << iLiterature << " " << temperature << " " << pressure << " " << mCH4 << " " << CH4solu << endl;
            cout << i << " " << iLiterature << " " << temperature << " " << pressure << " " << mCH4 << " " << CH4solu << endl;
        }
    }
    else if (iChoice == 25)
    {
        ifstream expData("CH4_Solu_Brine_Data.txt");
        ofstream calData("CH4_Brine_calc_phrqcNew_2018_12.txt");

        int iLiterature;
        double pressure, temperature, mCH4;
        if (!expData.eof())
        {
            expData.ignore(2000, '\n');
        }

        int i = 0;
        while (expData >> iLiterature >> temperature >> pressure >> mNaCl >> mCH4)
        {
            i++;
            CH4SolubilityInBrine ch4Water;
            double CH4solu = ch4Water.calCH4SolubilityInBrineNewMethod(temperature, pressure, mNaCl);
            calData << i << " " << iLiterature << " " << temperature << " " << pressure<<" "<<mNaCl << " " << mCH4 << " " << CH4solu << endl;
            cout << i << " " << iLiterature << " " << temperature << " " << pressure << " " << mNaCl << " " << mCH4 << " " << CH4solu << endl;
        }
    }
    else if (iChoice == 26)
    {
        ifstream expData("N2_solubility_in_water_Duan.txt");
        ofstream calData("N2_solubility_in_water_calc_phrqcNew_2018_12.txt");

        double pressure, temperature, mN2;
        if (!expData.eof())
        {
            expData.ignore(2000, '\n');
        }

        int i = 0;
        mNaCl = 0;
        while (expData >> pressure >> temperature >> mN2)
        {
            i++;
            N2SolubilityInBrine N2Water;
            double N2solu = N2Water.calN2SolubilityInBrineNewMethod(temperature, pressure, mNaCl);
            calData << i << " " << temperature << " " << pressure << " " << mN2 << " " << N2solu << endl;
            cout << i  << " " << temperature << " " << pressure << " "  << mN2 << " " << N2solu << endl;
        }
    }
    else if (iChoice == 27)
    {
        double pressure = 10.0;
        double temperature = 298.15;
        double mNaCl = 1.0;

        ofstream evapCalResults("evaporiteProcess.txt");

        CO2SoluInBrine co2Evaperite;
        co2Evaperite.CO2evaporationCal(temperature, pressure, mNaCl);
        GasWaterEquilibria * evap = co2Evaperite.getEvaporite();
        evapCalResults << "No. " << "ph " << "massH2O " << "mCO2  " << "mNa  " ;
        
        int sizeOfResults = (int) evap->map_water.size();

        map<string, minerals> * minerName = &(evap->map_mineral.at(1));
        map<string, minerals> ::iterator it0 = minerName->begin();
        for (; it0 != minerName->end(); it0++)
        {
            evapCalResults << it0->first<<" ";
        }
        evapCalResults << endl;
        for (int i = 1; i < sizeOfResults+1; i++)
        {
            double ph= evap->map_water.at(i).ph;
            double massofWater = evap->map_water.at(i).mass_of_H2O;
            double mCO2 = evap->map_water.at(i).masterSpecies_molality.find("C(4)")->second;
            double mNa = evap->map_water.at(i).masterSpecies_molality.find("Na")->second;
//            double yH2O = evap->map_gas.at(i).moleNumber.find("H2O")->second;
            map<string, minerals> * miner = &(evap->map_mineral.at(i));
            evapCalResults << i << " " << ph << " " << massofWater << " " << mCO2 << " " << mNa << " " ;
            map<string, minerals>::iterator it = miner->begin();
            for (; it != miner->end(); it++)
            {
                evapCalResults << it->second.final_mole<<"  ";
            }
            evapCalResults << endl;
        }
        
    }
    else if (iChoice == 28)
    {
        MagnesiteSolubility mgsSolubility;
        mgsSolubility.initialize();
        double pressureBar = 1;
        double temperatureK = 37.5 + 273.15;
        double CO2moleFrac = 1;
        double mNaCl = 2.182;
        double magMolality = 0.0;
        double Mg2Molality = 0.0;
        double pH = 7.0;
        cout << endl;
        ofstream outPut("magSolu.txt");
        
        mgsSolubility.calcMagnesiteSolubility(temperatureK, pressureBar, mNaCl, CO2moleFrac, magMolality, Mg2Molality, pH);
//        mgsSolubility.calcMagnesiteSolubilityNoCO2(temperatureK, pressureBar, mNaCl, magMolality);
        cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << magMolality << "  " << Mg2Molality <<" "<< pH << endl;
        
        /*temperatureK = 0 + 273.15;
        while (temperatureK < 101+ 273.15)
        {
            mgsSolubility.calcMagnesiteSolubility(temperatureK, pressureBar, mNaCl, CO2moleFrac, magMolality, Mg2Molality, pH);
            
            outPut << temperatureK << " " << magMolality << endl;
            temperatureK += 2;
        }
        outPut.close();*/
    }
    else if (iChoice == 29)
    {
        Nesquehonite nesSolubility;
        nesSolubility.initialize();
        double pressureBar = 100;// *1.01325;
        double temperatureK = 100+ 273.15;
        double CO2moleFrac = 1;
        double mNaCl = 0.0;
        double nesMolality = 0.0;
        double Mg2Molality = 0.0;
        double pH = 7.0;

        double nesSolubility_measure = 0.1196;
        double LogK_up = -4.0;
        double LogK_down = -6.2;

        ofstream outPut("nesSolubility.txt");
        //nesSolubility.calcNesquehoniteSolubility(temperatureK, pressureBar, mNaCl, CO2moleFrac, nesMolality, Mg2Molality, pH);
        ////        mgsSolubility.calcMagnesiteSolubilityNoCO2(temperatureK, pressureBar, mNaCl, magMolality);
        //cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << nesMolality << "  " << Mg2Molality << " " << pH << endl;
        ///*double logK = nesSolubility.searchLogK(temperatureK, pressureBar, mNaCl, CO2moleFrac, nesSolubility_measure, LogK_up, LogK_down);
        //cout << "LogK = " << logK << endl;*/

        while (mNaCl <= 6)
        {
            nesSolubility.calcNesquehoniteSolubility(temperatureK, pressureBar, mNaCl, CO2moleFrac, nesMolality, Mg2Molality, pH);
            outPut << mNaCl << " " << nesMolality << endl;
            mNaCl += 0.02;
        }
    }
    else if (iChoice == 30)
    {
        Lansfordite lansSolubility;
        lansSolubility.initialize();
        double pressureBar = 100;//* 1.01325;
        double temperatureK = 100 + 273.15;
        double CO2moleFrac = 1;
        double mNaCl = 0.0;
        double lansMolality = 0.0;
        double Mg2Molality = 0.0;
        double pH = 7.0;

        double lansSolubility_measure = 0.402;
        double LogK_up = -4.0;
        double LogK_down = -6.2;
        ofstream outPut("lansforditeSolu.txt");
        //lansSolubility.calcLansforditeSolubility(temperatureK, pressureBar, mNaCl, CO2moleFrac, lansMolality, Mg2Molality, pH);
        ////        mgsSolubility.calcMagnesiteSolubilityNoCO2(temperatureK, pressureBar, mNaCl, magMolality);
        //cout << temperatureK << "  " << pressureBar << "  " << mNaCl << "  " << lansMolality << "  " << Mg2Molality << " " << pH << endl;
        ///*double logK = lansSolubility.searchLogK(temperatureK, pressureBar, mNaCl, CO2moleFrac, lansSolubility_measure, LogK_up, LogK_down);
        //cout << "LogK = " << logK << endl;*/

        while (mNaCl <= 6)
        {
            lansSolubility.calcLansforditeSolubility(temperatureK, pressureBar, mNaCl, CO2moleFrac, lansMolality, Mg2Molality, pH);
            
            outPut << mNaCl << " " << lansMolality << endl;
            mNaCl += 0.2;
        }
    }
    else if (iChoice == 31)
    {
        MgCO3Coexit mgco3;
        mgco3.initialize();
        double pressureBar = 1;
        double temperatureK = 437.15;
        double CO2moleFrac = 1.0;
        double mNaCl = 0.0;
        double MgCO3molality = 0.0;
        double Mg2molality = 0.0;
        double pH = 7.0;
        vector<bool> minralIndex;
        ofstream outPut("phaseDiagram.txt");

        while (pressureBar < 1000)
        {
            temperatureK = 273.15;
            while (temperatureK < 480.15)
            {

                mgco3.calcMgCO3Equilibria(temperatureK, pressureBar, mNaCl, CO2moleFrac, MgCO3molality, Mg2molality, pH, minralIndex);
                outPut << pressureBar << " " << temperatureK << " " << minralIndex[0] << " " << minralIndex[1] << " " << minralIndex[2] << endl;
                temperatureK += 2;
            }
            outPut << endl;
            pressureBar += 5;
        }
    }
    else if (iChoice == 32)
    {
        MgCO3Coexit mgCO3;
        mgCO3.initialize();
        double temperatureK = 100 + 273.15;
        double pressureBar = 300.0;
        double CO2moleFrac = 1.0;
        double mNaCl = 2.0;
        double pH = 7.0;
        double mgCO3molality = 0;
        double mg2molality = 0.0;
        vector<bool> mineralIndex;
        mgCO3.calcMgCO3EquilibriaNoCO2(temperatureK, pressureBar, mNaCl, mgCO3molality);

        ofstream outPut("initialSolution.txt");
        ofstream outAfterCO2inj("solutionAfterCO2inj.txt");

        //outPut << endl;
        waterSolution_results *water = &mgCO3.m_water;
        map<string, double>::iterator it = water->species_molality.begin();
        outPut << "PH: " << water->ph<<endl;
        for (; it != water->species_molality.end(); it++)
        {
            outPut << it->first << " " << it->second << endl;
        }

        mgCO3.calcMgCO3Equilibria(temperatureK, pressureBar, mNaCl, CO2moleFrac, mgCO3molality, mg2molality, pH, mineralIndex);

        water = &mgCO3.m_water;
        it = water->species_molality.begin();
        outAfterCO2inj << "PH: " << water->ph<<endl;
        for (; it != water->species_molality.end(); it++)
        {
            outAfterCO2inj << it->first << " " << it->second << endl;
        }
    }
    else if (iChoice == 33)
    {
        CO2SoluInBrine co2solubility;
        double temperatureK = 453.15;
        double mNaCl = 4.0;
        double pressureBar = 1.0;
        ofstream outPutCO2solubility("CO2solu_MgCO3system.txt");
        while (pressureBar < 1000)
        {
            double CO2solu = co2solubility.calCO2Solubility(temperatureK, pressureBar, mNaCl, 60);
            outPutCO2solubility << temperatureK << " " << pressureBar <<"  "<< mNaCl << " " << CO2solu << endl;
       //     cout << temperatureK << " " << pressureBar << "  " << mNaCl << " " << CO2solu << endl;
            pressureBar += 5;
        }
    }
    else if (iChoice == 34)
    {
        KineticsExample example6;
        example6.initialize();
        example6.m_getResults();
        
        waterSolution_results water;
        water = example6.water;
        conditionChange newCondition;
        newCondition.deltaTime = 10519200;
        newCondition.mass_of_H2O = 1.0;
        newCondition.gasPhaseIn = false;
        newCondition.pressure = 1.01;
        newCondition.temperature = 25.0;
        newCondition.masterChange = water.masterSpecies_molality;
        for (int i = 0; i < 15; i++)
        {
            example6.m_stepRun(newCondition);
            example6.m_getResults();
            water = example6.water;
            newCondition.masterChange = water.masterSpecies_molality;
            double  results = water.masterSpecies_molality.find("Si")->second*1000;
            cout << (i + 1)*newCondition.deltaTime / 3.1536e7 << "  " << results << endl;
        }
    }
    else if (iChoice == 35)
    {
        double time = clock();
        QuestBatch test;
        double pressureBar = 1.0;
        LDBLE atm_1bar = 0.9869233;
        conditionChange condition;
        condition.temperature = 61.0;
        condition.pressure = 1.0;  // pressureBar*atm_1bar;
        condition.mass_of_H2O = 1.0;
        condition.deltaTime = 0;
   //     test.initialize(condition);
        map<string, double> gasMoleNumber;
        gasMoleNumber["CO2"] = 100.0;
        ofstream outputs("Quest_standalone_fug_98bar_corr6_1000yrs.out");
        test.initialize(condition, outputs);
        test.questBatchModel("Quest_standalone_fug_98bar_corr6_1000yrs.phrq", gasMoleNumber);
        double time_seconds = (clock() - time) / CLOCKS_PER_SEC;
        cout << "Time used: " << time_seconds << endl;
    }
    else if (iChoice == 36)
    {
        double pressureBar = 100.0;
        double temperature = 50.0+273.15;
        double gasMoleNumber = 10.0;
        H2S_solubility h2sSolu;
        double h2sMolality = h2sSolu.calH2SSolubilityInBrineNewMethod(temperature, pressureBar, 0.0);
        cout << h2sMolality << endl;
    }
    else
    {
        N2SolubilityInBrine N2solubility;
        N2solubility.initialize();
        double N2solu = N2solubility.calN2Solubility(temperatureK, pressureBar, mNaCl);
        cout << "N2 solubility in molality (mole/KgW): " << N2solu << " !" << endl;
        cout << "H2O solubility in N2 rich phase: " << N2solubility.m_yH2O << " !" << endl;
    }
#endif
}
