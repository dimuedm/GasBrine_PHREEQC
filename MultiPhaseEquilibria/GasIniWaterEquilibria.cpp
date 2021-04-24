
#include "../logError.h"
#include "CommonFunctions.h"
#include "GasIniWaterEquilibria.h"
GasIniWaterEquilibria::GasIniWaterEquilibria()
{
}

GasIniWaterEquilibria::~GasIniWaterEquilibria()
{
}

void GasIniWaterEquilibria::initializePhreeqc(string inputFileName)
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

  //  phreeqc.set_kinetics_time(-2, 0.0);

    phreeqc.set_ReSoC_intitial(true);
    
    phreeqc.run_simulations();
    
    phreeqc.Get_phrq_io()->clear_istream();
    phreeqc.gasComponentDatabase.initialize();
    phreeqc.initialReSoC = false;
    phreeqc.results_for_ReSoC(finalGas, finalWater, finalMineral, finalKineticMineral);
}

void GasIniWaterEquilibria::updateReaction(conditionChange &newCondition)
{
    phreeqc.run_simulation_timeStep(newCondition);
    phreeqc.results_for_ReSoC(finalGas, finalWater, finalMineral, finalKineticMineral);
}

void GasIniWaterEquilibria::generateGasCompName(map<string, double> gasMoleNumber)
{
    int numberOfGasComponents = gasMoleNumber.size();
    if (numberOfGasComponents <= 0)
    {
        cerr << "ERROR: no inputs of gas compnents! Check the inputs!" << endl;
    }

    int iGasComp = 0;
    bool H2Oin = false;
    map<string, double>::iterator it_Gas = gasMoleNumber.begin();
    string nameH2O = "H2O";
    for (;it_Gas != gasMoleNumber.end();it_Gas++)
    {
        if (it_Gas->first != nameH2O)
        {
            gasCompName.push_back(it_Gas->first);
            gasSequence[it_Gas->first] = iGasComp;
            iGasComp++;
        }
        else
        {
            gasCompName.push_back(nameH2O);
            gasSequence[nameH2O] = iGasComp;
            iGasComp++;
            H2Oin = true;
        }
    }

    if (!H2Oin)
    {
        gasCompName.push_back("H2O");
        gasSequence[nameH2O] = numberOfGasComponents;
    }
}

double GasIniWaterEquilibria::calPhaseEquilibria_RealH2Omass(map<string, double>& gasMoleNumber, double &a_massOfH2OInWater, map<string, double>& masterSpecies, double temperatureK, double pressureBar, double& totalMoleNumber)
{
    int numberOfGasComponents = (int)gasMoleNumber.size();
    vector<double> zFeed(numberOfGasComponents);
 //   gasCompName.resize(numberOfGasComponents);

    CommonFunctions RReqnSolver;

    const double moleWeightH2O = 18.105; //  g/mole
    double totalMole = 0.0;
    bool H2OFound = false;
    string nameH2O = "";
    double moleH2O = 0.0;
    double KH2O = KEquiliH2O(temperatureK, pressureBar);

    map<string, double> ::iterator it = gasMoleNumber.begin();
    for (; it != gasMoleNumber.end(); it++)
    {
        if (it->second < 0)
        {
            logError logerror;
            logerror.LOGERROR("GasWaterEquilibria: Negative mole fraction!");
            return 0;
        }
        if ((it->first == "H2O") || (it->first == "h2o") || (it->first == "H2o") || (it->first == "h2O"))
        {
            if (H2OFound)
            {
                logError logerror;
                logerror.LOGERROR("GasWaterEquilibria: More than 1 H2O component information are found!");
                return 0;
            }
            nameH2O = it->first;
            moleH2O = it->second;
            H2OFound = true;
        }

        totalMole += it->second;
    }

    if (totalMole <= 1.0e-12)
    {
        cerr << "ERROR: wrong defination of the problem!" << endl;
        exit(-1);
    }

    int iComponent = 0;
    it = gasMoleNumber.begin();
    for (; it != gasMoleNumber.end(); it++)
    {
        if (it->first != nameH2O)
        {
            zFeed[iComponent] = it->second / totalMole;
//            gasCompName[iComponent] = it->first;
            iComponent++;
        }
        else
        {
            zFeed[numberOfGasComponents - 1] = it->second / totalMole;
 //           gasCompName[numberOfGasComponents - 1] = nameH2O;
        }
    }

    if (!H2OFound)
    {
        zFeed.push_back(0.0);
        nameH2O = "H2O";
        
        gasMoleNumber[nameH2O] = 0.0;
        numberOfGasComponents++;
        moleH2O = 0.0;
    }

    double massOfH2O = a_massOfH2OInWater;
    double moleNumberOfH2O = massOfH2O * 1000.0 / moleWeightH2O;

    for (int i = 0; i < numberOfGasComponents - 1; i++)
    {
        zFeed[i] *= totalMole / (totalMole + moleNumberOfH2O);
    }
    zFeed[numberOfGasComponents - 1] = (gasMoleNumber.find(nameH2O)->second + moleNumberOfH2O) / (totalMole + moleNumberOfH2O);

    totalMoleNumber = totalMole + moleNumberOfH2O;

    LDBLE atm_1bar = 0.9869233;
    phreeqc.gasComponentDatabase.initialize();
    conditionChange  solutionCondition;
    solutionCondition.pressure = pressureBar;
    solutionCondition.pressure *= atm_1bar;
    solutionCondition.temperature = temperatureK - 273.15;
    solutionCondition.mass_of_H2O = massOfH2O;
    solutionCondition.masterChange = masterSpecies;

    int numberOfgasComposition_H2OFree = numberOfGasComponents - 1;
    iterations = 0;

    double gasMoleFraction;
    xAqueous.resize(numberOfGasComponents);
    yGas.resize(numberOfGasComponents);
    vector<double> gasCompMolality;
    gasCompMolality.resize(numberOfgasComposition_H2OFree);
    elemCoeff.resize(numberOfgasComposition_H2OFree);
    elemName.resize(numberOfgasComposition_H2OFree);
    for (int i = 0; i < numberOfgasComposition_H2OFree; i++)
    {
        xAqueous[i] = 0.01 * zFeed[i];
        gasCompMolality[i] = xAqueous[i] *1000.0/moleWeightH2O;

        elemCoeff[i] = getElemCoeff(gasCompName[i]);
        elemName[i] = getElemName(gasCompName[i]);
//        gasSequence[gasCompName[i]] = i;
    }
    xAqueous[numberOfGasComponents - 1] = 1.0 - 0.01 * (zFeed[numberOfGasComponents - 1]);

    for (int i = 0; i < numberOfgasComposition_H2OFree; i++)
    {
        string l_gasCompName = gasCompName[i];
        map<string, double > ::iterator it = solutionCondition.masterChange.find(l_gasCompName);
        if (it == solutionCondition.masterChange.end())
        {
            solutionCondition.masterChange[elemName[i]] = gasCompMolality[i] * elemCoeff[i];
        }
        else
        {
            it->second = gasCompMolality[i] * elemCoeff[i];
        }
    }

    double yH2OinitiEst = waterSaturationPressure(temperatureK) / pressureBar;
    for (int i = 0; i < numberOfgasComposition_H2OFree; i++)
    {
        yGas[i] = zFeed[i] / (1.0 - zFeed[numberOfgasComposition_H2OFree]) * (1.0 - yH2OinitiEst);
    }
    yGas[numberOfgasComposition_H2OFree] = yH2OinitiEst;

    vector<double> equiliConst;
    equiliConst.resize(numberOfGasComponents);
    vector <double> Kvalue;
    Kvalue.resize(numberOfGasComponents);

    map<string, speciesComp> gasCompMap;
    vector<speciesComp> gasCompVector;
    gasCompVector.resize(numberOfGasComponents - 1);
    for (int i = 0; i < numberOfGasComponents - 1; i++)
    {
        gasCompMap.insert(pair<string, speciesComp>(gasCompName[i], gasCompVector[i]));
    }

    if (numberOfGasComponents == 2)
    {
        EOSPhaseProp eosProp;
        eosProp.initialize(gasCompName, yGas, &phreeqc.gasComponentDatabase, GAS, temperatureK, pressureBar, PR78);
        eosProp.calcFugacity();

        int gasCompIndex = getGasCompIndex(gasCompName[0]);
        equiliConst[0] = calEquiliConst(gasCompIndex, temperatureK, pressureBar);

        phreeqc.calActivityCoeff(solutionCondition, gasCompMap);
        double KGas = gasCompMap[gasCompName[0]].activityCoeff * equiliConst[0] / exp(eosProp.logPhi[0]) / pressureBar;
        double KWater = KH2O / exp(eosProp.logPhi[1]) / pressureBar;

        yGas[1] = (1.0 - KGas) / (1.0 - KGas / KWater);
        yGas[0] = 1.0 - yGas[1];

        xAqueous[0] = (1.0 - 1.0 / KWater) / (1.0 - KGas / KWater);
        xAqueous[1] = 1.0 - xAqueous[0];

        if (zFeed[0] <= xAqueous[0]) //single aqueous phase!
        {
            gasMoleFraction = 0.0;
            xAqueous[0] = zFeed[0];
            xAqueous[1] = zFeed[1];
        }
        else if (zFeed[1] <= yGas[1]) //single gas phase!
        {
            gasMoleFraction = 1.0;
            yGas[0] = zFeed[0];
            yGas[1] = zFeed[1];
        }
        else
        {
            gasMoleFraction = (zFeed[1] - xAqueous[1]) / (yGas[1] - xAqueous[1]);
        }
        return gasMoleFraction;
    }
    else if (numberOfGasComponents > 2)
    {
        vector <double> Kiold(numberOfGasComponents);
        for (int i = 0; i < numberOfGasComponents; i++)
        {
            Kiold[i] = 100.0;

            if (i < numberOfgasComposition_H2OFree)
            {
                int gasIndex = getGasCompIndex(gasCompName[i]);
                equiliConst[i] = calEquiliConst(gasIndex, temperatureK, pressureBar);
            }
        }
        double gasMole = 0.5;
        double maxDiffKvalue = 0.0;

        iterations = 0;
        while (true)
        {
            iterations++;

            if (iterations >= 100)
            {
                logError logerror;
                logerror.LOGERROR("GasWaterEquilibria: Too many iterations!");
                return -1;
            }

            EOSPhaseProp eosProp;

            eosProp.initialize(gasCompName, yGas, &phreeqc.gasComponentDatabase, GAS, temperatureK, pressureBar, PR78);
            eosProp.calcFugacity();
            phreeqc.calActivityCoeff(solutionCondition, gasCompMap);

            for (int i = 0; i < numberOfgasComposition_H2OFree; i++)
            {
                Kvalue[i] = equiliConst[i] * gasCompMap[gasCompName[i]].activityCoeff / exp(eosProp.logPhi[i]) / pressureBar;
            }
            Kvalue[numberOfgasComposition_H2OFree] = KH2O / exp(eosProp.logPhi[numberOfgasComposition_H2OFree]) / pressureBar;

            maxDiffKvalue = 0.0;
            for (int i = 0; i < numberOfGasComponents; i++)
            {
                if (fabs((Kvalue[i] - Kiold[i]) / Kiold[i]) > maxDiffKvalue)
                {
                    maxDiffKvalue = fabs((Kvalue[i] - Kiold[i]) / Kiold[i]);
                }
            }
            if (maxDiffKvalue < 1.0e-5)
            {
                return gasMole;
            }

            RReqnSolver.calcGasMole(gasMole, zFeed, Kvalue);
            if (gasMole < 0 || gasMole>1.0)
            {
                return gasMole;
            }

            for (int i = 0; i < numberOfGasComponents; i++)
            {
                xAqueous[i] = zFeed[i] / (1.0 - gasMole + gasMole * Kvalue[i]);
                yGas[i] = Kvalue[i] * xAqueous[i];
            }
            for (int i = 0; i < numberOfGasComponents; i++)
            {
                Kiold[i] = Kvalue[i];
            }
        }
        return gasMole;
    }
    else
    {
        logError logerror;
        logerror.LOGERROR("Error: Component number!");
        return -1;
    }
}

double GasIniWaterEquilibria::calPhaseEquilibria(map<string, double> &gasMoleNumber, map<string, double> &masterSpecies, double temperatureK, double pressureBar)
{
 //  struct species * dissolvedGasComp = 

    int numberOfGasComponents = (int) gasMoleNumber.size();
    vector<double> zFeed(numberOfGasComponents);
    gasCompName.resize(numberOfGasComponents);

    CommonFunctions RReqnSolver;

    const double moleWeightH2O = 18.105; //  g/mole
    double totalMole = 0.0;
    bool H2OFound = false;
    string nameH2O = "";
    double moleH2O = 0.0;
    double KH2O = KEquiliH2O(temperatureK, pressureBar);

    map<string, double> ::iterator it = gasMoleNumber.begin();
    for (; it != gasMoleNumber.end(); it++)
    {
        if (it->second<0)
        {
            logError logerror;
            logerror.LOGERROR("GasWaterEquilibria: Negative mole fraction!");
            return 0;
        }
        if ((it->first == "H2O") || (it->first == "h2o") || (it->first == "H2o") || (it->first == "h2O"))
        {
            if (H2OFound)
            {
                logError logerror;
                logerror.LOGERROR("GasWaterEquilibria: More than 1 H2O component information are found!");
                return 0;
            }
            nameH2O = it->first;
            moleH2O = it->second;
            H2OFound = true;
        }

        totalMole += it->second;
    }

    int iComponent = 0;
    it = gasMoleNumber.begin();
    for (; it != gasMoleNumber.end(); it++)
    {
        if (it->first != nameH2O)
        {
            zFeed[iComponent] = it->second / totalMole;
            gasCompName[iComponent] = it->first;
            iComponent++;
        }
        else
        {
            zFeed[numberOfGasComponents - 1] = it->second / totalMole;
            gasCompName[numberOfGasComponents - 1] = nameH2O;
        }
    }

    if (!H2OFound)
    {
        zFeed.push_back(0.0);
        nameH2O = "H2O";
        gasCompName.push_back("H2O");
        gasMoleNumber[nameH2O] = 0.0;
        numberOfGasComponents++;
        moleH2O = 0.0;
    }

    double Nw = 55.508;
    for (int i = 0; i < numberOfGasComponents - 1; i++)
    {
        zFeed[i] *= totalMole / (totalMole + Nw);
    }
    zFeed[numberOfGasComponents - 1] = (gasMoleNumber[nameH2O] + Nw) / (totalMole + Nw);

    LDBLE atm_1bar = 0.9869233;
    phreeqc.gasComponentDatabase.initialize();
    conditionChange  solutionCondition;
    solutionCondition.pressure = pressureBar;
    solutionCondition.pressure *= atm_1bar;
    solutionCondition.temperature = temperatureK - 273.15;
    solutionCondition.mass_of_H2O = 1.0;

    /*update_timeStep.gasPhaseIn = TRUE;
    map <string, double> gasPhaseMoleNumber_update;
    gasPhaseMoleNumber_update = gasPhaseMoleNumber;
    gasPhaseMoleNumber_update.erase(nameH2O);*/
    double moleNumberOfH2OinGas_old = moleH2O;
    solutionCondition.masterChange = masterSpecies;

    int numberOfgasComposition_H2OFree = numberOfGasComponents - 1;
    iterations = 0;

    double gasMoleFraction;
    xAqueous.resize(numberOfGasComponents);
    yGas.resize(numberOfGasComponents);
    vector<double> gasCompMolality;
    gasCompMolality.resize(numberOfgasComposition_H2OFree);
    elemCoeff.resize(numberOfgasComposition_H2OFree);
    elemName.resize(numberOfgasComposition_H2OFree);
    for (int i = 0; i < numberOfgasComposition_H2OFree; i++)
    {
        xAqueous[i] = 0.01*zFeed[i];
        gasCompMolality[i] = xAqueous[i] * Nw;

        elemCoeff[i] = getElemCoeff(gasCompName[i]);
        elemName[i] = getElemName(gasCompName[i]);
    }
    xAqueous[numberOfGasComponents - 1] = 1.0 - 0.01*(zFeed[numberOfGasComponents - 1]);

    for (int i = 0; i < numberOfgasComposition_H2OFree; i++)
    {
        string l_gasCompName = gasCompName[i];
        map<string, double > ::iterator it = solutionCondition.masterChange.find(l_gasCompName);
        if (it == solutionCondition.masterChange.end())
        {
            solutionCondition.masterChange[elemName[i]] = gasCompMolality[i]*elemCoeff[i];
        }
        else
        {
            it->second = gasCompMolality[i] * elemCoeff[i];
        }
    }

    double yH2OinitiEst = waterSaturationPressure(temperatureK) / pressureBar;
    for (int i = 0; i < numberOfgasComposition_H2OFree; i++)
    {
        yGas[i] = zFeed[i] / (1.0 - zFeed[numberOfgasComposition_H2OFree])*(1.0 - yH2OinitiEst);
    }
    yGas[numberOfgasComposition_H2OFree] = yH2OinitiEst;

    vector<double> equiliConst;
    equiliConst.resize(numberOfGasComponents); 
    vector <double> Kvalue;
    Kvalue.resize(numberOfGasComponents);

    map<string, speciesComp> gasCompMap;
    vector<speciesComp> gasCompVector;
    gasCompVector.resize(numberOfGasComponents - 1);
    for (int i = 0; i < numberOfGasComponents - 1; i++)
    {
        gasCompMap.insert(pair<string, speciesComp>(gasCompName[i], gasCompVector[i]));
    }

    if (numberOfGasComponents == 2)
    {
        EOSPhaseProp eosProp;
        eosProp.initialize(gasCompName, yGas, &phreeqc.gasComponentDatabase, GAS, temperatureK, pressureBar, PR78);
        eosProp.calcFugacity();

        int gasCompIndex = getGasCompIndex(gasCompName[0]);
        equiliConst[0] = calEquiliConst(gasCompIndex, temperatureK, pressureBar);

        phreeqc.calActivityCoeff(solutionCondition, gasCompMap);
        double KGas = gasCompMap[gasCompName[0]].activityCoeff*equiliConst[0]/exp(eosProp.logPhi[0])/pressureBar;
        double KWater = KH2O / exp(eosProp.logPhi[1]) / pressureBar;

        yGas[1] = (1.0-KGas) / (1.0-KGas/KWater);
        yGas[0] = 1.0 - yGas[1];

        xAqueous[0] = (1.0-1.0/KWater) / (1.0-KGas/KWater);
        xAqueous[1] = 1.0 - xAqueous[0];

        if (zFeed[0] <= xAqueous[0]) //single aqueous phase!
        {
            gasMoleFraction = 0.0;
            xAqueous[0] = zFeed[0];
            xAqueous[1] = zFeed[1];
        }
        else if (zFeed[1] <= yGas[1]) //single gas phase!
        {
            gasMoleFraction = 1.0;
            yGas[0] = zFeed[0];
            yGas[1] = zFeed[1];
        }
        else
        {
            gasMoleFraction = (zFeed[1] - xAqueous[1]) / (yGas[1] - xAqueous[1]);
        }
        return gasMoleFraction;
    }
    else if (numberOfGasComponents>2)
    {
        vector <double> Kiold(numberOfGasComponents);
        for (int i = 0; i < numberOfGasComponents; i++)
        {
            Kiold[i] = 100.0;
        }
        double gasMole = 0.5;
        double maxDiffKvalue = 0.0;
        iterations = 0;
        while (true)
        {
            iterations++;

            if (iterations >= 100)
            {
                logError logerror;
                logerror.LOGERROR("GasWaterEquilibria: Too many iterations!");
                return -1;
            }

            EOSPhaseProp eosProp;

            eosProp.initialize(gasCompName, zFeed, &phreeqc.gasComponentDatabase, GAS, temperatureK, pressureBar, PR78);
            eosProp.calcFugacity();
            phreeqc.calActivityCoeff(solutionCondition, gasCompMap);

            for (int i = 0; i < numberOfgasComposition_H2OFree; i++)
            {
                Kvalue[i] = equiliConst[i] * gasCompMap[gasCompName[i]].activityCoeff / exp(eosProp.logPhi[i]) / pressureBar;
            }
            Kvalue[numberOfgasComposition_H2OFree] = KH2O / exp(eosProp.logPhi[numberOfgasComposition_H2OFree]) / pressureBar;

            maxDiffKvalue = 0.0;
            for (int i = 0; i < numberOfGasComponents; i++)
            {
                if (fabs((Kvalue[i] - Kiold[i]) / Kiold[i])>maxDiffKvalue)
                {
                    maxDiffKvalue = fabs((Kvalue[i] - Kiold[i]) / Kiold[i]);
                }
            }
            if (maxDiffKvalue < 1.0e-5)
            {
                return gasMole;
            }

            RReqnSolver.calcGasMole(gasMole, zFeed, Kvalue);
            if (gasMole<0 || gasMole>1.0)
            {
                return gasMole;
            }

            for (int i = 0; i < numberOfGasComponents; i++)
            {
                xAqueous[i] = zFeed[i] / (1.0 - gasMole + gasMole*Kvalue[i]);
                yGas[i] = Kvalue[i] * xAqueous[i];
            }
            for (int i = 0; i < numberOfGasComponents; i++)
            {
                Kiold[i] = Kvalue[i];
            }
        }
        return gasMole;
    }
    else
    {
        logError logerror;
        logerror.LOGERROR("Error: Component number!");
        return -1;
    }
}

double GasIniWaterEquilibria::KEquiliH2O(double temperature, double pressure)
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

double GasIniWaterEquilibria::waterSaturationPressure(double a_temperature)
{
    double a1, a2, a3, a4, a5, a6;
    a1 = -7.85951783;
    a2 = 1.84408259;
    a3 = -11.7866497;
    a4 = 22.6807411;
    a5 = -15.9618719;
    a6 = 1.80122502;

    double Pc, Tc, tau;
    Pc = 221.2;
    Tc = 647.3;
    tau = 1.0 - a_temperature / Tc;

    double ln_Ps_Pc = Tc / a_temperature*(a1*tau + a2*pow(tau, 1.5) + a3*tau*tau*tau + a4*pow(tau, 3.5) + a5*pow(tau, 4) + a6*pow(tau, 7.5));
    double Ps = Pc*exp(ln_Ps_Pc);
    return Ps;
}

double GasIniWaterEquilibria::calEquiliConst(int NO_forWaterSolu, double temperature, double pressure)
{
    switch (NO_forWaterSolu)
    {
    case 1:
        return calWaterEquiliConst(temperature, pressure);
        break;
    case 2:
        return calCO2EquiliConst(temperature, pressure);
        break;
    case 3:
        return calCH4EquiliConst(temperature, pressure);
        break;
    case 4:
        return calN2EquiliConst(temperature, pressure);
        break;
    case 5:
        return calH2SEquiliConst(temperature, pressure);
        break;
    case 6:
        return calO2EquiliConst(temperature, pressure);
        break;
    default:
        return bigValue;
        break;
    }
}

int GasIniWaterEquilibria::getGasCompIndex(string gasCompName)
{
    if (gasCompName == "CO2" || gasCompName == "co2" || gasCompName == "Co2" || gasCompName == "cO2")
    {
        return 2;
    }
    else if (gasCompName == "CH4" || gasCompName == "ch4" || gasCompName == "cH4" || gasCompName == "Ch4")
    {
        return 3;
    }
    else if (gasCompName == "N2" || gasCompName == "n2")
    {
        return 4;
    }
    else if (gasCompName == "H2S" || gasCompName == "h2s" || gasCompName == "H2s" || gasCompName == "h2S"||gasCompName=="H2Sg")
    {
        return 5;
    }
    else if (gasCompName == "O2" || gasCompName == "o2")
    {
        return 6;
    }
    else
    {
        return 0;
    }
}

double GasIniWaterEquilibria::calWaterEquiliConst(double temperature, double pressure)
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

double  GasIniWaterEquilibria::calCH4EquiliConst(double temperature, double pressure)
{
    double a1 = -16.3978647450248;
    double a2 = 0.0325711411111918;
    double a3 = 9468.16073642105;
    double a4 = -2.65734173526389e-5;
    double a5 = -1435285.11188553;
    double a6 = -0.0132253629347427;
    double a7 = 3.13045055468136e-5;
    double a8 = 2.26047549807049;
    double a9 = -2.90633722703171e-8;
    double a10 = 3.44079701140235e-9;
    double a11 = -6.11808900183315e-10;
    double equiliConstValue = a1 + a2*temperature + a3 / temperature + a4*temperature*temperature + a5 / temperature / temperature + a6*pressure + a7*pressure*temperature + a8*pressure / temperature + a9*pressure*temperature*temperature + a10*pressure*pressure * temperature + a11*pressure*pressure*pressure;

    equiliConstValue = exp(equiliConstValue)*55.508;
    return equiliConstValue;
}

double  GasIniWaterEquilibria::calCO2EquiliConst(double temperature, double pressure)
{
    double c1 = 0.2301854e2;
    double c2 = -0.36540569e-1;
    double c3 = -0.18366895e4;
    double c4 = 0.20330876e-4;
    double c5 = -0.39072384e6;
    double c6 = -0.58269326e-1;
    double c7 = 0.15061716e-3;
    double c8 = 0.78086969e1;
    double c9 = -0.13013307e-6;
    double c10 = 0.11145375e-8;
    double c11 = -0.13073985e-9;

    double d1 = 0.11625439e4;
    double d2 = -0.1335038e1;
    double d3 = -0.44522815e6;
    double d4 = 0.5732654e-3;
    double d5 = 0.64139318e8;
    double d6 = -0.12735549e0;
    double d7 = 0.21720184e-3;
    double d8 = 0.25529557e2;
    double d9 = -0.12542533e-6;
    double d10 = 0.0;
    double d11 = 0.0;

    double equiliConstValue;

    if (temperature < 503.0)
    {
        equiliConstValue = c1 + c2*temperature + c3 / temperature + c4*temperature*temperature + c5 / temperature / temperature + c6*pressure + c7*pressure*temperature + c8*pressure / temperature + c9*pressure*temperature*temperature + c10*pressure*pressure * temperature + c11*pressure*pressure*pressure;
    }
    else if (temperature > 523.0)
    {
        equiliConstValue = d1 + d2*temperature + d3 / temperature + d4*temperature*temperature + d5 / temperature / temperature + d6*pressure + d7*pressure*temperature + d8*pressure / temperature + d9*pressure*temperature*temperature + d10*pressure*pressure * temperature + d11*pressure*pressure*pressure;
    }
    else
    {
        double t1 = 503.0;
        double t2 = 523.0;
        double a503 = c1 + c2*t1 + c3 / t1 + c4*t1*t1 + c5 / t1 / t1 + c6*pressure + c7*pressure*t1 + c8*pressure / t1 + c9*pressure*t1*t1 + c10*pressure*pressure * t1 + c11*pressure*pressure*pressure;
        double a523 = d1 + d2*t2 + d3 / t2 + d4*t2*t2 + d5 / t2 / t2 + d6*pressure + d7*pressure*t2 + d8*pressure / t2 + d9*pressure*t2*t2 + d10*pressure*pressure * t2 + d11*pressure*pressure*pressure;
        equiliConstValue = a503 + (a523 - a503)*(temperature - 503.0) / 20.0;
    }
    equiliConstValue = exp(equiliConstValue);
    equiliConstValue *= 55.5080;
    return equiliConstValue;
}

double GasIniWaterEquilibria::calH2SEquiliConst(double temperature, double pressure)
{
    /*double c1 = -825.351676920731;
    double c2 = 1.46600846686058;
    double c3 = 207875.144655533;
    double c4 = -0.000966632169841103;
    double c5 = -19617368.7486012;
    double c6 = 0.823684902091498;
    double c7 = -0.00222376062074741;
    double c8 = -101.113160193415;
    double c9 = 1.96421386690219e-6;
    double c10 = 3.77937696024575e-8;
    double c11 = 2.92076340659316e-9;
    double equiliConstValue = c1 + c2*temperature + c3 / temperature + c4*temperature*temperature + c5 / temperature / temperature + c6*pressure + c7*pressure*temperature + c8*pressure / temperature + c9*pressure*temperature*temperature + c10*pressure*pressure * temperature + c11*pressure*pressure*pressure;
   */
  //  double c1 = 42.564957;
  //  double c2 = -0.086260377;
  //  double c3 = -6084.3775;
  //  double c4 = 0.000068714437;
  //  double c5 = -102.76849;
  //  double c6 = 0.00084482895;
  //  double c7 = -1.0590768;
  //  double c8 = 0.0035665902;

  ////  double equiliConstValue = c1 + c2 * temperature + c3 / temperature + c4 * temperature * temperature + c5 / (680 - temperature) + c6 * pressure + c7 * pressure / (680 - temperature) + c8 * pressure * pressure / temperature;
  //  double equiliConstValue = c1 + c2 * temperature + c3 / temperature + c4 * temperature * temperature + c5 / (680 - temperature) + c6  + c7  / (680 - temperature) + c8  / temperature;
  //  equiliConstValue += 41.84 * 0.02 * (pressure - 1) / R / temperature;

    double c1 = 2.29279459095917E+01;
    double c2 = -3.01815918516077E-02;
    double c3 = -3.86194515547495E+03;
    double c4 = 1.28343525311545E-05;
    double c5 = -5.07658723195908E-03;
    double c6 = 3.00740134449227E+00;
    double c7 = 2.04029590002629E+00;
    double c8 = -1.71227175529834E-03;
    double c9 = 5.21402829195479E-06;
    double equiliConstValue = c1 + c2 * temperature + c3 / temperature + c4 * temperature * temperature + c5 * pressure + c6 / pressure + c7 * pressure / temperature + c8 * pressure * pressure / temperature + c9 * pressure * pressure;

    equiliConstValue = exp(equiliConstValue)*55.508;
    return equiliConstValue;
}

double  GasIniWaterEquilibria::calN2EquiliConst(double temperature, double pressure)
{
    /*double c1 = -23.093813e0;
    double c2 = 0.56048525e-1;
    double c3 = 0.98808898e4;
    double c4 = -0.51091621e-4;
    double c5 = -0.13220298e7;
    double c6 = -0.49542866e-3;
    double c7 = 0.12698747e-5;
    double c8 = 0.51411144e0;
    double c9 = -0.64733978e-4;
    double equiliConstValue = c1 + c2*temperature + c3 / temperature + c4*temperature*temperature + c5 / temperature / temperature + c6*pressure + c7*pressure*temperature + c8*pressure / temperature + c9*pressure*pressure / temperature;
   */

    double equiliConstValue = -58.453 + 1.818e-3 * temperature + 3199 / temperature + 17.909 * log10(temperature) - 27460 / temperature / temperature;
    equiliConstValue *= -2.3;
    equiliConstValue += 41.84 * 0.073 * (pressure - 1) / R / temperature;
    equiliConstValue = exp(equiliConstValue)*55.508 ;
    return equiliConstValue;
}

double  GasIniWaterEquilibria::calO2EquiliConst(double temperature, double pressure)
{
  //  (void)temperature; (void)pressure;

    double equiliConstValue = -7.5001 + 7.8981e-003 * temperature + 2.0027e+005 / temperature / temperature;
    equiliConstValue *= -2.3;
    equiliConstValue += 41.84 * 0.05 * (pressure - 1) / R / temperature;
    equiliConstValue = exp(equiliConstValue) * 55.508;

    return equiliConstValue;
}

string GasIniWaterEquilibria::getElemName(string a_gasCompName)
{
    if (a_gasCompName == "CO2" || a_gasCompName == "co2" || a_gasCompName == "Co2" || a_gasCompName == "cO2")
    {
        if (finalWater.masterSpecies_molality.find("C(4)") != finalWater.masterSpecies_molality.end())
        {
            return "C(4)";
        }
        else if (finalWater.masterSpecies_molality.find("C(+4)") != finalWater.masterSpecies_molality.end())
        {
            return "C(+4)";
        }
        else if (finalWater.masterSpecies_molality.find("C") != finalWater.masterSpecies_molality.end())
        {
            return "C";
        }
        else
        {
            cerr << "ERROR: cannot find master species name of CO2! Please check the database!" << endl;
            exit(-1);
            return "";
        }
    }
    else if (a_gasCompName == "O2" || a_gasCompName == "o2")
    {
        return "O(0)";
    }
    else if (a_gasCompName == "N2" || a_gasCompName == "n2")
    {
        return "N(0)";
    }
    else if (a_gasCompName == "CH4" || a_gasCompName == "ch4" || a_gasCompName == "cH4" || a_gasCompName == "Ch4")
    {
        return "C(-4)";
    }
    else if (a_gasCompName == "H2S" || a_gasCompName == "h2s" || a_gasCompName == "H2s" || a_gasCompName == "h2S")
    {
        return "S(-2)";
    }
    else if (a_gasCompName == "H2Sg" || a_gasCompName == "h2sg" || a_gasCompName == "H2sg" || a_gasCompName == "h2Sg")
    {
        return "Sg";
    }
    else
    {
        logError logerror;
        logerror.LOGERROR("Unknown gas component!");
        return "";
    }
}

double GasIniWaterEquilibria::getElemCoeff(string a_gasCompName)
{
    if (a_gasCompName == "CO2" || a_gasCompName == "co2" || a_gasCompName == "Co2" || a_gasCompName == "cO2")
    {
        return 1.0;
    }
    else if (a_gasCompName == "O2" || a_gasCompName == "o2")
    {
        return 2.0;
    }
    else if (a_gasCompName == "N2" || a_gasCompName == "n2")
    {
        return 2.0;
    }
    else if (a_gasCompName == "CH4" || a_gasCompName == "ch4" || a_gasCompName == "cH4" || a_gasCompName == "Ch4")
    {
        return 1.0;
    }
    else if (a_gasCompName == "H2S" || a_gasCompName == "h2s" || a_gasCompName == "H2s" || a_gasCompName == "h2S")
    {
        return 1.0;
    }
    else if (a_gasCompName == "H2Sg" || a_gasCompName == "h2sg" || a_gasCompName == "H2sg" || a_gasCompName == "h2Sg")
    {
        return 1.0;
    }
    else
    {
        logError logerror;
        logerror.LOGERROR("Unknown gas component!");
        return -1.0;
    }
}