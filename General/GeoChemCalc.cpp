#include "GeoChemCalc.h"

GeoChemCalc::GeoChemCalc()
{
    m_gasIn_ = false;
}

GeoChemCalc::~GeoChemCalc()
{
}

//void GeoChemCalc::readIn(string a_inputFile)
//{
//    m_geochemInputs_.
//    m_geochemInputs_.readInputs()
//}

void GeoChemCalc::initilize(map<string, vector<string> >* a_gcmData)
{    
    m_geochemInputs_.readInputs(a_gcmData);
    
    string l_phrqcInputFile = *m_geochemInputs_.m_getPhrqcInputFile();
    m_gasWater_.initializePhreeqc(l_phrqcInputFile);

    m_inputWater_ = m_gasWater_.finalWater;
    
    map<string, minerals> ::iterator it_eqmin = m_gasWater_.finalMineral.begin();
    for (;it_eqmin != m_gasWater_.finalMineral.end();it_eqmin++)
    {
        string l_mineralName = it_eqmin->first;
        double l_moleNumber = it_eqmin->second.final_mole;
        m_inputEqMinerals_[l_mineralName] = l_moleNumber;
    }

    map<string, kineticPhases>::iterator it_kinMineral = m_gasWater_.finalKineticMineral.begin();
    for (;it_kinMineral != m_gasWater_.finalKineticMineral.end();it_kinMineral++) 
    {
        string l_minealName = it_kinMineral->first;
        double l_moleNumber = it_kinMineral->second.final_mole;
        m_inputKineticMinerals_[l_minealName] = l_moleNumber;

        map<string, double> ::iterator it_sameEqmineral = m_inputEqMinerals_.find(l_minealName);
        if (it_sameEqmineral != m_inputEqMinerals_.end())
        {
            m_inputEqMinerals_.erase(it_sameEqmineral);
        }
    }    

    map<string, double>* l_gasMoleFraction = m_geochemInputs_.m_getGasMoleFraction();
    int l_gasSize = (int)l_gasMoleFraction->size();
    int l_gasListSize = (int)m_geochemInputs_.m_getMoleFractionList()->size();
    double l_totalFractionSummation = 0.0;
    if (l_gasSize <= 0 && l_gasListSize<=0)
    {
        cout << "Warning: there is no gas inputs!" << endl;
        return;
    }

    if (l_gasSize > 0)
    {
        map<string, double>::iterator it_Gas = l_gasMoleFraction->begin();
        for (; it_Gas != l_gasMoleFraction->end(); it_Gas++)
        {
            if (it_Gas->second < 0)
            {
                cout << "Warning: negative gas mole fraction! Set to 0.0. Check the inputs!" << endl;
                it_Gas->second = 0.0;
            }
            l_totalFractionSummation += it_Gas->second;
        }

        if (l_totalFractionSummation <= 0)
        {
            cout << "Warning: gas mole fraction summation is 0.0. Check the inputs!" << endl;
            m_gasIn_ = false;
            return;
        }

        it_Gas = l_gasMoleFraction->begin();
        m_totalGasMoleNumber_ = m_geochemInputs_.m_getGasMoleNumber();
        if (m_totalGasMoleNumber_ < 0.0)
        {
            cout << "Warning: total gas mole number is 0.0. Check the inputs!" << endl;
            m_gasIn_ = false;
            return;
        }
        for (; it_Gas != l_gasMoleFraction->end(); it_Gas++)
        {
            double l_gasMoleNumber = it_Gas->second / l_totalFractionSummation;
            m_gasMoleNumber_[it_Gas->first] = l_gasMoleNumber * m_totalGasMoleNumber_;
        }
        m_inputGasMoleNumber_ = m_gasMoleNumber_;
        m_gasWater_.generateGasCompName(m_gasMoleNumber_);
    }

    m_gasIn_ = true;
    m_redox_ = m_geochemInputs_.m_getRedox();

}

void GeoChemCalc::calculateWithIniWater(string a_inputFileName)
{

    m_outPutFileName_ = a_inputFileName;
    
    if (m_geochemInputs_.m_isKinetics)
    {
        m_condition_.masterChange = m_gasWater_.finalWater.masterSpecies_molality;
        m_condition_.gasPhaseIn = false;
        m_condition_.mass_of_H2O = 1.0;
        m_massOfH2OInWater = 1.0;

 //       m_condition_.mineralsUpdate = m_inputEqMinerals_;

        kineticReaction();
        return;
    }
    else if (m_geochemInputs_.m_getTempPressList()->size() > 0)
    {
        //pair<double,double> l_temperatureK_pressureBar = *m_geochemInputs_.m_getTempPressList()->begin();
        list<pair<double, double>>* l_temperaure_pressureList = m_geochemInputs_.m_getTempPressList();
        m_outPutFileName_ += "_temperature_pressure.out";

        vector<string>* l_masterPrint = m_geochemInputs_.m_getMasterPrint();
        vector<string>* l_secondaryPrint = m_geochemInputs_.m_getSecondaryPrint();
        vector<string>* l_mineralPrint = m_geochemInputs_.m_getMineralPrint();
        vector<string>* l_waterMolaltiyPrint = m_geochemInputs_.m_waterMolalityPrint();
        vector<string>* l_gasMoleFractionPrint = m_geochemInputs_.m_gasMoleFractionPrint();
        vector<string>* l_gasMoleNumberPrint = m_geochemInputs_.m_gasMoleNumbers();

        m_outPut_.open(m_outPutFileName_.c_str());
        m_outPut_ << "Temperature(K) " << "Pressure(bar) ";

        for (int i = 0;i < l_masterPrint->size();i++)
        {
            m_outPut_ << l_masterPrint->at(i) << " ";
        }
        for (int i = 0;i < l_secondaryPrint->size();i++)
        {
            m_outPut_ << l_secondaryPrint->at(i) << " ";
        }
        for (int i = 0;i < l_mineralPrint->size();i++)
        {
            m_outPut_ << l_mineralPrint->at(i) << " ";
        }
        for (int i = 0;i < l_waterMolaltiyPrint->size();i++)
        {
            m_outPut_ << l_waterMolaltiyPrint->at(i) << " ";
        }
        /*for (int i = 0;i < l_gasMoleFractionPrint->size();i++)
        {
            m_outPut_ << l_gasMoleFractionPrint->at(i) << " ";
        }*/
        for (int i = 0;i < l_gasMoleNumberPrint->size();i++)
        {
            m_outPut_ << l_gasMoleNumberPrint->at(i) << " ";
        }

        m_outPut_ << endl;

        for (pair<double, double> l_temperatureK_pressureBar : *l_temperaure_pressureList)
        {
            double l_temperatureK = l_temperatureK_pressureBar.first;
            double l_pressureBar = l_temperatureK_pressureBar.second;
            m_calc_(l_temperatureK, l_pressureBar);
        }
        m_outPut_ << endl;
        return;
    }
    else if ( m_geochemInputs_.m_getMoleFractionList()->size()==0 && (m_geochemInputs_.m_getPressure()->size() > 1 || ((m_geochemInputs_.m_getPressure()->size() == 1) && (m_geochemInputs_.m_getTemperature()->size() == 1))))
    {
        double l_temperatureK = *m_geochemInputs_.m_getTemperature()->begin();
        list<double>* l_pressureList = m_geochemInputs_.m_getPressure();
        m_outPutFileName_ += "_pressure.out";

        vector<string>* l_masterPrint = m_geochemInputs_.m_getMasterPrint();
        vector<string>* l_secondaryPrint = m_geochemInputs_.m_getSecondaryPrint();
        vector<string>* l_mineralPrint = m_geochemInputs_.m_getMineralPrint();
        vector<string>* l_waterMolaltiyPrint = m_geochemInputs_.m_waterMolalityPrint();
        vector<string>* l_gasMoleFractionPrint = m_geochemInputs_.m_gasMoleFractionPrint();
        vector<string>* l_gasMoleNumberPrint = m_geochemInputs_.m_gasMoleNumbers();
        if (l_temperatureK < 273.15 || l_temperatureK >533.15)
        {
            cerr << "ERROR: Temperature out of range. Please check the inputs!" << endl;
            exit(-1);
        }

        m_outPut_.open(m_outPutFileName_.c_str());
        m_outPut_ << "Temperature(K) " << "Pressure(bar) ";

        for (int i = 0;i < l_masterPrint->size();i++)
        {
            m_outPut_ << l_masterPrint->at(i) << " ";
        }
        for (int i = 0;i < l_secondaryPrint->size();i++)
        {
            m_outPut_ << l_secondaryPrint->at(i) << " ";
        }
        for (int i = 0;i < l_mineralPrint->size();i++)
        {
            m_outPut_ << l_mineralPrint->at(i) << " ";
        }
        for (int i = 0;i < l_waterMolaltiyPrint->size();i++)
        {
            m_outPut_ << l_waterMolaltiyPrint->at(i) << " ";
        }
        /*for (int i = 0;i < l_gasMoleFractionPrint->size();i++)
        {
            m_outPut_ << l_gasMoleFractionPrint->at(i) << " ";
        }*/
        for (int i = 0;i < l_gasMoleNumberPrint->size();i++)
        {
            m_outPut_ << l_gasMoleNumberPrint->at(i) << " ";
        }

        m_outPut_ << endl;

        for (double l_pressureBar : *l_pressureList)
        {
            m_calc_(l_temperatureK, l_pressureBar);
        }
        m_outPut_ << endl;
        m_outPut_.close();
        return;
    }
    else if (m_geochemInputs_.m_getTemperature()->size() > 1)
    {
        double l_pressureBar = *m_geochemInputs_.m_getPressure()->begin();
        list<double>* l_temperatureList = m_geochemInputs_.m_getTemperature();
        m_outPutFileName_ += "temperature.txt";

        vector<string>* l_masterPrint = m_geochemInputs_.m_getMasterPrint();
        vector<string>* l_secondaryPrint = m_geochemInputs_.m_getSecondaryPrint();
        vector<string>* l_mineralPrint = m_geochemInputs_.m_getMineralPrint();
        vector<string>* l_waterMolaltiyPrint = m_geochemInputs_.m_waterMolalityPrint();
        vector<string>* l_gasMoleFractionPrint = m_geochemInputs_.m_gasMoleFractionPrint();
        vector<string>* l_gasMoleNumberPrint = m_geochemInputs_.m_gasMoleNumbers();
        if (l_pressureBar < 1.0 || l_pressureBar > 1000.0)
        {
            cerr << "ERROR: Pressure out of range. Please check the inputs!" << endl;
            exit(-1);
        }

        m_outPut_.open(m_outPutFileName_.c_str());
        m_outPut_ << "Temperature(K) " << "Pressure(bar) ";

        for (int i = 0; i < l_masterPrint->size(); i++)
        {
            m_outPut_ << l_masterPrint->at(i) << " ";
        }
        for (int i = 0; i < l_secondaryPrint->size(); i++)
        {
            m_outPut_ << l_secondaryPrint->at(i) << " ";
        }
        for (int i = 0; i < l_mineralPrint->size(); i++)
        {
            m_outPut_ << l_mineralPrint->at(i) << " ";
        }
        for (int i = 0;i < l_waterMolaltiyPrint->size();i++)
        {
            m_outPut_ << l_waterMolaltiyPrint->at(i) << " ";
        }
        /*for (int i = 0;i < l_gasMoleFractionPrint->size();i++)
        {
            m_outPut_ << l_gasMoleFractionPrint->at(i) << " ";
        }    */    
        for (int i = 0;i < l_gasMoleNumberPrint->size();i++)
        {
            m_outPut_ << l_gasMoleNumberPrint->at(i) << " ";
        }

        m_outPut_ << endl;

        for (double l_temperatureK : *l_temperatureList)
        {
            m_calc_(l_temperatureK, l_pressureBar);
        }
        m_outPut_ << endl;
        m_outPut_.close();
        return;
    }
    else if (m_geochemInputs_.m_getMoleFractionList()->size()>0)
    {
        double l_pressureBar = *m_geochemInputs_.m_getPressure()->begin();
        double l_temperatureK = *m_geochemInputs_.m_getTemperature()->begin();
        vector<string>* l_masterPrint = m_geochemInputs_.m_getMasterPrint();
        vector<string>* l_secondaryPrint = m_geochemInputs_.m_getSecondaryPrint();
        vector<string>* l_mineralPrint = m_geochemInputs_.m_getMineralPrint();
        vector<string>* l_waterMolaltiyPrint = m_geochemInputs_.m_waterMolalityPrint();
//        vector<string>* l_gasMoleFractionPrint = m_geochemInputs_.m_gasMoleFractionPrint();
        vector<string>* l_gasMoleNumberPrint = m_geochemInputs_.m_gasMoleNumbers();

        m_outPutFileName_ += "_moleFractions.txt";

        if (l_pressureBar < 1.0 || l_pressureBar > 1000.0)
        {
            cerr << "ERROR: Pressure out of range. Please check the inputs!" << endl;
            exit(-1);
        }

        m_outPut_.open(m_outPutFileName_.c_str());
        m_outPut_ << "Temperature(K) " << "Pressure(bar) ";

        for (int i = 0; i < m_geochemInputs_.m_getGasComponentName()->size(); i++)
        {
            m_outPut_ << m_geochemInputs_.m_getGasComponentName()->at(i) << " ";
        }

        for (int i = 0; i < l_masterPrint->size(); i++)
        {
            m_outPut_ << l_masterPrint->at(i) << " ";
        }
        for (int i = 0; i < l_secondaryPrint->size(); i++)
        {
            m_outPut_ << l_secondaryPrint->at(i) << " ";
        }
        for (int i = 0; i < l_mineralPrint->size(); i++)
        {
            m_outPut_ << l_mineralPrint->at(i) << " ";
        }
        for (int i = 0; i < l_waterMolaltiyPrint->size(); i++)
        {
            m_outPut_ << l_waterMolaltiyPrint->at(i) << " ";
        }
        /*for (int i = 0;i < l_gasMoleFractionPrint->size();i++)
        {
            m_outPut_ << l_gasMoleFractionPrint->at(i) << " ";
        }    */
        for (int i = 0; i < l_gasMoleNumberPrint->size(); i++)
        {
            m_outPut_ << l_gasMoleNumberPrint->at(i) << " ";
        }

        m_outPut_ << endl;

        m_totalGasMoleNumber_ = m_geochemInputs_.m_getGasMoleNumber();
        if (m_totalGasMoleNumber_ < 0.0)
        {
            cout << "Warning: total gas mole number is 0.0. Check the inputs!" << endl;
            m_gasIn_ = false;
            return;
        }

        int l_NoCalc = 0;
        for (vector<double> l_gasMoleFrac : *m_geochemInputs_.m_getMoleFractionList())
        {
            map<string, double> l_gasMoleFraction;
            double l_totalGasMoleFraction = 0;

            for (int i = 0; i < l_gasMoleFrac.size(); i++)
            {
                l_totalGasMoleFraction += l_gasMoleFrac[i];
            }

            if (l_totalGasMoleFraction < 1.0e-9)
            {
                cerr << "ERROR: GCM: Invalid input for gas mole fractions! Please check the inputs!" << endl;
                exit(-1);
                return;
            }

            for (int i = 0; i < l_gasMoleFrac.size(); i++)
            {
                double l_molefraction_ = l_gasMoleFrac[i] / l_totalGasMoleFraction;
                double l_moleNumber_ = l_molefraction_ * m_totalGasMoleNumber_;
                l_gasMoleFraction.insert(map<string, double>::value_type(m_geochemInputs_.m_getGasComponentName()->at(i), l_molefraction_));
                m_inputGasMoleNumber_[m_geochemInputs_.m_getGasComponentName()->at(i)] = l_moleNumber_;
            }

            if (l_NoCalc++ == 0)
            {
                m_gasWater_.generateGasCompName(m_inputGasMoleNumber_);
            }

            m_calc_(l_temperatureK, l_pressureBar);
        }
        m_outPut_ << endl;
        m_outPut_.close();
        return;
    }
    else
    {
        cerr << "ERROR: out of considerations. Please check the inputs!" << endl;
        m_outPut_.close();
        exit(-1);
    }
}

void GeoChemCalc::calculate(string a_inputFileName)
{
    m_condition_.masterChange = m_gasWater_.finalWater.masterSpecies_molality;
    m_condition_.gasPhaseIn = false;
    m_condition_.mass_of_H2O = 1.0;
    m_massOfH2OInWater = 1.0;
    m_outPutFileName_ = a_inputFileName;
    if (m_geochemInputs_.m_isKinetics)
    {
        kineticReaction();
        return;
    }
    else if (m_geochemInputs_.m_getTempPressList()->size() > 0)
    {
        //pair<double,double> l_temperatureK_pressureBar = *m_geochemInputs_.m_getTempPressList()->begin();
        list<pair<double,double>>* l_temperaure_pressureList = m_geochemInputs_.m_getTempPressList();
        m_outPutFileName_ += "_temperature_pressure.out";

        vector<string>* l_masterPrint = m_geochemInputs_.m_getMasterPrint();
        vector<string>* l_secondaryPrint = m_geochemInputs_.m_getSecondaryPrint();
        vector<string>* l_mineralPrint = m_geochemInputs_.m_getMineralPrint();
        vector<string>* l_waterMolaltiyPrint = m_geochemInputs_.m_waterMolalityPrint();
        vector<string>* l_gasMoleFractionPrint = m_geochemInputs_.m_gasMoleFractionPrint();

        m_outPut_.open(m_outPutFileName_.c_str());
        m_outPut_<<"Temperature(K) " << "Pressure(bar) ";

        for (int i = 0;i < l_masterPrint->size();i++)
        {
            m_outPut_ << l_masterPrint->at(i) << " ";
        }
        for (int i = 0;i < l_secondaryPrint->size();i++)
        {
            m_outPut_ << l_secondaryPrint->at(i) << " ";
        }
        for (int i = 0;i < l_mineralPrint->size();i++)
        {
            m_outPut_ << l_mineralPrint->at(i) << " ";
        }
        for (int i = 0;i < l_waterMolaltiyPrint->size();i++)
        {
            m_outPut_ << l_waterMolaltiyPrint->at(i) << " ";
        }
        for (int i = 0;i < l_gasMoleFractionPrint->size();i++)
        {
            m_outPut_ << l_gasMoleFractionPrint->at(i) << " ";
        }

        m_outPut_ << endl;

        for (pair<double, double> l_temperatureK_pressureBar : *l_temperaure_pressureList)
        {
            double l_temperatureK = l_temperatureK_pressureBar.first;
            double l_pressureBar = l_temperatureK_pressureBar.second;
            if (l_pressureBar < 1.0 || l_pressureBar > 1500.0)
            {
                cerr << "ERROR: Pressure out of range. Please check the inputs!" << endl;
                exit(-1);
            }

            if (l_temperatureK < 273.15 || l_temperatureK >533.15)
            {
                cerr << "ERROR: Temperature out of range. Please check the inputs!" << endl;
                exit(-1);
            }

            LDBLE atm_1bar = 0.9869233;
            m_condition_.pressure = l_pressureBar * atm_1bar;
            m_condition_.temperature = l_temperatureK - 273.15;
            updateGasSolubility(l_temperatureK, l_pressureBar);
            m_gasWater_.updateReaction(m_condition_);

            map<string, double> eqMineral;
            map<string, minerals>::iterator it_min = m_gasWater_.finalMineral.begin();
            for (;it_min != m_gasWater_.finalMineral.end();it_min++)
            {
                string l_name = it_min->first;
                double l_mineralMoleNumber = it_min->second.final_mole;
                eqMineral[l_name] = l_mineralMoleNumber;
            }
            m_condition_.eqMineralChange = eqMineral;
            m_print_(l_temperatureK, l_pressureBar);
        }
        m_outPut_ << endl;
        m_outPut_.close();
        return;
    }
    else if (m_geochemInputs_.m_getPressure()->size() > 1 || ((m_geochemInputs_.m_getPressure()->size()==1)&&(m_geochemInputs_.m_getTemperature()->size()==1)))
    {
        double l_temperatureK =* m_geochemInputs_.m_getTemperature()->begin();
        list<double>* l_pressureList = m_geochemInputs_.m_getPressure();
        m_outPutFileName_ += "_pressure.out";
        
        vector<string>* l_masterPrint = m_geochemInputs_.m_getMasterPrint();
        vector<string>* l_secondaryPrint = m_geochemInputs_.m_getSecondaryPrint();
        vector<string>* l_mineralPrint = m_geochemInputs_.m_getMineralPrint();
        vector<string>* l_waterMolaltiyPrint = m_geochemInputs_.m_waterMolalityPrint();
        vector<string>* l_gasMoleFractionPrint = m_geochemInputs_.m_gasMoleFractionPrint();
        if (l_temperatureK < 273.15 || l_temperatureK >533.15)
        {
            cerr << "ERROR: Temperature out of range. Please check the inputs!" << endl;
            exit(-1);
        }

        m_outPut_.open(m_outPutFileName_.c_str());
        m_outPut_ << "Pressure(bar) ";

        for (int i = 0;i < l_masterPrint->size();i++)
        {
            m_outPut_ << l_masterPrint->at(i) << " ";
        }
        for (int i = 0;i < l_secondaryPrint->size();i++)
        {
            m_outPut_ << l_secondaryPrint->at(i) << " ";
        }
        for (int i = 0;i < l_mineralPrint->size();i++)
        {
            m_outPut_ << l_mineralPrint->at(i) << " ";
        }
        for (int i = 0;i < l_waterMolaltiyPrint->size();i++)
        {
            m_outPut_ << l_waterMolaltiyPrint->at(i) << " ";
        }
        for (int i = 0;i < l_gasMoleFractionPrint->size();i++)
        {
            m_outPut_ << l_gasMoleFractionPrint->at(i) << " ";
        }

        m_outPut_ << endl;

        for (double l_pressureBar : *l_pressureList)
        {
            if (l_pressureBar < 1.0 || l_pressureBar > 1500.0)
            {
                cerr << "ERROR: Pressure out of range. Please check the inputs!" << endl;
                exit(-1);
            }
        //    m_massOfH2OInWater = 1.0;
            LDBLE atm_1bar = 0.9869233;
            m_condition_.pressure = l_pressureBar*atm_1bar;
            m_condition_.temperature = l_temperatureK - 273.15;
            updateGasSolubility(l_temperatureK, l_pressureBar);
            m_gasWater_.updateReaction(m_condition_);

            map<string, double> eqMineral;
            map<string, minerals>::iterator it_min = m_gasWater_.finalMineral.begin();
            for (;it_min != m_gasWater_.finalMineral.end();it_min++)
            {
                string l_name = it_min->first;
                double l_mineralMoleNumber = it_min->second.final_mole;
                eqMineral[l_name] = l_mineralMoleNumber;
            }
            m_condition_.eqMineralChange = eqMineral;
            m_print_(l_pressureBar);
        }
        m_outPut_.close();
        return;
    }
    else if (m_geochemInputs_.m_getTemperature()->size() > 1)
    {
        double l_pressureBar = *m_geochemInputs_.m_getPressure()->begin();
        list<double>* l_temperatureList = m_geochemInputs_.m_getTemperature();        
        m_outPutFileName_ +=  "temperature.txt";

        vector<string>* l_masterPrint = m_geochemInputs_.m_getMasterPrint();
        vector<string>* l_secondaryPrint = m_geochemInputs_.m_getSecondaryPrint();
        vector<string>* l_mineralPrint = m_geochemInputs_.m_getMineralPrint();
        vector<string>* l_waterMolaltiyPrint = m_geochemInputs_.m_waterMolalityPrint();
        vector<string>* l_gasMoleFractionPrint = m_geochemInputs_.m_gasMoleFractionPrint();
        if (l_pressureBar < 1.0 || l_pressureBar > 1000.0)
        {
            cerr << "ERROR: Pressure out of range. Please check the inputs!" << endl;
            exit(-1);
        }

        m_outPut_.open(m_outPutFileName_.c_str());
        m_outPut_ << "Temperature(K) ";

        for (int i = 0; i < l_masterPrint->size(); i++)
        {
            m_outPut_ << l_masterPrint->at(i) << " ";
        }
        for (int i = 0; i < l_secondaryPrint->size(); i++)
        {
            m_outPut_ << l_secondaryPrint->at(i) << " ";
        }
        for (int i = 0; i < l_mineralPrint->size(); i++)
        {
            m_outPut_ << l_mineralPrint->at(i) << " ";
        }    
        for (int i = 0;i < l_waterMolaltiyPrint->size();i++)
        {
            m_outPut_ << l_waterMolaltiyPrint->at(i) << " ";
        }
        for (int i = 0;i < l_gasMoleFractionPrint->size();i++)
        {
            m_outPut_ << l_gasMoleFractionPrint->at(i) << " ";
        }
        m_outPut_ << endl;

        for (double l_temperatureK : *l_temperatureList)
        {
            if (l_temperatureK < 273.15 || l_temperatureK > 533.15)
            {
                cerr << "ERROR: Temperature out of range. Please check the inputs!" << endl;
                exit(-1);
            }
    //        m_massOfH2OInWater = 1.0;
            LDBLE atm_1bar = 0.9869233;
            m_condition_.temperature = l_temperatureK-273.15;
            m_condition_.pressure = l_pressureBar * atm_1bar;
            updateGasSolubility(l_temperatureK, l_pressureBar);
            m_gasWater_.updateReaction(m_condition_);

            map<string, double> eqMineral;
            map<string, minerals>::iterator it_min = m_gasWater_.finalMineral.begin();
            for (;it_min != m_gasWater_.finalMineral.end();it_min++)
            {
                string l_name = it_min->first;
                double l_mineralMoleNumber = it_min->second.final_mole;
                eqMineral[l_name] = l_mineralMoleNumber;
            }
            m_condition_.eqMineralChange = eqMineral;
            m_print_(l_temperatureK);
        }
        m_outPut_.close();
        return;
    }
    else
    {
        cerr << "ERROR: out of considerations. Please check the inputs!" << endl;
        m_outPut_.close();
        exit(-1);
    }
}

void GeoChemCalc::kineticReaction()
{
    list<double>* l_timeList = m_geochemInputs_.m_getSimTime(); // Time Unit: Day
    double l_temperatureK = *m_geochemInputs_.m_getTemperature()->begin();
    double l_pressureBar = *m_geochemInputs_.m_getPressure()->begin();
    LDBLE atm_1bar = 0.9869233;
    m_condition_.temperature = l_temperatureK - 273.15;
    m_condition_.pressure = l_pressureBar * atm_1bar;

    m_outPutFileName_ += "_kinetic_time.out";

    vector<string>* l_masterPrint = m_geochemInputs_.m_getMasterPrint();
    vector<string>* l_secondaryPrint = m_geochemInputs_.m_getSecondaryPrint();
    vector<string>* l_mineralPrint = m_geochemInputs_.m_getMineralPrint();

    m_outPut_.open(m_outPutFileName_.c_str());
    m_outPut_ << "Time(day) ";

    for (int i = 0; i < l_masterPrint->size(); i++)
    {
        m_outPut_ << l_masterPrint->at(i) << " ";
    }
    for (int i = 0; i < l_secondaryPrint->size(); i++)
    {
        m_outPut_ << l_secondaryPrint->at(i) << " ";
    }
    for (int i = 0; i < l_mineralPrint->size(); i++)
    {
        m_outPut_ << l_mineralPrint->at(i) << " ";
    }
    m_outPut_ << endl;

    double lastTime = 0.0;
   /* map<string, minerals>::iterator it_min = m_gasWater_.finalMineral.begin();
    for (; it_min != m_gasWater_.finalMineral.end(); it_min++)
    {
        m_condition_.mineralsUpdate[it_min->first] = it_min->second.final_mole;       
    }*/

    m_condition_.mineralsUpdate = m_inputEqMinerals_;


    if (m_geochemInputs_.m_getPrintInitial())
    {
        m_print_(-1);
    }
//    kineticReaction_withTime(0, l_temperatureK, l_pressureBar);
    for (double time : *l_timeList)
    {
        double deltaTime;
        
        deltaTime = time - lastTime;
        
        double reactionTimeInSecond = deltaTime * 24.0 * 3600.0;
        
        kineticReaction_withTime(reactionTimeInSecond, l_temperatureK, l_pressureBar);
        lastTime = time;

        m_print_(time);
    }
}

void GeoChemCalc::kineticReaction_withTime(double a_Time_inSecond, double a_temperatureK, double a_pressureBar)
{
  //  double l_localPressure = 212;
    updateGasSolubility(a_temperatureK, a_pressureBar);
    m_condition_.deltaTime = a_Time_inSecond;
    m_gasWater_.updateReaction(m_condition_);

    m_massOfH2OInWater = m_gasWater_.finalWater.mass_of_H2O / 1.0;
 //   m_massOfH2OInWater =  1.0;

    /*
    It looks there is no need to update mineral information, because the mineral change is already stored in class phreeqc.
    */

 //   map<string, double> eqMineral;
    map<string, minerals>::iterator it_min = m_gasWater_.finalMineral.begin();
    for (; it_min != m_gasWater_.finalMineral.end(); it_min++)
    {
        m_condition_.mineralsUpdate[it_min->first] = it_min->second.final_mole;
    }
//    m_condition_.eqMineralChange = eqMineral;

    map<string, double> kinMineral;
    map<string, kineticPhases>::iterator it_kin = m_gasWater_.finalKineticMineral.begin();
    for (; it_kin != m_gasWater_.finalKineticMineral.end();it_kin++)
    {
        string l_name = it_kin->first;
        double l_kinMineralMoleNumber = it_kin->second.final_mole;
        kinMineral[l_name] = l_kinMineralMoleNumber;
    }
    m_condition_.kineticMineralChange = kinMineral;
    
}

void GeoChemCalc::updateGasSolubilityWithIniWater(double a_temperatureK, double a_pressureBar, conditionChange &a_condition, map<string, double> &a_gasMoleNumber)
{
    if (a_temperatureK >= 373.15)
    {
        double l_Psaturation = waterSaturationPressure(a_temperatureK);
        if (a_pressureBar < l_Psaturation)
        {
            cerr << "ERROR: pressure is lower than water saturation pressure!" << endl;
            exit(-1);
        }
    }
    else
    {
        if (a_pressureBar < 1.0)
        {
            cerr << "ERROR: pressrure is too low!" << endl;
            exit(-1);
        }
    }
//    m_condition_.masterChange = m_inputWater_.masterSpecies_molality;

    map<string, double>::iterator it_gas = a_gasMoleNumber.begin();// m_gasMoleNumber_.begin();
    for (; it_gas != a_gasMoleNumber.end(); it_gas++)
    {
        string l_gasName = it_gas->first;

        if (l_gasName == "H2O")
        {
            continue;
        }

        string l_element = m_gasWater_.getElemName(l_gasName);
        double l_coeff = m_gasWater_.getElemCoeff(l_gasName);

        double l_dissolvedGasMoleNumber = a_condition.masterChange.find(l_element)->second*a_condition.mass_of_H2O / l_coeff;
        it_gas->second += l_dissolvedGasMoleNumber;
    }

    double totalMolenumbers = 0;

//    map<string, double> l_gasMoleNumber = m_inputGasMoleNumber_;

    double gasMoleFraction = m_gasWater_.calPhaseEquilibria_RealH2Omass(a_gasMoleNumber, m_massOfH2OInWater, a_condition.masterChange, a_temperatureK, a_pressureBar, totalMolenumbers);

    m_massOfH2OInWater = totalMolenumbers * (1.0 - gasMoleFraction) * m_gasWater_.xAqueous.back() * 18.105 / 1000.0;
    
    for (int i = 0; i < (int)a_gasMoleNumber.size();i++)
    {
        double l_gasComponentMole = totalMolenumbers * gasMoleFraction * m_gasWater_.yGas[i];
        a_gasMoleNumber.find(m_gasWater_.gasCompName[i])->second = l_gasComponentMole;

        if (m_gasWater_.gasCompName[i] == "H2O")
        {
            continue;
        }
        a_condition.masterChange.find(m_gasWater_.getElemName(m_gasWater_.gasCompName[i]))->second = totalMolenumbers * (1.0 - gasMoleFraction) * m_gasWater_.xAqueous[i] * m_gasWater_.getElemCoeff(m_gasWater_.gasCompName[i]) / m_massOfH2OInWater;
    }

    a_condition.mass_of_H2O = m_massOfH2OInWater;
}

void GeoChemCalc::updateGasSolubility(double a_temperatureK ,double a_pressureBar)
{
    double l_temperature = a_temperatureK;
    double l_pressure = a_pressureBar;
    if (a_temperatureK >= 373.15)
    {
        double l_Psaturation = waterSaturationPressure(a_temperatureK);
        if (a_pressureBar < l_Psaturation)
        {
            cerr << "ERROR: pressure is lower than water saturation pressure!" << endl;
            exit(-1);
        }
    }
    else
    {
        if (a_pressureBar < 1.0)
        {
            cerr << "ERROR: pressrure is too low!" << endl;
            exit(-1);
        }
    }
    m_condition_.masterChange = m_gasWater_.finalWater.masterSpecies_molality;

    map<string, double>::iterator it_gas = m_gasMoleNumber_.begin();
    for (; it_gas != m_gasMoleNumber_.end(); it_gas++)
    {
        string l_gasName = it_gas->first;

        if (l_gasName == "H2O")
        {
            continue;
        }

        string l_element = m_gasWater_.getElemName(l_gasName);
        double l_coeff = m_gasWater_.getElemCoeff(l_gasName);

        double l_dissolvedGasMoleNumber = m_gasWater_.finalWater.masterSpecies_moleNumber.find(l_element)->second / l_coeff;
        it_gas->second += l_dissolvedGasMoleNumber;
    }

    double totalMolenumbers = 0;
    double gasMoleFraction = m_gasWater_.calPhaseEquilibria_RealH2Omass(m_gasMoleNumber_, m_massOfH2OInWater, m_condition_.masterChange, l_temperature, l_pressure, totalMolenumbers);
    
    m_massOfH2OInWater = totalMolenumbers * (1.0 - gasMoleFraction) * m_gasWater_.xAqueous.back()*18.105/1000.0;
    cout << "mass of H2O in water: " << m_massOfH2OInWater <<",  Gas mole fraction: "<<gasMoleFraction << endl;
    for (int i = 0; i < (int)m_gasMoleNumber_.size();i++)
    {
        double l_gasComponentMole = totalMolenumbers * gasMoleFraction * m_gasWater_.yGas[i];
        m_gasMoleNumber_.find(m_gasWater_.gasCompName[i])->second = l_gasComponentMole;

        if(m_gasWater_.gasCompName[i]=="H2O") 
        {
            continue;
        }
        m_condition_.masterChange.find(m_gasWater_.getElemName(m_gasWater_.gasCompName[i]))->second = totalMolenumbers * (1.0 - gasMoleFraction) * m_gasWater_.xAqueous[i] * m_gasWater_.getElemCoeff(m_gasWater_.gasCompName[i]) / m_massOfH2OInWater;
    }

    m_condition_.mass_of_H2O = m_massOfH2OInWater;
//    m_condition_.mass_of_H2O = 1.0;
}

void GeoChemCalc::m_print_(double a_firstColume)
{
    vector<string>* l_masterPrint = m_geochemInputs_.m_getMasterPrint();
    vector<string>* l_secondaryPrint = m_geochemInputs_.m_getSecondaryPrint();
    vector<string>* l_mineralPrint = m_geochemInputs_.m_getMineralPrint();
    vector<string>* l_printWaterMolality = m_geochemInputs_.m_waterMolalityPrint();
    vector<string>* l_printGasMolefraction = m_geochemInputs_.m_gasMoleFractionPrint();
    if (a_firstColume > -0.000001)
    {
        m_outPut_ << a_firstColume << " ";
    }
    else
    {
        m_outPut_ << "Initial ";
    }

    for (int i = 0; i < l_masterPrint->size(); i++)
    {
        string l_name = l_masterPrint->at(i);
        if (l_name == "PH" || l_name == "pH" || l_name == "Ph" || l_name == "ph")
        {
            m_outPut_ << m_gasWater_.finalWater.ph << " ";
            continue;
        }
        map<string, double> ::iterator it_master = m_gasWater_.finalWater.masterSpecies_molality.find(l_name);
        if (it_master != m_gasWater_.finalWater.masterSpecies_molality.end())
        {
            m_outPut_ << it_master->second << " ";
            continue;
        }
        else
        {
            cerr << "ERROR: cannot find the master species " << l_name << ". Please check the inputs and database file!" << endl;
            exit(-1);
        }

    }
    for (int i = 0; i < l_secondaryPrint->size(); i++)
    {
        string l_name = l_secondaryPrint->at(i);

        map<string, double> ::iterator it_secondary = m_gasWater_.finalWater.species_molality.find(l_name);
        if (it_secondary != m_gasWater_.finalWater.species_molality.end())
        {
            m_outPut_ << it_secondary->second << " ";
            continue;
        }
        else
        {
            m_outPut_ << "0  ";
            /*cerr << "ERROR: cannot file the secondary species " << l_name << ". Please check the inputs and database file!" << endl;
            exit(-1);*/
        }
    }
    for (int i = 0; i < l_mineralPrint->size(); i++)
    {
        string l_name = l_mineralPrint->at(i);
        map<string, minerals>::iterator it = m_gasWater_.finalMineral.find(l_name);
        if (it != m_gasWater_.finalMineral.end())
        {
            m_outPut_ << it->second.final_mole << " ";
            continue;
        }
        map<string, kineticPhases>::iterator it_kin = m_gasWater_.finalKineticMineral.find(l_name);
        if (it_kin != m_gasWater_.finalKineticMineral.end())
        {
            m_outPut_ << it_kin->second.final_mole << " ";
            continue;
        }
        else
        {
            cerr << "ERROR: cannot find the mineral " << l_name << ". Please check the inputs and database file!" << endl;
            exit(-1);
        }
    }

    for (int i = 0;i < l_printWaterMolality->size();i++)
    {
        int l_componentNumber = (int) m_gasWater_.xAqueous.size();
        double l_waterMoleFraction = m_gasWater_.xAqueous[l_componentNumber - 1];
        string l_name = l_printWaterMolality->at(i);
        int l_index = m_gasWater_.gasSequence.find(l_name)->second;
        double l_molality = m_gasWater_.xAqueous[l_index] / l_waterMoleFraction * 1000 / 18.105;

        m_outPut_ << l_molality << " ";
    }
    
    for (int i = 0;i < l_printGasMolefraction->size();i++)
    {
        string l_name = l_printGasMolefraction->at(i);
        int l_index = m_gasWater_.gasSequence.find(l_name)->second;
        double l_moleFraction = m_gasWater_.yGas[l_index];

        m_outPut_ << l_moleFraction << " ";
    }

    m_outPut_ << endl;
}

void GeoChemCalc::m_print_(double a_firstColume, double a_secondColume)
{
    vector<string>* l_masterPrint = m_geochemInputs_.m_getMasterPrint();
    vector<string>* l_secondaryPrint = m_geochemInputs_.m_getSecondaryPrint();
    vector<string>* l_mineralPrint = m_geochemInputs_.m_getMineralPrint();
    vector<string>* l_printWaterMolality = m_geochemInputs_.m_waterMolalityPrint();
    vector<string>* l_printGasMolefraction = m_geochemInputs_.m_gasMoleFractionPrint();
    vector<string>* l_printGasMoleNumber = m_geochemInputs_.m_gasMoleNumbers();

    m_outPut_ << a_firstColume << " "<<a_secondColume<<" ";

    if (m_geochemInputs_.m_getMoleFractionList()->size() > 0)
    {
        for (int i = 0; i < m_geochemInputs_.m_getGasComponentName()->size(); i++)
        {
            string l_gasCompName = m_geochemInputs_.m_getGasComponentName()->at(i);
            m_outPut_ << m_inputGasMoleNumber_.find(l_gasCompName)->second / m_totalGasMoleNumber_ << " ";
        }
    }

    for (int i = 0; i < l_masterPrint->size(); i++)
    {
        string l_name = l_masterPrint->at(i);
        if (l_name == "PH" || l_name == "pH" || l_name == "Ph" || l_name == "ph")
        {
            m_outPut_ << m_gasWater_.finalWater.ph << " ";
            continue;
        }
        map<string, double> ::iterator it_master = m_gasWater_.finalWater.masterSpecies_molality.find(l_name);
        if (it_master != m_gasWater_.finalWater.masterSpecies_molality.end())
        {
            m_outPut_ << it_master->second << " ";
            continue;
        }
        else
        {
            cerr << "ERROR: cannot file the master species " << l_name << ". Please check the inputs and database file!" << endl;
            exit(-1);
        }

    }
    for (int i = 0; i < l_secondaryPrint->size(); i++)
    {
        string l_name = l_secondaryPrint->at(i);

        map<string, double> ::iterator it_secondary = m_gasWater_.finalWater.species_molality.find(l_name);
        if (it_secondary != m_gasWater_.finalWater.species_molality.end())
        {
            m_outPut_ << it_secondary->second << " ";
            continue;
        }
        else
        {
            m_outPut_ << "0  ";
            /*cerr << "ERROR: cannot file the secondary species " << l_name << ". Please check the inputs and database file!" << endl;
            exit(-1);*/
        }
    }
    for (int i = 0; i < l_mineralPrint->size(); i++)
    {
        string l_name = l_mineralPrint->at(i);
        map<string, minerals>::iterator it = m_gasWater_.finalMineral.find(l_name);
        if (it != m_gasWater_.finalMineral.end())
        {
            m_outPut_ << it->second.final_mole << " ";
            continue;
        }
        map<string, kineticPhases>::iterator it_kin = m_gasWater_.finalKineticMineral.find(l_name);
        if (it_kin != m_gasWater_.finalKineticMineral.end())
        {
            m_outPut_ << it_kin->second.final_mole << " ";
            continue;
        }
        else
        {
            /*cerr << "ERROR: cannot find the mineral " << l_name << ". Please check the inputs and database file!" << endl;
            exit(-1);*/
            m_outPut_ << "0.0  ";
            continue;
        }
    }

    for (int i = 0;i < l_printWaterMolality->size();i++)
    {
        int l_componentNumber = (int)m_gasWater_.xAqueous.size();
        if (l_componentNumber <= 0)
        {
            m_outPut_ << "0  ";
            continue;
        }
        double l_waterMoleFraction = m_gasWater_.xAqueous[l_componentNumber - 1];
        string l_name = l_printWaterMolality->at(i);
        int l_index = m_gasWater_.gasSequence.find(l_name)->second;
        double l_molality = m_gasWater_.xAqueous[l_index] / l_waterMoleFraction * 1000 / 18.105;

        m_outPut_ << l_molality << " ";
    }

    /*for (int i = 0;i < l_printGasMolefraction->size();i++)
    {
        string l_name = l_printGasMolefraction->at(i);
        if (m_gasWater_.yGas.size() <= 0)
        {
            m_outPut_ << "0  ";
            continue;
        }
        int l_index = m_gasWater_.gasSequence.find(l_name)->second;
        double l_moleFraction = m_gasWater_.yGas[l_index];

        m_outPut_ << l_moleFraction << " ";
    }*/

    for (int i = 0;i < l_printGasMoleNumber->size();i++)
    {
        string l_name = l_printGasMoleNumber->at(i);
        
        map<string, double> ::iterator it = m_gasMoleNumber_.find(l_name);
        if (it == m_gasMoleNumber_.end())
        {
            m_outPut_ << "0  ";
            continue;
        }
        double l_moleNumber = it->second;

        m_outPut_ << l_moleNumber << " ";
    }

    m_outPut_ << endl;
}

double GeoChemCalc::waterSaturationPressure(double a_temperature)
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

    double ln_Ps_Pc = Tc / a_temperature * (a1 * tau + a2 * pow(tau, 1.5) + a3 * tau * tau * tau + a4 * pow(tau, 3.5) + a5 * pow(tau, 4) + a6 * pow(tau, 7.5));
    double Ps = Pc * exp(ln_Ps_Pc);
    return Ps;
}

void GeoChemCalc::m_calc_(double a_temperatureK, double a_pressureBar)
{
    if (a_pressureBar < 1.0 || a_pressureBar > 1500.0)
    {
        cerr << "ERROR: Pressure out of range. Please check the inputs!" << endl;
        exit(-1);
    }

    if (a_temperatureK < 273.15 || a_temperatureK >533.15)
    {
        cerr << "ERROR: Temperature out of range. Please check the inputs!" << endl;
        exit(-1);
    }

    LDBLE atm_1bar = 0.9869233;
    conditionChange l_condition;
    l_condition.pressure = a_pressureBar * atm_1bar;
    l_condition.temperature = a_temperatureK - 273.15;
    l_condition.masterChange = m_inputWater_.masterSpecies_molality;
    l_condition.gasPhaseIn = false;
    l_condition.mass_of_H2O = 1.0;
    m_massOfH2OInWater = 1.0;

    l_condition.kineticMineralChange = m_inputKineticMinerals_;
    l_condition.eqMineralChange = m_inputEqMinerals_;
    l_condition.mineralsUpdate = m_inputEqMinerals_;
    //    updateGasSolubility(l_temperatureK, l_pressureBar);
    map<string, double> l_gasMoleNumber = m_inputGasMoleNumber_;

    if (m_totalGasMoleNumber_ <= 1.0e-12)
    {
        m_gasWater_.updateReaction(l_condition);
    }
    else
    {
        int l_numberOfIterations = 0;
        while (true)
        {
            std::cout << "Iterations - " << l_numberOfIterations << endl;
            updateGasSolubilityWithIniWater(a_temperatureK, a_pressureBar, l_condition, l_gasMoleNumber);
            m_gasWater_.updateReaction(l_condition);
            l_condition.mass_of_H2O = m_massOfH2OInWater;
            double l_maxChange = 0;
            map<string, double> ::iterator it_master = m_gasWater_.finalWater.masterSpecies_molality.begin();
            for (;it_master != m_gasWater_.finalWater.masterSpecies_molality.end();it_master++)
            {
                string l_masterName = it_master->first;
                map<string, double> ::iterator it_conditonMaster = l_condition.masterChange.find(l_masterName);
                if (it_conditonMaster == l_condition.masterChange.end())
                {
                    cerr << "ERROR: cannot find the master species: " << l_masterName << endl;
                    exit(-1);
                }
                if (it_master->second > 0)
                {
                    bool l_smallValue = it_conditonMaster->second < 1.0e-9 && it_master->second < 1.0e-9;
                    if (!l_smallValue && l_maxChange < fabs(it_conditonMaster->second - it_master->second) / it_master->second)
                    {
                        l_maxChange = fabs(it_conditonMaster->second - it_master->second) / it_master->second;
                        //                    cout << l_masterName << endl;
                    }
                }
                else
                {
                    if (it_conditonMaster->second > 0 && l_maxChange < fabs(it_conditonMaster->second - it_master->second) / 1e-5)
                    {
                        l_maxChange = fabs(it_conditonMaster->second - it_master->second) / 1e-5;
                        //                    cout << l_masterName << endl;
                    }
                }
                it_conditonMaster->second = it_master->second;
            }

            if (l_maxChange < 1.0e-4)
            {
                break;
            }

            map<string, minerals>::iterator it_eqMineral = m_gasWater_.finalMineral.begin();
            for (;it_eqMineral != m_gasWater_.finalMineral.end();it_eqMineral++)
            {
                string l_eqMineralName = it_eqMineral->first;
                map<string, double>::iterator it_conditionMineral = l_condition.eqMineralChange.find(l_eqMineralName);
                if (it_conditionMineral == l_condition.eqMineralChange.end())
                {
                    /*    cerr << "ERROR: cannot find the mineral: " << l_eqMineralName << endl;
                        exit(-1);*/
                    l_condition.eqMineralChange[l_eqMineralName] = it_eqMineral->second.final_mole;
                    continue;
                }
                it_conditionMineral->second = it_eqMineral->second.final_mole;
            }

            map<string, kineticPhases> ::iterator it_kinMineral = m_gasWater_.finalKineticMineral.begin();
            for (;it_kinMineral != m_gasWater_.finalKineticMineral.end();it_kinMineral++)
            {
                string l_kinMineralName = it_kinMineral->first;
                map<string, double> ::iterator it_conditionKinMineral = l_condition.kineticMineralChange.find(l_kinMineralName);
                if (it_conditionKinMineral == l_condition.kineticMineralChange.end())
                {
                    l_condition.kineticMineralChange[l_kinMineralName] = it_kinMineral->second.final_mole;
                    continue;
                    /*cerr << "ERROR: cannot find the mineral: " << l_kinMineralName << endl;
                    exit(-1);*/
                }
                it_conditionKinMineral->second = it_kinMineral->second.final_mole;
            }

            //                m_gasWater_.finalWater.masterSpecies_molality
            l_numberOfIterations++;
            if (!m_redox_)
            {
                break;
            }
        }
        cout << "Number of iterations: " << l_numberOfIterations << endl;
        cout << "mass of H2O in water: " << m_massOfH2OInWater << endl;
    }

    m_gasMoleNumber_ = l_gasMoleNumber;
    cout << "Mass of solution: " << m_gasWater_.finalWater.density * m_gasWater_.finalWater.volume << endl;
    cout << "Mass of H2O: " << m_gasWater_.finalWater.mass_of_H2O << endl;
    cout << "Geochem calculation finished!" << endl;
    m_print_(a_temperatureK, a_pressureBar);
//    m_outPut_ << endl;
}