#pragma once
#include "../Framework/Base.h"
#include "../UtilIO.h"

class GeoChemInputs
{
public:
    GeoChemInputs();
    ~GeoChemInputs();

    bool readInputs(map<string, vector<string> >* a_dataBlock);
    bool m_isKinetics;
    int m_printType;
    string* m_getDataBaseFile();
    string* m_getPhrqcInputFile();
    list<double>* m_getSimTime(); //Unit: Day;
    list<double>* m_getPressure();
    list<double>* m_getTemperature();
    list<pair<double, double>>* m_getTempPressList();
    list < vector<double> >* m_getMoleFractionList();

    vector<string>* m_getGasComponentName();

    vector<string>* m_getMasterPrint();
    vector<string>* m_getSecondaryPrint();
    vector<string>* m_getMineralPrint();
    vector<string>* m_waterMolalityPrint();
    vector<string>* m_gasMoleFractionPrint();
    vector<string>* m_gasMoleNumbers();

    map<string, double>* m_getGasMoleFraction();
    double m_getGasMoleNumber();
    bool   m_getRedox();
    bool   m_getPrintInitial();
private:

    bool m_readInputs_(vector<string>* a_gcmData);
    bool m_readPhrqcInputFile_(string& a_phrqc_input);

    string m_databaseFile_;
    string m_phrqc_inputFile_;
    map<int, string> m_phrqcInputBlock_;
    list<double> m_simTime_;
    list<double> m_pressureList_;
    list<double> m_temperatureList_;
    list<pair<double, double>> m_temperature_pressureList_;
    list<vector<double>> m_gasMoleFractionList_;
    vector<string> m_gasComponentName_;
    map<string, double> m_gasMoleFraction_;

    vector<string> m_printMaster_;
    vector<string> m_printSecondary_;
    vector<string> m_printMineral_;
    vector<string> m_printWaterMolality_;
    vector<string> m_printGasMolefraction_;
    vector<string> m_printGasMoleNumbers_;
    
    double m_gasMoleNumber_;

    bool m_timeListIn_;

    bool m_redox_;

    bool m_printInitialEquil_; // The equilibrium between water and minerals. Usually, gas is not included.
};
