#pragma once
#include "../MultiPhaseEquilibria/GasIniWaterEquilibria.h"
#include "GeoChemInputs.h"
class GeoChemCalc
{
public:
    GeoChemCalc();
    ~GeoChemCalc();

//    void readIn(string a_inputFileName);
    void initilize(map<string, vector<string> >* a_gcmData);
    void calculate(string a_inputFileName);
    void calculateWithIniWater(string a_inputFileName);

private:
    void kineticReaction();
    void kineticReaction_withTime(double a_time, double a_temperatureK, double a_pressureBar);
    void updateGasSolubility(double a_temperatureK, double a_pressureBar);
    void updateGasSolubilityWithIniWater(double a_temperatureK, double a_pressureBar, conditionChange &a_conditon, map<string, double > &a_gasMoleNumber);
    void m_calc_(double a_temperatureK, double a_pressureBar);
    map<string, double> m_gasMoleNumber_;
    double m_massOfH2OInWater; // Kg

    GasIniWaterEquilibria m_gasWater_;
    GeoChemInputs m_geochemInputs_;

    conditionChange m_condition_;
    bool m_gasIn_;
    string m_outPutFileName_;

    ofstream m_outPut_;

    void m_print_(double a_firstColume);
    void m_print_(double a_firstColume, double a_secondColume);
    double waterSaturationPressure(double a_temperatureK);

    waterSolution_results m_inputWater_;
    map<string, double>   m_inputEqMinerals_;
    map<string, double>   m_inputKineticMinerals_;
    map<string, double>   m_inputGasMoleNumber_;
    double m_totalGasMoleNumber_;
    bool m_redox_;
    
};

