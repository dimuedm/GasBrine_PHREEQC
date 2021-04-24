#ifndef UTILITIES_H
#define UTILITIES_H
#include <iostream>
#include <vector>
using namespace std;

const double g_constant = 9.81;
const double Pa_to_bar = 1.0/100000.0;
const double cent_to_unit = 1.0/1000.0;
const double unit_to_kilo = 1.0 / 1000.0;
const double day_to_sec = 24 * 60 * 60;
const double sec_to_day = 1.0/(24*60*60);
const double mD_to_m2 = 9.869233e-16;
const double mD_Pas_to_m2_bars = 9.869233e-11; //mD_to_m2/Pa_to_bar
const double mD_Pas_to_m2_barDay = 9.869233e-11/ sec_to_day; //mD_to_m2/Pa_to_bar

double minimumOf(double, double);

double maximumOf(double, double);

double sumArray(double [], int);

double maximumIn(vector<double > &);

class dataInterp
{
public:
    dataInterp();
    ~dataInterp();
    void m_updateData(vector<double > &, vector<double > &);
    double m_getValue(double &);
    double m_getValue(double &, double &);
private:

    double m_linearInterp_(double &, double &, double &, double &, double &);
    double m_linearInterp_(double &, double &, double &, double &, double &, double &);
    vector<double > x_, y_;
};
#endif // !UTILITIES_H





