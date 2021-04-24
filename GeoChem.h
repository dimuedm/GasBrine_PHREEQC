/*
   Interface class for geochemistry cacluation.
   Created by Jun Li, 2017-12-27
   Dimue Company, Wuhan, China
*/

#ifndef GEOCHEM_H_
#define GEOCHEM_H_

#include "Framework\Base.h"
#include "Phreeqc_ReSoC.h"
//#include "Phreeqc.h"

using namespace std;

class Phreeqc;

class GeoChem
{
public:
    GeoChem();
    ~GeoChem();

    GeoChem(const GeoChem &);
    GeoChem& operator =(GeoChem &obj1);
//    void initialize(istream *a_inputStream, istream *a_dataStream);
//    void initialize(string *a_inputstring, string *a_databaseFileName);
    void initialize(string *a_inputstring, istream *a_dataStream);


    void initialize_database(istream *a_databaseStream);
    void initialize(string *a_inputString);
    void m_assign(Phreeqc *a_phreeqc);

    void save_solution();
    void print_solution(int n_user);

    void update_timeStep(conditionChange &a_newCondition);
    bool geoChemExist;
    Phreeqc* m_getPhreeqc();
//    int m_getMasterNumber();
    std::map<int, std::string> * getMasterSpecies();
    std::map<std::string, double>  *getMasterSpeciesMolality();
    waterSolution_results * getWater();
    std::map<std::string, minerals> *getMineral();

private:
    Phreeqc *m_phreeqc_;
    gasPhase_results m_gas_;
    waterSolution_results m_water_;
    std::map<std::string, minerals> m_mineral_;
//    int m_numberOfMasterSpecies_;

};






#endif