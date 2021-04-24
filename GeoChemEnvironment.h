/*
To setup the common infitial informatin for geochemistry.
Created by Jun Li, 2017-12-27
Dimue company, Wuhan, China.
*/


#ifndef GEOCHEMENVIRONMENT_H_
#define GEOCHEMENVIRONMENT_H_

#include <map>
//#include "MeshData\MeshData.h"
#include "Framework\Util\UtilLog.h"
#include "Framework\Util\UtilFile.h"
//using namespace std;
#include "Input\UtilIO.h"


struct zone
{
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
};

class GeoChemEnvironment
{
public:
    GeoChemEnvironment();
    ~GeoChemEnvironment();

    bool m_GeochemEvironmentSetup(map <string, vector<string> > * a_pinputDataBlock_);
    /*map<int, istream*> * m_getPhrqcIstream();*/
    istream ** m_getDatabaseStream();
    map<int, string>* m_getPrqcInputData();
    string * m_getDatabaseFileName();
    map<int, zone> * m_getZones();
    int m_getDefaultZoneIndex();
    vector<int>* m_getInjectedSolution();

private:
    bool m_readInputs_(map <string, vector<string> > * a_pinputDataBlock_);
    bool m_readGCMDataBlock_(vector<string> * a_v_GCMData_);
    vector<int> m_zoneIndex_;
    int m_zoneNumber_;
    int m_defaultZoneIndex_;
    vector<int> m_injectedSolution_;
    map <int, zone> m_mp_Zones_;
    string m_phrqc_inputFile_;
    bool m_readPhrqcInputFile_(string &a_phrqc_input);
    map<int, string> m_phrqcInputBlock_;
//    map<int, istream*> m_phrqcIstram_;
    string m_databaseFile_;
//    string m_databaseFileData_;
    istream *m_databaseStream_;
    
};




#endif