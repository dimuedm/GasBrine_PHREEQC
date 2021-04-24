/*
Gas component database for cubic EOS
Created by Jun Li, Dimue company, Wuhan.
Created: 2018-4-2
Modifited: 2018-4-2
*/

#ifndef GASCOMPONENTDATABASE_H_
#define GASCOMPONENTDATABASE_H_
#include <map>
#include "Component.h"
using namespace std;

class GasComponentDataBase
{
public:
    GasComponentDataBase();
    ~GasComponentDataBase();

public:
    void initialize();
    map<string, Component>      *getData();
    vector< vector<double> >    *getKij();
    int                          getSize();
    int                          isInitialized = 0;

private:
    map<string, Component>       m_gasDatabase_;
    vector< vector<double> >     m_kij_;
    int                          m_sizeOfDatabase_;
};




#endif