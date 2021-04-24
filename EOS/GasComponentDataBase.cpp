
#include "GasComponentDataBase.h"
#include "EnumDef.h"
using namespace std;

GasComponentDataBase::GasComponentDataBase()
{
}

GasComponentDataBase::~GasComponentDataBase()
{
}

void GasComponentDataBase::initialize()
{
    isInitialized = 1;
    m_sizeOfDatabase_ = 6;

    Component tempComponent;

    //H2O
    tempComponent.Name = "H2O";
    tempComponent.Tc = 647.3;
    tempComponent.Pc = 221.2;
    tempComponent.Omega = 0.3434;
    tempComponent.MoleWeight = 18.015;
    tempComponent.Zcrit = 0.307;
    tempComponent.s_cb = 0;
    tempComponent.number = 0;
    
    m_gasDatabase_["H2O"] = tempComponent;

    //CO2
    tempComponent.Name = "CO2";
    tempComponent.Tc = 304.20999;
    tempComponent.Pc = 73.83048;
    tempComponent.Omega = 0.2236;
    tempComponent.MoleWeight = 44.01;
    tempComponent.Zcrit = 0.2742;
    tempComponent.s_cb = -0.0817;
    tempComponent.number = 1;

    m_gasDatabase_["CO2"] = tempComponent;

    //CH4
    tempComponent.Name = "CH4";
    tempComponent.Tc = 190.56;
    tempComponent.Pc = 45.99042;
    tempComponent.Omega = 0.0115;
    tempComponent.MoleWeight = 16.043;
    tempComponent.Zcrit = 0.2884;
    tempComponent.s_cb = -0.1595;
    tempComponent.number = 2;

    m_gasDatabase_["CH4"] = tempComponent;

    //N2
    tempComponent.Name = "N2";
    tempComponent.Tc = 126.2;
    tempComponent.Pc = 33.99961;
    tempComponent.Omega = 0.0377;
    tempComponent.MoleWeight = 28.014;
    tempComponent.Zcrit = 0.2916;
    tempComponent.s_cb = -0.1927;
    tempComponent.number = 3;

    m_gasDatabase_["N2"] = tempComponent;

    //O2
    tempComponent.Name = "O2";
    tempComponent.Tc = 154.6;
    tempComponent.Pc = 49.8;
    tempComponent.Omega = 0.021;
    tempComponent.MoleWeight = 32.0;
    tempComponent.Zcrit = 0.307;
    tempComponent.s_cb = 0;
    tempComponent.number = 4;

    m_gasDatabase_["O2"] = tempComponent;

    //H2S
    tempComponent.Name = "H2S";
    tempComponent.Tc = 373.53;
    tempComponent.Pc = 89.62908;
    tempComponent.Omega = 0.1081; // 0.0942;
    tempComponent.MoleWeight = 34.082;
    tempComponent.Zcrit = 0.2831;
    tempComponent.s_cb = -0.1288;
    tempComponent.number = 5;

    m_gasDatabase_["H2S"] = tempComponent;

    //H2Sg
    tempComponent.Name = "H2Sg";
    tempComponent.Tc = 373.53;
    tempComponent.Pc = 89.62908;
    tempComponent.Omega = 0.1081;//0.0942;
    tempComponent.MoleWeight = 34.082;
    tempComponent.Zcrit = 0.2831;
    tempComponent.s_cb = -0.1288;
    tempComponent.number = 5;

    m_gasDatabase_["H2Sg"] = tempComponent;

    m_kij_.resize(m_sizeOfDatabase_);
    for (int i = 0; i < m_sizeOfDatabase_; i++)
    {
        (m_kij_[i]).resize(m_sizeOfDatabase_);
        m_kij_[i][i] = 0.0;
    }

    m_kij_[H2Og][CO2g] = 0.1896; m_kij_[H2Og][CH4g] = 0.485;  m_kij_[H2Og][N2g] = 0.32547/*0.4778*/; m_kij_[H2Og][O2g] = 0.20863/*0.48*/;  m_kij_[H2Og][H2Sg] = 0.105;
    m_kij_[CO2g][H2Og] = 0.1896; m_kij_[CO2g][CH4g] = 0.105;  m_kij_[CO2g][N2g] = -0.007/*0.0*/;    m_kij_[CO2g][O2g] = 0.1140/*0.0*/;   m_kij_[CO2g][H2Sg] = 0.0974;
    m_kij_[CH4g][H2Og] = 0.485;  m_kij_[CH4g][CO2g] = 0.105;  m_kij_[CH4g][N2g] = 0.025;  m_kij_[CH4g][O2g] = 0.025; m_kij_[CH4g][H2Sg] = 0.07;
    m_kij_[N2g][H2Og] = 0.32547/*0.4778*/;  m_kij_[N2g][CO2g] = -0.007/*0.0*/;     m_kij_[N2g][CH4g] = 0.025;  m_kij_[N2g][O2g] = -0.0119/*0.0*/;    m_kij_[N2g][H2Sg] = 0.13;
    m_kij_[O2g][H2Og] = 0.20863/*0.48*/;    m_kij_[O2g][CO2g] = 0.1140/*0.0*/;     m_kij_[O2g][CH4g] = 0.025;  m_kij_[O2g][N2g] = -0.0119/*0.0*/;    m_kij_[O2g][H2Sg] = 0.13;
    m_kij_[H2Sg][H2Og] = 0.105;  m_kij_[H2Sg][CO2g] = 0.0974; m_kij_[H2Sg][CH4g] = 0.07;  m_kij_[H2Sg][N2g] = 0.13;  m_kij_[H2Sg][O2g] = 0.13;
}

map<string, Component> * GasComponentDataBase::getData()
{
    return &m_gasDatabase_;
}

vector< vector<double> > * GasComponentDataBase::getKij()
{
    return &m_kij_;
}

int GasComponentDataBase::getSize()
{
    return m_sizeOfDatabase_;
}