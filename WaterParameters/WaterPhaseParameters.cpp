
#include "WaterPhaseParameters.h"


WaterPhaseParameters::WaterPhaseParameters()
{
   
}

WaterPhaseParameters::~WaterPhaseParameters()
{
}

//void WaterPhaseParameters::initialize(PVTenvironment *environment_, double temperature, double pressure)
//{
//	(void)pressure;
//
//    environment = environment_;
//    iComponentNumber = environment->getComponentNumber();
//
//    if (temperature != environment->m_temperature)
//    {
//        environment->updateCubicParam(temperature);
//    }
//
//    components = environment->getComponents();
//
//    activityCoeff.resize(iComponentNumber);
//    equilibriumConst.resize(iComponentNumber);
//}

//void WaterPhaseParameters::getEquiliConst(double temperature, double pressure)
//{
//    if (environment->m_getWaterPropModel())
//    {
//        return;
//    }
//
//    int NO_forWaterSolu;
//    for (int i = 0; i < iComponentNumber; i++)
//    {
//        NO_forWaterSolu = (*components)[i].NO_forWaterSolu;
//        equilibriumConst[i] = calEquiliConst(NO_forWaterSolu, temperature, pressure);
//    }
//}

//void WaterPhaseParameters::getActivityCoefficient(double molalityNaCl, double temperature, double pressure)
//{
//    if (environment->m_getWaterPropModel())
//    {
//        return;
//    }
//
//    double lamdaTemp, zetaTemp;
//    int NO_forWaterSolu;
//    for (int i = 0; i < iComponentNumber; i++)
//    {
//        NO_forWaterSolu = (*components)[i].NO_forWaterSolu;
//        lamdaTemp = calActivityCoefficient_lamda(NO_forWaterSolu, temperature, pressure);
//        zetaTemp = calActivityCoefficient_zeta(NO_forWaterSolu, temperature, pressure);
//        activityCoeff[i] = 2.0*molalityNaCl*lamdaTemp + molalityNaCl*molalityNaCl*zetaTemp;
//        activityCoeff[i] = exp(activityCoeff[i]);
//    }
//}

double WaterPhaseParameters::specificVolume(int NO_forWaterSolu, double temperature, double pressure, double molalityNaCl)
{
    switch (NO_forWaterSolu)
    {
    case 2://CO2
        return calCO2SpecificVolume(temperature, molalityNaCl);
        break;
    case 3://CH4
        return calCH4SpecificVolume(temperature, pressure, molalityNaCl);
        break;
    case 4://N2
        return calN2SpecificVolume(temperature, pressure, molalityNaCl);
        break;
    case 5://H2S
        return calH2SSpecificVolume(temperature, pressure, molalityNaCl);
        break;
    default:
        return 1.0;
        break;
    }
}

// From 
double WaterPhaseParameters::calCO2SpecificVolume(double temperature, double molalityNaCl)
{
    double a1 = -0.636988853022e2;
    double a2 = 0.704083840800e0;
    double a3 = -0.175530451238e-2;
    double a4 = 0.162606713882e-5;
    double a5 = -0.902263067104e2;
    double a6 = 0.668916598035e0;
    double a7 = -0.160072919102e-2;
    double a8 = 0.115890935719e-5;
    double a9 = -0.136682312625e2;
    double a10 = 0.106931009961e0;
    double a11 = -0.293911238637e-3;
    double a12 = 0.306116540808e-6;
    double a13 = 0.668319293470e1;
    double a14 = -0.527779110962e-1;
    double a15 = 0.139806267407e-3;
    double a16 = -0.126734373947e-6;

    double t2 = temperature*temperature;
    double t3 = t2*temperature;
    double vphi = a1 + a2*temperature + a3*t2 + a4*t3 + (a5 + a6*temperature + a7*t2 + a8*t3)*molalityNaCl + (a9 + a10*temperature + a11*t2 + a12*t3)*molalityNaCl*molalityNaCl + (a13 + a14*temperature + a15*t2 + a16*t3)*molalityNaCl*molalityNaCl*molalityNaCl;
    return vphi;
}

double WaterPhaseParameters::calCH4SpecificVolume(double temperature, double pressure, double molalityNaCl)
{// Unit:cm3/mole
    double a6 = -0.0132253629347427e0;
    double a7 = 3.13045055468136e-5;
    double a8 = 2.26047549807049e0;
    double a9 = -2.90633722703171e-8;
    double a10 = 3.44079701140235e-9;
    double a11 = -6.11808900183315e-10;
    double dmu_dp = a6 + a7*temperature + a8 / temperature + a9*temperature*temperature + a10*2.0*pressure*temperature + 3.0*a11*pressure*pressure;

    double dLamda_dp = 8.13889426952425e-5 + 2.0*2.30245272398387e-11*pressure * temperature - 2.34735334659631e-7*temperature;

    double vphi = dmu_dp + dLamda_dp*molalityNaCl*2.e0;
    vphi *= R*temperature*10.e0; //1.0*e-5 is to transfer unit to cm3/mole;

    return vphi;
}

double WaterPhaseParameters::calH2SSpecificVolume(double temperature, double pressure, double molalityNaCl)
{
    double c6 = 0.823684902091498e0;
    double c7 = -0.00222376062074741e0;
    double c8 = -101.113160193415e0;
    double c9 = 1.96421386690219e-6;
    double c10 = 3.77937696024575e-8;
    double c11 = 2.92076340659316e-9;
    double dmu_dp = c6 + c7*temperature + c8 / temperature + c9*temperature*temperature + 2.0*c10*pressure* temperature + 3.0*c11*pressure*pressure;

    double b4 = -0.00770031630821377e0;
    double b5 = -4.36067311771693e-8;
    double b6 = 2.50009500063678e-5;
    double dLamda_dp = b4 + 2.0*b5*pressure * temperature + b6*temperature;

    double vphi = dmu_dp + dLamda_dp*molalityNaCl*2.e0;
    vphi *= R*temperature*10.e0; //1.0*e-5 is to transfer unit to cm3/mole;

    return vphi;
}

double WaterPhaseParameters::calN2SpecificVolume(double temperature, double pressure, double molalityNaCl)
{
    double c6 = -0.49542866e-3;
    double c7 = 0.12698747e-5;
    double c8 = 0.51411144e0;
    double c9 = -0.64733978e-4;
    double dmu_dp = c6 + c7*temperature + c8 / temperature + 2.0e0*c9*pressure / temperature;
   
    double dLamda_dp = -0.13711527e-4 + 2.0e0*0.71037217e-5*pressure / temperature;

    double vphi = dmu_dp + dLamda_dp*molalityNaCl*2.e0;
    vphi *= R*temperature*10.e0; //1.0*e-5 is to transfer unit to cm3/mole;

    return vphi;
}

double WaterPhaseParameters::calActivityCoefficient_lamda(int NO_forWaterSolu, double temperature, double pressure)
{
    double lamdaTemple;
    double c1 = -2.4434074e0;
    double c2 = 0.36351795e-2;
    double c3 = 0.44747364e3;
    double c6 = -0.13711527e-4;
    double c9 = 0.71037217e-5;

    double a1 = -0.718015778827693e0;
    double a2 = 0.00101441901569244e0;
    double a3 = 166.242391217727e0;
    double a4 = 8.13889426952425e-5;
    double a5 = 2.30245272398387e-11;
    double a6 = -2.34735334659631e-7;

    double b1 = 1.33924585177518;
    double b2 = -0.00203344940110559;
    double b3 = -203.201235035567;
    double b4 = -0.00770031630821377;
    double b5 = -4.36067311771693e-8;
    double b6 = 2.50009500063678e-5;

    switch (NO_forWaterSolu)
    {
    case 1://H2O
        return 1.0;
        break;
    case 2://CO2
        lamdaTemple = -0.31312239e0 + 0.5532647e-3 * temperature + 0.75844401e2 / temperature - 0.18950519e-3 * pressure + 0.71628762e-6 * pressure*temperature - 0.1458572e-9 * pressure*pressure * temperature;
        return lamdaTemple;
        break;
    case 3://CH4
        //     lamdaTemple = -5.7066455e-1 + 7.2997588e-4 * temperature + 1.5176903e2 / temperature + 3.1927112e-5 * pressure - 1.642651e-5 * pressure / temperature;
        lamdaTemple = a1 + a2*temperature + a3 / temperature + a4*pressure + a5*pressure*pressure * temperature + a6*pressure*temperature;
        return lamdaTemple;
        break;
    case 4://N2
        lamdaTemple = c1 + c2*temperature + c3 / temperature + c6*pressure + c9*pressure*pressure / temperature;
        return lamdaTemple;
        break;
    case 5://H2S
        //    lamdaTemple = 1.03658689e0 - 1.1784797e-3 * temperature - 1.7754826e2 / temperature - 4.5313285e-4 * pressure + 47.75165e0*pressure / temperature / temperature;
        lamdaTemple = b1 + b2*temperature + b3 / temperature + b4*pressure + b5*pressure*pressure * temperature + b6*pressure*temperature;
        return lamdaTemple;
        break;
    default:
        return 1.0;
        break;
    }
}

double WaterPhaseParameters::calActivityCoefficient_zeta(int NO_forWaterSolu, double temperature, double pressure)
{
    double zetaTemple = 0.0;
    switch (NO_forWaterSolu)
    {
    case 1://H2O
        return 1.0;
        break;
    case 2://CO2
        zetaTemple = 0.34096802e-2 - 0.27671084e-4 * temperature - 0.83847525e-7 * pressure*temperature + 0.34225403e-10 * pressure*pressure * temperature;
        return zetaTemple;
        break;
    case 3://CH4
        //        zetaTemple = -2.9990084e-3;
        zetaTemple = -0.00165439125284811e0;
        return zetaTemple;
        break;
    case 4://N2
        zetaTemple = -0.58071053e-2;
        return zetaTemple;
    case 5://H2S
        //      zetaTemple = -0.010274152e0;
        zetaTemple = 0.00254763081450856e0;
        return zetaTemple;
        break;
    default:
        return 1.0;
        break;
    }
}

double WaterPhaseParameters::calEquiliConst(int NO_forWaterSolu, double temperature, double pressure)
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
    default:
        return 1.e10;
        break;
    }
}



double WaterPhaseParameters::calWaterEquiliConst(double temperature, double pressure)
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

double  WaterPhaseParameters::calCH4EquiliConst(double temperature, double pressure)
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

double  WaterPhaseParameters::calCO2EquiliConst(double temperature, double pressure)
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

double WaterPhaseParameters::calH2SEquiliConst(double temperature, double pressure)
{
    double c1 = -825.351676920731;
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
    equiliConstValue = exp(equiliConstValue)*55.508;
    return equiliConstValue;
}

double  WaterPhaseParameters::calN2EquiliConst(double temperature, double pressure)
{
    double c1 = -23.093813e0;
    double c2 = 0.56048525e-1;
    double c3 = 0.98808898e4;
    double c4 = -0.51091621e-4;
    double c5 = -0.13220298e7;
    double c6 = -0.49542866e-3;
    double c7 = 0.12698747e-5;
    double c8 = 0.51411144e0;
    double c9 = -0.64733978e-4;
    double equiliConstValue = c1 + c2*temperature + c3 / temperature + c4*temperature*temperature + c5 / temperature / temperature + c6*pressure + c7*pressure*temperature + c8*pressure / temperature + c9*pressure*pressure / temperature;
    equiliConstValue = exp(equiliConstValue)*55.508e0;
    return equiliConstValue;
}

double  WaterPhaseParameters::calO2EquiliConst(double temperature, double pressure)
{
	(void)temperature; (void)pressure;

    return 0;
}