#include <assert.h>
#include <ostream>

#include "EOSPhaseProp.h"
#include "EnumDef.h"
#include "Component.h"
#include "GasComponentDataBase.h"

#include "../Framework/Util/UtilLog.h"
#include "../Framework/Util/UtilString.h"

using namespace std;


EOSPhaseProp::EOSPhaseProp()
{

}

EOSPhaseProp::~EOSPhaseProp()
{

}

void EOSPhaseProp::initialize(vector<string> componentFormula, vector<double> zFeed_, GasComponentDataBase *a_componentDatabase_, PhaseType a_phase_, double temperature_, double pressure_, EOSChoice a_eosType_)
{
    temperature = temperature_;
    pressure = pressure_;
    eosType = a_eosType_;
    Phase = a_phase_;

    if (componentFormula.size() != zFeed_.size())
    {
        LOG(ERROR) << "EOSPhaseProp: component formula size and zFeed size are not equal!";
        exit(-1);
    }

    ComponentNumber = componentFormula.size();

    componentProperty = a_componentDatabase_->getData();

    if (componentProperty->size()<=0)
    {
        LOG(ERROR) << "EOSPhaseProp: Empty component database input!";
        exit(-1);
    }

    components.resize(ComponentNumber);
    for (int i = 0; i < ComponentNumber; i++)
    {
        UtilString stringProcess;
        //string tempFormula = stringProcess.toUpper(componentFormula[i]);
        string tempFormula = componentFormula[i];
        map<string, Component>::iterator it = componentProperty->find(tempFormula);
        if (it != componentProperty->end())
        {
            components[i] = it->second;
            components[i].calcEOSProperty(temperature, eosType);
        }
        else
        {
            LOG(ERROR) << "EOSPhaseProp: component is not existing in database!";
            exit(-1);
        }
    }

    da_dx.resize(ComponentNumber);
    db_dx.resize(ComponentNumber);
    dc_dx.resize(ComponentNumber);
    dA_dx.resize(ComponentNumber);
    dB_dx.resize(ComponentNumber);

 
   // kij=environment->getBinaryCoeff();
    kij.resize(ComponentNumber);
    for (int i = 0; i < ComponentNumber; i++)
    {
        kij[i].resize(ComponentNumber);
        for (int j = 0; j < ComponentNumber; j++)
        {
            if (i == j)
            {
                kij[i][j] = 0;
            }
            else
            { 
                int i_database = components[i].number;
                int j_database = components[j].number;
                kij[i][j] = (*a_componentDatabase_->getKij())[i_database][j_database];
            }
        }
    }

    zFeed.resize(ComponentNumber);
    for (int i = 0; i < ComponentNumber; i++)
    {
        zFeed[i] = zFeed_[i];
    }

    fugacity.resize(ComponentNumber);
    df_dp.resize(ComponentNumber);
    logPhi.resize(ComponentNumber);

    //initialize the df_dx, dlogPhi_dx;
    df_dx.resize(ComponentNumber);
    dlogPhi_dx.resize(ComponentNumber);
    for (int i = 0; i < ComponentNumber; i++)
    {
        df_dx[i].resize(ComponentNumber);
        dlogPhi_dx[i].resize(ComponentNumber);
    }

    EOSChoice eos = a_eosType_;

    if (eos == PR76 || eos == PR78)
    {
        u = 2.0;
        w = -1.0;
    }
    else if (eos == RK||eos==SRK)
    {
        u = 1.0;
        w = 0;
    }
    else
    {
        LOG(ERROR) << "Unknown input EOS choice" << endl;
        exit(-1);
    }

    calMoleWeight();

    cal_abc();
    calZ();
}

void EOSPhaseProp::calMoleWeight()
{
    int i;
    moleWeight = 0.0;
    for (i = 0; i < ComponentNumber; i++)
    {
        moleWeight += zFeed[i] * (components)[i].MoleWeight;
    }
}

void EOSPhaseProp::calDensity(bool derivative)
{

    moleVolume = Z*R*temperature / pressure - c;
    densityMole = 1.0 / moleVolume*1.0e5;  // 10e-5*m3/mole for moleVolume; mole/m3 for densityMole;
    moleVolume *= 1.e-2; //  Unit:L/mole

    densitySI = densityMole * moleWeight / 1000.0; //Unit: KG/M3

    if (derivative)
    {
        double rtp = -R*temperature*1.e5 / (pressure*moleVolume*moleVolume);

        dDensity_dP = rtp*(-Z / pressure + dZ_dp)   * moleWeight / 1000.e0;
        dDensity_dT = -1.0 / moleVolume / moleVolume*(Z*R / pressure + R*temperature / pressure*dZ_dT)    * moleWeight / 1000.e0  *1.e5;

        dDensity_dx.resize(ComponentNumber);
        for (int i = 0; i < ComponentNumber; i++)
        {
            //Modified: Jun Li 05/09/2017
            //dDensity_dx[i] = rtp*(dZ_dx[i] - dc_dx[i])    * moleWeight / 1000.e0;   //This is how Cao's code does!
            dDensity_dx[i] = (R*temperature / pressure*dZ_dx[i] - dc_dx[i])*(-1.0 / moleVolume / moleVolume) *1.e5 *moleWeight / 1000.0;
        }
    }


}

void EOSPhaseProp::calViscosity(double density) //the input density should be mole density;
{
    int i;
    double viscAt1atm = 0.0, temple = 0.0;

    calComponentViscosity();

    for (i = 0; i < ComponentNumber; i++)
    {
        viscAt1atm += zFeed[i] * sqrt((components)[i].MoleWeight)*viscosityComponent[i];
        temple += zFeed[i] * sqrt((components)[i].MoleWeight);
    }
    viscAt1atm /= temple;

    calPseudoCriticalProperties();
    double densityRelative = density / densityPseudoCritical;

    temple = 0.1023 + 0.023364*densityRelative + 0.058533*densityRelative*densityRelative - 0.040758*pow(densityRelative, 3) + 0.0093324*pow(densityRelative, 4);
    temple = pow(temple, 4) - 1.0e-4;
    temple /= pow(Tc_pseudo , 1.0 / 6.0) / sqrt(moleWeight) / pow(Pc_pseudo , 2.0 / 3.0);

    viscosity = viscAt1atm + temple;

    //Currntly 0;
    dViscosity_dP = 0.0;
    dViscosity_dT = 0.0;
    dViscosity_dx.resize(ComponentNumber);
    for (i = 0; i < ComponentNumber; i++)
    {
        dViscosity_dx[i] = 0.0;
    }
}

void EOSPhaseProp::calComponentViscosity()
{
    int i;
    double viscParam;
    double Tr;
    viscosityComponent.resize(ComponentNumber);

    for (i = 0; i < ComponentNumber; i++)
    {
        viscParam = pow((components)[i].Tc, 1.0 / 6.0) / sqrt((components)[i].MoleWeight) / pow((components)[i].Pc, 2.0 / 3.0);
        Tr = temperature / (components)[i].Tc;

        if (Tr <= 1.5)
        {
            viscosityComponent[i] = 3.4e-4*pow(Tr, 0.94) / viscParam;
        }
        else
        {
            viscosityComponent[i] = 1.778e-4*pow(4.58*Tr - 1.67, 0.625) / viscParam;
        }
    }
}

inline void EOSPhaseProp::calPseudoCriticalProperties()
{
    Pc_pseudo = 0.0;
    Tc_pseudo = 0.0;
    Vc_pseudo = 0.0;

    for (int i = 0; i < ComponentNumber; i++)
    {
        Pc_pseudo += zFeed[i] * (components)[i].Pc;
        Tc_pseudo += zFeed[i] * (components)[i].Tc;
        Vc_pseudo += zFeed[i] * ((components)[i].Zcrit*R*(components)[i].Tc) / (components)[i].Pc;
    }
    densityPseudoCritical = 1.0 / Vc_pseudo*1.0e5; //RA 17-08-2017
}

// Cubic model is implemented!
void EOSPhaseProp::calcFugacity(bool calcDireviatives)
{
    double B_inv = 1.0 / B;
    double b_inv = 1.0 / b;
    double u2_4w = sqrt(1.0*(u*u-4*w));
    double ZB1_inv = 1.0 / (2.0*Z + B*(u - u2_4w));
    double ZB2_inv = 1.0 / (2.0*Z + B*(u + u2_4w));
    double LOGZB2_ZB1 = log(ZB1_inv / ZB2_inv);
    double ZB_log = log(Z - B);
    double A_B = A*B_inv;
    double A_BB = A_B*B_inv;
    double ABuw = A_B / u2_4w;
    double a_inv = 1.0 / a;
    double bi_b;
    int i;

    
// Calculate f
    for (i = 0; i < ComponentNumber; i++)
    {
        bi_b = (components)[i].bCubic*b_inv;
        logPhi[i] = bi_b*(Z - 1) - ZB_log + ABuw*(bi_b - da_dx[i] * a_inv)*LOGZB2_ZB1;
        fugacity[i] = pressure*zFeed[i] * exp(logPhi[i]);
    }

// Calculate df_dp
    double ZminusB_inv = 1.0 / (Z - B);
    for (i = 0; i < ComponentNumber; i++)
    {
        bi_b = (components)[i].bCubic*b_inv;
        df_dp[i] = bi_b*dZ_dp - (dZ_dp - dB_dp)*ZminusB_inv;
        df_dp[i] += (bi_b - da_dx[i] * a_inv) / u2_4w*((dA_dp / B - A_BB*dB_dp)*LOGZB2_ZB1 + A_B*((2 * dZ_dp + dB_dp*(u + u2_4w))*ZB2_inv - (2 * dZ_dp + dB_dp*(u - u2_4w))*ZB1_inv));
        df_dp[i] *= fugacity[i];
        df_dp[i] += fugacity[i] / pressure;
    }

// Calculate df_dx
    int j;
    for (i = 0; i < ComponentNumber; i++)
    {
        bi_b = (components)[i].bCubic*b_inv;
        for (j = 0; j < ComponentNumber; j++)
        {
            dlogPhi_dx[i][j] = bi_b*(dZ_dx[j] - (Z - 1.0)*b_inv*db_dx[j]) - (dZ_dx[j] - dB_dx[j])*ZminusB_inv;
            dlogPhi_dx[i][j] += (bi_b - da_dx[i] * a_inv) / u2_4w*((dA_dx[j] * B_inv - A_BB*dB_dx[j])*LOGZB2_ZB1 + A_B*((2 * dZ_dx[j] + dB_dx[j] * (u + u2_4w))*ZB2_inv - (2 * dZ_dx[j] + dB_dx[j] * (u - u2_4w))*ZB1_inv));
            dlogPhi_dx[i][j] += A_B / u2_4w*LOGZB2_ZB1*(-bi_b*b_inv*db_dx[j] + da_dx[i] * da_dx[j] * (a_inv*a_inv) - 2 * (components)[i].aCubic_sqrt*(components)[j].aCubic_sqrt*(1.0 - (kij)[i][j])*a_inv);
            df_dx[i][j] = dlogPhi_dx[i][j] * fugacity[i]; 
            if (i == j) df_dx[i][j] += fugacity[i] / zFeed[i];
        }
    }
}

void EOSPhaseProp::calcFugacity()
{
    double B_inv = 1.0 / B;
    double b_inv = 1.0 / b;
    double u2_4w = sqrt(1.0*(u*u - 4 * w));
    double ZB1_inv = 1.0 / (2.0*Z + B*(u - u2_4w));
    double ZB2_inv = 1.0 / (2.0*Z + B*(u + u2_4w));
    double LOGZB2_ZB1 = log(ZB1_inv / ZB2_inv);
    double ZB_log = log(Z - B);
    double A_B = A*B_inv;
    double A_BB = A_B*B_inv;
    double ABuw = A_B / u2_4w;
    double a_inv = 1.0 / a;
    double bi_b;
    int i;

    // Calculate f
    for (i = 0; i < ComponentNumber; i++)
    {
        bi_b = (components)[i].bCubic*b_inv;
        logPhi[i] = bi_b*(Z - 1) - ZB_log + ABuw*(bi_b - da_dx[i] * a_inv)*LOGZB2_ZB1;
        fugacity[i] = pressure*zFeed[i] * exp(logPhi[i]);
    }
}

void EOSPhaseProp::cal_abc()
{
    a = b = c = 0;
    int i, j;
    double daij_dT;
    for (i = 0; i < ComponentNumber; i++)
    {
        da_dx[i] = 0.0;
        double xa_sqrt = zFeed[i] * (components)[i].aCubic_sqrt;
        for (j = 0; j < ComponentNumber; j++)
        {
            double temp = zFeed[j] * (1.0 - (kij)[i][j])*(components)[j].aCubic_sqrt;
            da_dx[i] += temp;
            a += temp*xa_sqrt;
        }
        b += zFeed[i] * (components)[i].bCubic;
        c += zFeed[i] * (components)[i].cCubic;

        //calculate a, b, c derivatives;
        da_dx[i] *= 2 * (components)[i].aCubic_sqrt;
        db_dx[i] = (components)[i].bCubic;
        dc_dx[i] = (components)[i].cCubic;

        dA_dx[i] = da_dx[i] * pressure / (R*R*temperature*temperature);
        dB_dx[i] = db_dx[i] * pressure / R / temperature;
    }

    A = a*pressure/(R*R*temperature*temperature);
    B = b*pressure / R / temperature;
    
    dA_dp = A / pressure;
    dB_dp = B / pressure;
    wB = w*B;
    wBB = wB*B;

    // Calulate da_dT
    da_dT = 0.0;
    for (i = 0; i < ComponentNumber; i++)
    {
        for (j = 0; j < ComponentNumber; j++)
        {
            daij_dT = (1.0 - (kij)[i][j]) / 2.0*(((components)[i].aCubic) * ((components)[j].daCubic_dT) + ((components)[j].aCubic)*((components)[i].daCubic_dT)) / sqrt(((components)[j].aCubic) * ((components)[i].aCubic));
            da_dT += zFeed[i] * zFeed[j] * daij_dT;
        }
    }

    dA_dT = A / a*da_dT - 2.0*A / temperature;
    dB_dT = -B / temperature;

    //d2a_dxidT
    d2a_dxidT.resize(ComponentNumber);
    vector<double> d2a_dxidT1(ComponentNumber), d2a_dxidT2(ComponentNumber);
    for (i = 0; i < ComponentNumber; i++)
    {
        d2a_dxidT[i] = 0.0;

        for (j = 0; j < ComponentNumber; j++)
        {
            d2a_dxidT1[i] += zFeed[j] * (1.0 - (kij)[i][j])*(components)[j].aCubic_sqrt;
        }
        d2a_dxidT1[i] *= (components)[i].daCubic_dT / (components)[i].aCubic_sqrt;

        for (j = 0; j < ComponentNumber; j++)
        {
            d2a_dxidT2[i] += zFeed[j] * (1.0 - (kij)[i][j])*(components)[j].daCubic_dT / (components)[j].aCubic_sqrt;
        }
        d2a_dxidT2[i] *= (components)[i].aCubic_sqrt;

        d2a_dxidT[i] = d2a_dxidT1[i] + d2a_dxidT2[i];
    }
}


void EOSPhaseProp:: calDT()
{
    int i;
    
    //dlogPhi_dT
    dlogPhi_dT.resize(ComponentNumber);
    double ABuw = A / B / sqrt(u*u - 4*w);
    double ZBuw1 = 2 * Z + B*(u + sqrt(u*u - 4 * w));
    double ZBuw2 = 2 * Z + B*(u - sqrt(u*u - 4 * w));
    double logZBuw1_ZBuw2 = log( ZBuw1 /  ZBuw2 );
    for (i = 0; i < ComponentNumber; i++)
    {
        dlogPhi_dT[i] = (components)[i].bCubic / b*dZ_dT;
        dlogPhi_dT[i] += -1.0 / (Z - B)*(dZ_dT-dB_dT);
        dlogPhi_dT[i] += dA_dT / A*ABuw*((components)[i].bCubic / b - da_dx[i] / a)*logZBuw1_ZBuw2;
        dlogPhi_dT[i] += (-1.0) / B*dB_dT*ABuw*((components)[i].bCubic / b - da_dx[i] / a)*logZBuw1_ZBuw2;
        dlogPhi_dT[i] += ABuw*(-1.0)*(d2a_dxidT[i] / a - da_dT / a / a*da_dx[i])*logZBuw1_ZBuw2;
        dlogPhi_dT[i] += ABuw*((components)[i].bCubic / b - da_dx[i] / a) *ZBuw2 / ZBuw1*((2 * dZ_dT + dB_dT*(u + sqrt(u*u - 4 * w))) / ZBuw2 + (-1.0)*ZBuw1 / ZBuw2 / ZBuw2*(2 * dZ_dT + dB_dT*(u - sqrt(u*u - 4 * w))));
    }

    //df_dT;
    df_dT.resize(ComponentNumber);
    for (i = 0; i < ComponentNumber; i++)
    {
        df_dT[i] = fugacity[i] * dlogPhi_dT[i];
    }
}



void EOSPhaseProp::calZ()
{
    double coeff2, coeff1, coeff0;

    coeff2 = -(1 + B - u*B);
    coeff1 = A + w*B*B - u*B*B - u*B;
    coeff0 = -(A*B + w*B*B + w*B*B*B);
    double a1 = -coeff2*coeff2 / 3.0 + coeff1;
    double b1 = coeff0 + 2.0*coeff2*coeff2*coeff2 / 27.0 - coeff2*coeff1 / 3.0;
    double ja = b1*b1 / 4.0 + a1*a1*a1 / 27.0;

    double coeff2_3 = coeff2 / 3.0;

    if (ja > 0)
    {
        double ja_sqrt = sqrt(ja);
        double A1, B1;

        A1 = pow(fabs(-b1 / 2.0 + ja_sqrt), 0.33333333333333333333333333333333333333);
        B1 = pow(fabs(-b1 / 2.0 - ja_sqrt), 0.33333333333333333333333333333333333333);

        if (-b1 / 2.0 + ja_sqrt < 0) A1 = -A1;
        if (-b1 / 2.0 - ja_sqrt < 0) B1 = -B1;
        Z = A1 + B1 - coeff2_3;
    }
    else if (ja <= 0)
    {
        double z1, z2, z3;
        double fai = acos(-0.5*b1/sqrt(-a1*a1*a1/27.0));
        double a1_3_sqrt = 2 * sqrt(-a1*0.333333333333333333333333);
        double fai_3 = fai*0.333333333333333333333333;
        z1 = a1_3_sqrt*cos(fai_3) - coeff2_3;
        z2 = a1_3_sqrt*cos(fai_3 + 0.6666666666666666666666667*PI) - coeff2_3;
        z3 = a1_3_sqrt*cos(fai_3 + 1.3333333333333333333333333*PI) - coeff2_3;

        if (Phase == OIL)
        {
             Z = FindMin(z1, z2, z3);					// min(Z) for oil phase
        }
        else if (Phase == GAS)
        {
             Z = FindMax(z1, z2, z3);				// max(Z) for gas phase
        }
        else
        { 
            assert(false); 
        }
    }

    //Calculate dZ_dp
// NOTE: PR76 EOS is used!
//    double u = 2, w = -1;
    double Z_2 = Z*Z;
    double ZZpZq_inv = 1.0 / (3.0*Z_2+2*coeff2*Z+coeff1);
    dZ_dp = -(Z_2*(u - 1)*dB_dp + Z*(dA_dp - u*dB_dp + 2 * (w - u)*B*dB_dp) - (A*dB_dp + B*dA_dp + 2 * wB*dB_dp + 3 * wBB*dB_dp))*ZZpZq_inv;


    //Caclulate dZ_dT
    double p = coeff2;
    double q = coeff1;
    double r = coeff0;
    double dp_dT = (u-1)*dB_dT;
    double dq_dT = dA_dT + 2 * w*B*dB_dT - 2 * u*B*dB_dT - u*dB_dT;
    double dr_dT = -(B*dA_dT + A*dB_dT + 2 * w*B*dB_dT + 3 * w*B*B*dB_dT);
    dZ_dT = -(Z*Z*dp_dT+Z*dq_dT+dr_dT) / (3*Z*Z+2*p*Z+q);

// Calculate dZ_dx;
    int i;
    dZ_dx.resize(ComponentNumber);
    for (i = 0; i < ComponentNumber; i++)
    {
        dZ_dx[i] = Z_2*(u - 1)*dB_dx[i] + Z*(dA_dx[i] - u*dB_dx[i] + 2 * (w - u)*B*dB_dx[i]) - (A*dB_dx[i] + B*dA_dx[i] + 2 * wB*dB_dx[i] + 3 * wBB*dB_dx[i]);
        dZ_dx[i] *= -ZZpZq_inv;
    }

}

inline double EOSPhaseProp::FindMax(double x1, double x2, double x3)
{
    double max = x1;
    if (x2>max) max = x2;
    if (x3>max) max = x3;

    return max;
}

inline double EOSPhaseProp::FindMin(double x1, double x2, double x3)
{
    const double tem = 1.0e20;
    if (x1<0.0) x1 = tem;
    if (x2<0.0) x2 = tem;
    if (x3<0.0) x3 = tem;


    double min = x1;
    if (x2<min) min = x2;
    if (x3<min) min = x3;
    return min;
}
