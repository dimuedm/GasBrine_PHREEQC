#include <iostream>
#include <fstream>
#include "CommonFunctions.h"

using namespace std;

CommonFunctions::CommonFunctions()
{

}

CommonFunctions::~CommonFunctions()
{
}

void CommonFunctions::linearSoverLU(vector<vector<double> > &A, vector<double> &b, vector<double> &x)
{
    nMatrixSize = (int)b.size();

    if (nMatrixSize != A.size() || nMatrixSize != x.size())
    {
        cerr << "Error defination of linear equation problem!!!" << endl;
        exit(-1);
    }
    index.resize(nMatrixSize);

    LUdcmp(A);

    int i, ii = 0, ip, j;
    double sum;

    for (i = 0; i < nMatrixSize; i++)
    {
        x[i] = b[i];
    }

    for (i = 0; i < nMatrixSize; i++)
    {
        ip = index[i];
        sum = x[ip];

        x[ip] = x[i];

        if (ii != 0)
        {
            for (j = ii - 1; j < i; j++)
            {
                sum -= lu[i][j] * x[j];
            }
        }
        else if (sum != 0.0)
        {
            ii = i + 1;
        }

        x[i] = sum;
    }

    for (i = nMatrixSize - 1; i >= 0; i--)
    {
        sum = x[i];
        for (j = i + 1; j < nMatrixSize; j++)
        {
            sum -= lu[i][j] * x[j];
        }

        x[i] = sum / lu[i][i];
    }

}

void CommonFunctions::LUdcmp(vector<vector<double> > &A)
{
    const double TINY = 1.0e-40;
    int i, j, k, imax = -1;
    double max, temp;
    vector<double> scalFacor(nMatrixSize);

    lu = A;

    bool rowInterchange = false;

    for (i = 0; i < nMatrixSize; i++)
    {
        max = 0.0;
        for (j = 0; j < nMatrixSize; j++)
        {
            temp = fabs(lu[i][j]);
            if (temp > max)
            {
                max = temp;
            }
        }
        if (max == 0.0)
        {
            cerr << "Singular matrix in LUcmp!!!" << endl;
            exit(-1);
        }

        scalFacor[i] = max;
    }

    for (k = 0; k < nMatrixSize; k++)
    {
        max = 0.0;
        for (i = k; i < nMatrixSize; i++)
        {
            temp = scalFacor[i] * fabs(lu[i][k]);
            if (temp > max)
            {
                max = temp;
                imax = i;
            }
        }

        if (k != imax)
        {
            for (j = 0; j < nMatrixSize; j++)
            {
                temp = lu[imax][j];
                lu[imax][j] = lu[k][j];
                lu[k][j] = temp;
            }

            rowInterchange = true;
            scalFacor[imax] = scalFacor[k];
        }
        index[k] = imax;

        if (lu[k][k] == 0.0)
        {
            lu[k][k] = TINY;
        }

        for (i = k + 1; i < nMatrixSize; i++)
        {
            lu[i][k] /= lu[k][k];
            temp = lu[i][k];
            for (j = k + 1; j < nMatrixSize; j++)
            {
                lu[i][j] -= temp*lu[k][j];
            }
        }
    }
}

bool CommonFunctions::calcGasMole(double &beta_GasMole, vector<double> &zFeed, vector<double> &Kvalue)
{
    double RReqn0, RReqn1;
    double beta0 = 0.0, beta1 = 1.0;
    int i;

    if (zFeed.size() != Kvalue.size())
    {
        cerr << "Error size of composition or K value!" << endl;
    }

    RReqn0 = calRReqn(beta0, zFeed, Kvalue);
    RReqn1 = calRReqn(beta1, zFeed, Kvalue);

    if (RReqn0 <= 0.0)
    {
        beta_GasMole = -1.0;
        return false;
    }

    if (RReqn1 >= 0.0)
    {
        beta_GasMole = 2.0;
        return false;
    }

    //Find tighter boundary..
    double betaMin = 0.0, betaMax = 1.0;
    double LowerBound, HigherBound;
    int ComponentNumber = (int)zFeed.size();

    for (i = 0; i < ComponentNumber; i++)
    {
        if (Kvalue[i] > 1.0)
        {
            LowerBound = (Kvalue[i] * zFeed[i] - 1.0) / (Kvalue[i] - 1.0);
            if (betaMin < LowerBound) betaMin = LowerBound;
        }

        if (Kvalue[i] < 1.0)
        {
            HigherBound = (1.0 - zFeed[i]) / (1.0 - Kvalue[i]);
            if (betaMax > HigherBound) betaMax = HigherBound;
        }
    }

    beta_GasMole = 0.5*(betaMax + betaMin);

    //Find beta_GasMole by Newton iteration!
    const double eps_beta = 1.0e-6;
    const int MaxIteration_beta = 10000;
    int iterNum = 0;
    double RReqn_Current, RReqnPrim_Current;
    double beta_GasMole_new;
    while (true)
    {
        if (iterNum > MaxIteration_beta)
        {
            //          if (fabs(RReqn_Current) < 1.e-4)
            //         {
            return true;
            //         }
            //cerr << "Too many iterations for gas mole fraction calculation!!!" << endl;
            //exit(-1);
        }

        RReqn_Current = calRReqn(beta_GasMole, zFeed, Kvalue);
        //        cout << RReqn_Current << endl;

        if (fabs(RReqn_Current) <= eps_beta)
        {
            return true;
        }

        if (RReqn_Current > 0)
        {
            betaMin = beta_GasMole;
        }
        else
        {
            betaMax = beta_GasMole;
        }

        iterNum += 1;

        RReqnPrim_Current = calRReqnPrim(beta_GasMole, zFeed, Kvalue);

        //      cout << "RReqn_Current = " << RReqn_Current << endl;

        if (RReqnPrim_Current == 0.0)
        {
            cerr << "RReqnPrim is 0!!!" << endl;
            exit(1);
        }
        beta_GasMole_new = beta_GasMole + (-RReqn_Current / RReqnPrim_Current);
        //    cout << "beta_Gasmole = " << beta_GasMole_new << endl;
        if (beta_GasMole_new<betaMax && beta_GasMole_new>betaMin)
        {
            beta_GasMole = beta_GasMole_new;
        }
        else
        {
            beta_GasMole = 0.5*(betaMin + betaMax);
        }
    }
}

double CommonFunctions::calRReqn(double &beta, vector<double>&zFeed, vector<double> &Kvalue)
{
    double RReqn = 0.0;
    if (Kvalue.size() != zFeed.size())
    {
        cerr << "Error size of composition or K value!" << endl;
    }

    int i;
    int ComponentNumber = (int)zFeed.size();
    for (i = 0; i < ComponentNumber; i++)
    {
        RReqn += zFeed[i] * (Kvalue[i] - 1.0) / (1.0 - beta + beta*Kvalue[i]);
    }

    return RReqn;
}

double CommonFunctions::calRReqnPrim(double &beta, vector<double>&zFeed, vector<double>&Kvalue)
{
    double RReqnPrim = 0.0;
    if (Kvalue.size() != zFeed.size())
    {
        cerr << "Error size of composition or K value!" << endl;
    }

    int i;
    int ComponentNumber = (int)zFeed.size();
    for (i = 0; i < ComponentNumber; i++)
    {
        RReqnPrim -= zFeed[i] * (Kvalue[i] - 1.0)*(Kvalue[i] - 1.0) / (1.0 - beta + beta*Kvalue[i]);
    }

    return RReqnPrim;
}