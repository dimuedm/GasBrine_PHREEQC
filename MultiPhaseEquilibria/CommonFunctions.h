//Deal with PVT input; get input data from InputManager.
//Created by Jun Li, Dimue company, Wuhan.
//Created: 2017-4-25
//Modified: 2017-7-14, function nonlinearSolve_Newton is deleted!
//Modified: 2017-7-17, function calGasMole - if RReqn_Current<1e-4 and iteration number>10000, we assume the function will return true!

#ifndef COMMONFUNCTIONS_H
#define COMMONFUNCTIONS_H
#include <vector>


using namespace std;


class CommonFunctions
{
public:
    CommonFunctions();
    ~CommonFunctions();

    //    void nonlinearSolve_Newton(vector<vector<double> > &Jacob, vector<double> &function);
    void linearSoverLU(vector<vector<double> > &A, vector<double> &b, vector<double> &x);
    vector<vector<double> > lu;
    int nMatrixSize;
    vector<int> index;

    bool calcGasMole(double &beta_GasMole, vector<double> &zFeed, vector<double> &Kvalue);
private:
    void LUdcmp(vector<vector<double> > &A);
    double calRReqn(double &beta_gasMole, vector<double> &zFeed, vector<double> & Kvalue);
    double calRReqnPrim(double &beta, vector<double>&zFeed, vector<double>&Kvalue);

};





#endif