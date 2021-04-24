#include "UtilitiesF.h"
#include "../Util/UtilLog.h"
double minimumOf(double S1, double S2)
{
	if (S1 < S2)
		return S1;
	else
		return S2;
}

double maximumIn(vector<double > &vec)
{
    double maxV = vec[0];
    int vecSize = vec.size();
    for (int vi = 1; vi < vecSize; vi++)
    {
        maxV = maximumOf(vec[vi],maxV);
    }
    return maxV;
}

double maximumOf(double S1, double S2)
{
	if (S1 > S2)
		return S1;
	else
		return S2;
}

double sumArray(double Array[], int sizeArray)
{
	double sumA = 0.0;
	for (int i = 0; i < sizeArray; i++)
	{
		sumA += Array[i];
	}
	return sumA;
}

dataInterp::dataInterp() {}
dataInterp::~dataInterp() {}

void dataInterp::m_updateData(vector<double > &a_x, vector<double > &a_y)
{
    x_ = a_x; y_ = a_y;
}

double dataInterp::m_getValue(double &xi)
{
    double yi = NULL;
    if ((x_[1] > x_[0] && x_[0] > xi) || (x_[1] < x_[0] && x_[0] < xi))
    {
        LOG(TRACE) << "input point is out of first bound " << endl;
        yi = y_[0];
    }
    else if ((x_.back() > x_.at(x_.size() - 2) && xi > x_.back()) || (x_.back() < x_.at(x_.size() - 2) && xi < x_.back()))
    {
        LOG(TRACE) << "input point is out of second bound " << endl;
        yi = y_.back();
    }
    else
    {
        int xSize = x_.size();
        for (int i = 0; i < xSize - 1; i++)
        {
            if (xi > x_[i] && xi < x_[i + 1])
            {
                yi = m_linearInterp_(x_[i], y_[i], x_[i + 1], y_[i + 1], xi);
                break;
            }
            else if (xi <= x_[i] && xi > x_[i + 1])
            {
                yi = m_linearInterp_(x_[i], y_[i], x_[i + 1], y_[i + 1], xi);
                break;
            }
            else if (xi == x_[i])
            {
                yi = y_[i];
                break;
            }
            else if (xi == x_[i+1])
            {
                yi = y_[i+1];
                break;
            }
        }
    }
    return yi;
}

double dataInterp::m_linearInterp_(double &x1, double &y1, double &x2, double &y2, double &xi)
{
    double yi = (xi - x1)*(y2 - y1) / (x2 - x1) + y1;
    return yi;
}

double dataInterp::m_getValue(double &xi, double &dyiDxi)
{
    double yi = NULL;
    if ((x_[1] > x_[0] && x_[0] > xi) || (x_[1] < x_[0] && x_[0] < xi))
    {
        cerr << "input point is out of first bound " << endl;
        yi = y_[0];
        dyiDxi = 0.0;
    }
    else if ((x_.back() > x_.at(x_.size() - 2) && xi > x_.back()) || (x_.back() < x_.at(x_.size() - 2) && xi < x_.back()))
    {
        cerr << "input point is out of second bound " << endl;
        yi = y_.back();
        dyiDxi = 0.0;
    }
    else
    {
        int xSize = x_.size();
        for (int i = 0; i < xSize - 1; i++)
        {
            if (xi >= x_[i] && xi < x_[i + 1])
            {
                yi = m_linearInterp_(x_[i], y_[i], x_[i + 1], y_[i + 1], xi, dyiDxi);
                break;
            }
            else if (xi <= x_[i] && xi > x_[i + 1])
            {
                yi = m_linearInterp_(x_[i], y_[i], x_[i + 1], y_[i + 1], xi, dyiDxi);
                break;
            }
        }
    }
    return yi;
}

double dataInterp::m_linearInterp_(double &x1, double &y1, double &x2, double &y2, double &xi, double &dyiDxi)
{
    double yi = (xi - x1)*(y2 - y1) / (x2 - x1) + y1;
    dyiDxi = (y2 - y1) / (x2 - x1);
    return yi;
}

//vector<double > stringstremToDouble()
//{
//
//}
