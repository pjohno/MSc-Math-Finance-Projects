#pragma once
#include <vector>

namespace MSC_PROJECTS
{
    
    int valueBermudanPutOption(const std::vector<double> &Si,std::vector<double> &Vi,double &theta,double X,double r,double sigma,double dT,int N,double xi,double tol);
    
}
