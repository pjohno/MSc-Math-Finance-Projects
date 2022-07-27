#pragma once

#include <vector>

namespace MSC_PROJECTS
{
    int solveStorageOption(int M,int N,double Qmax,double Xmax,std::vector<double> &Q,std::vector<double> &X,
                           double r,double p,double sigma,std::vector<std::vector<double>> &v,int qMax,double tol);
}
