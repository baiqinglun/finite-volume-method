#ifndef _TOOL_H
#define _TOOL_H

#include<iostream>
#include<vector>
using vector_1d = std::vector<double>;
using vector_2d = std::vector<vector_1d>;

namespace myTool
{
    void initVector_2d(vector_2d &a,const int x,const int y);
    void printVector_2d(vector_2d &a,const std::string msg);
    vector_1d solveGuauss(vector_2d &a,const int n);
    void swap(vector_2d &a,const int k,const int pivrow,const int n);
    void printVector_1d(vector_1d &a);
    void solveBodyHeart(vector_2d &bodyHeart,const double  W,const double H,const int x_N,const int y_N);
}

#endif