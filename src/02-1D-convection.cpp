/********************************************************************************
* @author: Bqi Qinglun
* @date: 2023/5/10 14:19
* @version: 1.0
* @description: 一维对流问题：使用中心差分格式求解对流项
********************************************************************************/
#include<iostream>
#include "tool/tool.h"
#include "matplot/matplot.h"

namespace plt = matplot;

class ODConvection
{
private:
    int n;
    double eps,
            L,
            u,
            dL,
            gamma,
            rho,
            dy,
            phia,
            phib,
            w;
    vector_1d b,x;
    vector_2d a;

public:
    ODConvection();
    ~ODConvection();
    void solution();
    void plotImg();
};

ODConvection::ODConvection() : n(50),eps(1e-6),L(1.0),u(2.5),gamma(0.1),rho(1),dy(1.0),phia(1.0),phib(0),w(1.5)
{
    dL = L / n;
}

ODConvection::~ODConvection()
{

}

/********************************************************************************
 * @fname: solution
 * @brief: 求解一维对流问题
 * @param: void
 * @return: void
 * @birth: created by Dablelv on bql
********************************************************************************/
void ODConvection::solution()
{
    // 初始矩阵
    myTool::initVector_2d(a,n,n);
    b.resize(n);

    for (int i = 0; i < n; ++i) {
        if (0 == i)
        {
            a[i][0] = rho * u / 2 + gamma / dL + 2 * gamma / dL;
            a[i][1] = -gamma/ dL + rho * u / 2;
            b[i] = (rho * u * phia) + 2 * gamma / dL * phia;
        }
        else if(i == (n - 1))
        {
            a[i][i] = 2 * gamma / dL + (gamma / dL - rho * u / 2);
            a[i][i-1] = -gamma / dL - rho * u / 2;
            b[i] = 2 * gamma / dL * phib - rho * u * phib;
        }
        else
        {
            a[i][i-1] = - gamma / dL - rho * u * dy / 2;
            a[i][i+1] = - gamma / dL + rho * u * dy / 2;
            a[i][i] = 2 * gamma / dL;
        }
    }
//    myTool::printVector_2d(a,"组装后");
//    myTool::printVector_1d(b,"b");

    x.resize(n);
    // RelaxationGaussSeidel和GaussSeidel适用于对角占优矩阵，用在本案例不合适。本案例使用高斯消去法直接求解
//    myTool::RelaxationGaussSeidel(a,b,x,eps,w);
//    myTool::GaussSeidel(a,b,x,eps);
    x = myTool::solveGuauss(a,b,n);
//    myTool::printVector_1d(x,"x");
}


/********************************************************************************
 * @fname: plotImg
 * @brief: 使用scatter和plot在一张图上同时绘制，需要使用hold()函数
 * @param: void
 * @return: void
 * @birth: created by Dablelv on bql
********************************************************************************/
void ODConvection::plotImg()
{
    vector_1d X = plt::linspace(0,L,n);
//    vector_1d X1 = plt::linspace(0,L,200);
    vector_1d Y(n,0);
    for (int i = 0; i < n; ++i) {
        Y[i] = (exp(rho * u * X[i] / gamma) -1) / (exp(rho * u * L / gamma) -1) * (phib - phia) + phia;
    }
    auto f = plt::figure();     // 创建画布
    auto ax = f->current_axes();    // 创建坐标系
    // 加上这两个有问题，图有时显示有时不显示
//    ax->title("Contrast Numerical and Analysis");
//    ax->legend({"Numerical ","Analysis"});
    ax->scatter(X,x);
    ax->hold(true);
    ax->plot(X,Y);
    ax->hold(false);
    plt::show();
}

int main()
{
    ODConvection odConvection;
    odConvection.solution();
    odConvection.plotImg();
    return 0;
}