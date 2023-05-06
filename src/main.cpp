/**
 * @brief: 二维热传导
 * 
 * @birth: created by Dablelv on bql at 2023-05-05
 */
#include<iostream>
#include "tool/tool.h"
#include "matplot/matplot.h"
#include "math.h"

namespace plt = matplot;

class TDHeat
{
private:
    int x_N,
        y_N;
    double  W,
            H,
            k,
            q,
            T_top,
            dx,
            dy;
    vector_2d a;
    vector_2d bodyHeart;
    vector_1d T,b;
    vector_2d Z;
public:
    TDHeat();
    ~TDHeat(){};
    void input();
    void output();
    void construction();
    void plotImg();
};

TDHeat::TDHeat():x_N(3),y_N(4),W(0.3),H(0.4),k(1000),q(500000),T_top(100)
{

}

/**
 * @fname: input
 * @brief: 输入参数
 * @param: void
 * @return: void
 * @birth: created by Dablelv on bql
 */
void TDHeat::input()
{
    this->output();
    int isInput = 0;
    std::cout << "是否输入参数：";
    std::cin >> isInput;
    if(!isInput) return;
    std::cout << "宽度W = "; std::cin >> W;
    std::cout << "高度H = "; std::cin >> H;
    std::cout << "x方向网格数x_N = "; std::cin >> x_N;
    std::cout << "y方向网格数y_N = "; std::cin >> y_N;
    std::cout << "传热q = "; std::cin >> q;
    std::cout << "热导率k = "; std::cin >> k;
    std::cout << "顶部温度T_top = "; std::cin >> T_top;
    std::cout << "输入后参数为：";
    this->output();
}

/**
 * @fname: output
 * @brief: 输出参数
 * @param: void
 * @return: void
 * @birth: created by Dablelv on bql
 */
void TDHeat::output()
{
    std::cout << "\n默认参数如下:" << std::endl;
    std::cout << "\t宽度：" << W << std::endl;
    std::cout << "\t高度：" << H << std::endl;
    std::cout << "\tx方向网格数：" << x_N << std::endl;
    std::cout << "\ty方向网格数：" << y_N << std::endl;
    std::cout << "\t传热：" << q << std::endl;
    std::cout << "\t热导率：" << k << std::endl;
    std::cout << "\t顶部温度：" << T_top << std::endl;
}

/**
 * @fname: construction
 * @brief: 构造矩阵
 * @param: void
 * @return: void
 * @birth: created by Dablelv on bql
 */
void TDHeat::construction()
{
    dx = W / x_N;
    dy = H / y_N;

    // 第一类单元（内部单元5、6、9、10、13、14、17、18）
    double gDiff1_x = dy / dx;
    double gDiff1_y = dx / dy;
    double a1_x = - k * gDiff1_x;
    double a1_y = - k * gDiff1_y;
    double ac1 = - 2 * (a1_x + a1_y);

    // 第二类单元（绝热边1、2、7、11、15、19）
    double ac2_x = - (2 * a1_x + a1_y);
    double ac2_y = - (2 * a1_y + a1_x);

    // 第三类单元（恒温21、22）
    double gDiff3_y_wall = 2 * dx / dy;
    double a3_y_wall = - k * gDiff3_y_wall;
    double ac3 = - (2 * a1_x + a1_y ) - a3_y_wall;
    double Flux3_y_wall = a3_y_wall * T_top;
    double b3 = - Flux3_y_wall;

    // 第四类单元（恒热流4、8、12、16）
    double ac4 = -(2 * a1_y + a1_x);
    double Flux4_x_wall = - q * dy; 
    double b4 = -Flux4_x_wall;

    // 第五类单元（1恒热流+1绝热+2常规 0）
    double ac5 = -(a1_x + a1_y);
    double b5 = -Flux4_x_wall;

    // 第六类单元（2绝热+常规 3）
    double ac6 = -(a1_x + a1_y);
    
    // 第七类单元（1绝热+1恒温+2常规 23）
    double ac7 = -(a1_x + a1_y + a3_y_wall);
    double b7 = - Flux3_y_wall;

    // 第八类单元（1恒温+1恒热流 20）
    double ac8 = -(a1_x + a1_y + a3_y_wall);
    double b8 = -(Flux3_y_wall + Flux4_x_wall);

    // 初始化矩阵
    myTool::initVector_2d(a,x_N * y_N,y_N * x_N);
    myTool::printVector_2d(a,"初始化矩阵");

    // 开始组装矩阵
    // 左侧 i%x_N == 0
    // 右侧 (i+1)%x_N == 0
    // 上侧 floor(i/x_N) == (y_N-1)
    // 下侧 floor(i/x_N) == 0
    b.resize(x_N * y_N);
    for(int i = 0; i < x_N * y_N; ++i)
    {
        if((i%x_N != 0) && ((i+1)%x_N != 0) && (floor(i/x_N) != (y_N-1)) && (floor(i/x_N) != 0))
        {
            std::cout << "第一类单元：" << i << std::endl;
            a[i][i] = ac1;
            a[i][i-1] = a1_x;
            a[i][i+1] = a1_x;
            a[i][i+x_N] = a1_y;
            a[i][i-x_N] = a1_y;
            // a[i][x_N*y_N] = 0;
            b[i] = 0;
        }else if((floor(i/x_N) == 0) && (i!=0) && (i!=(x_N-1)))
        {
            std::cout << "第二类单元1：" << i << std::endl;
            a[i][i] = ac2_x;
            a[i][i-1] = a1_x;
            a[i][i+1] = a1_x;
            a[i][i+x_N] = a1_y;
            // a[i][x_N*y_N] = 0;
            b[i] = 0;
        }else if(((i+1)%x_N == 0) && (i!=(x_N-1)) && (floor(i/x_N) != (y_N-1)))
        {
            std::cout << "第二类单元2：" << i << std::endl;
            a[i][i] = ac2_y;
            a[i][i-1] = a1_x;
            a[i][i+x_N] = a1_y;
            a[i][i-x_N] = a1_y;
            a[i][x_N*y_N] = 0;
            b[i] = 0;
        }else if(floor(i/x_N) == (y_N-1) && (i%x_N != 0) && ((i+1)%x_N != 0))
        {
            std::cout << "第三类单元：" << i << std::endl;
            a[i][i] = ac3;
            a[i][i+1] = a1_x;
            a[i][i-1] = a1_x;
            a[i][i-x_N] = a1_y;
            // a[i][x_N*y_N] = b3;
            b[i] = b3;
        }else if((i%x_N == 0) && (floor(i/x_N) != 0) && (floor(i/x_N) != (y_N-1)))
        {
            std::cout << "第四类单元：" << i << std::endl;
            a[i][i] = ac4;
            a[i][i+1] = a1_x;
            a[i][i-x_N] = a1_y;
            a[i][i+x_N] = a1_y;
            // a[i][x_N*y_N] = b4;
            b[i] = b4;
        }else if((floor(i/x_N) == 0) && (i%x_N == 0))
        {
            std::cout << "第五类单元：" << i << std::endl;
            a[i][i] = ac5;
            a[i][i+1] = a1_x;
            a[i][i+x_N] = a1_y;
            // a[i][x_N*y_N] = b5;
            b[i] = b5;
        }else if(floor(i/x_N) == 0 && ((i+1)%x_N == 0))
        {
            std::cout << "第六类单元：" << i << std::endl;
            a[i][i] = ac6;
            a[i][i+x_N] = a1_y;
            a[i][i-1] = a1_x;
            b[i] = 0;
        }else if(((i+1)%x_N == 0) && (floor(i/x_N) == (y_N-1)))
        {
            std::cout << "第七类单元：" << i << std::endl;
            a[i][i] = ac7;
            a[i][i-1] = a1_x;
            a[i][i-x_N] = a1_y;
            // a[i][x_N*y_N] = b7;
            b[i] = b7;
        }else if((i%x_N == 0) && (floor(i/x_N) == (y_N-1)))
        {
            std::cout << "第八类单元：" << i << std::endl;
            a[i][i] = ac8;
            a[i][i+1] = a1_x;
            a[i][i-x_N] = a1_y;
            // a[i][x_N*y_N] = b8;
            b[i] = b8;
        }else
        {
            std::cout << "出错" << std::endl;
        }
    };
    myTool::printVector_2d(a,"填充后");
    // T = myTool::solveGuauss(a,x_N*y_N);
    T.resize(x_N * y_N);
    double eps = 1e-6;
    myTool::GaussSeidel(a,b,T,eps);
    myTool::printVector_1d(T);
    myTool::solveBodyHeart(bodyHeart,W,H,x_N,y_N);
}

void TDHeat::plotImg()
{
    // myTool::initVector_2d(Z,x_N,y_N);
    // myTool::printVector_2d(Z,"Z初始化");
    vector_1d x = plt::linspace(0,W,x_N);
    vector_1d y = plt::linspace(0,H,y_N);
    for(int i = 0; i < x_N ; ++i)
    {
        x[i] = bodyHeart[i][0];
    };
    for(int i = 0; i < y_N; ++i)
    {
        y[i] = bodyHeart[i*x_N][1];
    };
    auto [X,Y] = plt::meshgrid(x,y);

    Z.resize(y_N);
    for(int i = 0; i < y_N; ++i)
    {
        Z[i].resize(x_N);
    };
    int iindex;
    int jindex;
    for(int k = 0; k < x_N * y_N; ++k)
    {
        iindex = k % x_N;
        jindex = floor(k/x_N);
        Z[jindex][iindex] = T[k];
    };
    // int iindex;
    // int jindex;
    // for(int k = 0; k < x_N*y_N; ++k)
    // { 
    //     // for(int m = 0; m < y_N * x_N; ++m)
    //     // {
    //     //     for(int n = 0; n < y_N * x_N; ++n)
    //     //     {
    //     //         if((m<(k+1)*x_N) && (n<(k+1)*y_N))
    //     //         {
    //     //             Z[m][n] = T[k];
    //     //         }
    //     //     };
    //     // };
    //     iindex = k % x_N;
    //     jindex = floor(k/x_N);
    //     for(int i = jindex*x_N; i < (jindex+1)*x_N; ++i)
    //     {
    //         for(int j = iindex*y_N; j < (iindex+1)*y_N; ++j)
    //         {
    //             Z[i][j] = T[k];
    //         };
    //     };
    // };
    myTool::printVector_2d(Z,"温度矩阵");
    myTool::printVector_2d(X,"X轴");
    myTool::printVector_2d(Y,"Y轴");
    std::vector<double> tem = {10.0,20,30};
    plt::contourf(X,Y,Z,100);
    plt::colorbar()
        .limits({100,300});
        // .tick_values({30,90,150,210,270})
        // .ticklabels({"cold","cool","neutral","warm","hot"})
        // .label("温度等高线图");
    plt::show();
}

int main(){
    TDHeat tdheat;
    tdheat.input();
    tdheat.construction();
    tdheat.plotImg();
    return 0;
}