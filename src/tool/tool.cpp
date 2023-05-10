#include "tool.h"
#include "math.h"

namespace myTool
{
    /**
     * @fname: initVector_2d
     * @brief: 初始化二维矩阵
     * @param: vector_2d &a,const int x,const int y
     * @return: void
     * @birth: created by Dablelv on bql
     */
    void initVector_2d(vector_2d &a,const int x,const int y)
    {
        a.resize(x);
        for(int i = 0; i < y; ++i)
        {
            a[i].resize(y);
        };
    }

    /**
     * @fname: printVector_2d
     * @brief: 打印二维矩阵
     * @param: vector_2d &a,const std::string msg
     * @return: void
     * @birth: created by Dablelv on bql
     */
    void printVector_2d(vector_2d &a,const std::string msg)
    {
        std::cout << "\n" << msg << "矩阵为" << std::endl;
        for(int i = 0; i < a.size(); ++i)
        {
            for(int j = 0; j < a[i].size(); ++j)
            {
                std::cout << "\t" << a[i][j];
            };
            std::cout << std::endl;
        };
    }

    /**
     * @fname: solveGuauss
     * @brief: 使用高斯列主元求解矩阵
     * @param: vector_2d &a,const int n
     * @return: vector_1d
     * @birth: created by Dablelv on bql
     */
    vector_1d solveGuauss(vector_2d &a,vector_1d &b,const int n)
    {
        for(int k = 0;k< n-1;k++){
            // 选择主元
            double max = fabs(a[k][k]);
            int pivrow = k;
            double ratio;
            for(int i=k+1;i<n;i++)
            {
                if(fabs(a[i][k]) > max)
                {
                    max = fabs(a[i][k]);
                    pivrow = i;
                    // swap
                    swap(a,b,k,pivrow,n);
                }
                // 消元
                ratio = a[i][k] / a[k][k];
                for(int j = k;j<n;j++)
                {
                    a[i][j] -= ratio * a[k][j];
                }
                b[i] -= ratio * b[k];
                a[i][k] = 0;
            }
        }
//        myTool::printVector_2d(a,"消元后");
//        myTool::printVector_1d(b,"b消元后");
        
        vector_1d T(n,0);
        T[n-1] = b[n-1] / a[n-1][n-1];
        double sum;
        for(int i = n-2;i>=0;i--){
            sum = 0.0;
            for(int j = i+1;j<n;j++){
                sum += a[i][j] * T[j];
            }
            T[i] = (b[i] - sum) / a[i][i];
        }
        return T;
    }

    /**
     * @fname: swap
     * @brief: 交换矩阵的行
     * @param: vector_2d &a,const int k,const int pivrow,const int n
     * @return: void
     * @birth: created by Dablelv on bql
     */
    void swap(vector_2d &a,vector_1d &b,const int k,const int pivrow,const int n)
    {
        double temp;
        temp = b[k];
        b[k] = b[pivrow];
        b[pivrow] = temp;

        for(int j = k;j<n;j++){
            temp = a[k][j];
            a[k][j] = a[pivrow][j];
            a[pivrow][j] = temp;
        }
    }

    /**
     * @fname: printVector_1d
     * @brief: 打印一维矩阵
     * @param: vector_1d &T
     * @return: void
     * @birth: created by Dablelv on bql
     */
    void printVector_1d(vector_1d &T,const std::string msg)
    {
        std::cout << "\n" <<  msg << "矩阵为：" << std::endl;
        for(int i = 0; i < T.size(); ++i)
        {
            std::cout << T[i] << "\t";
        };
    }

    /**
     * @fname: solveBodyHeart
     * @brief: 求解每个网格单元的体心
     * @param: vector_2d &bodyHeart,const double  W,const double H,const int x_N,const int y_N
     * @return: void
     * @birth: created by Dablelv on bql
     */
    void solveBodyHeart(vector_2d &bodyHeart,const double  W,const double H,const int x_N,const int y_N)
    {
        double dx = W / x_N;
        double dy = H / y_N;
        bodyHeart.resize(x_N * y_N);
        for(int i = 0; i < bodyHeart.size(); ++i)
        {
            bodyHeart[i].resize(2);
        };
        int iindex;
        int jindex;
        for(int i = 0; i < bodyHeart.size(); ++i)
        {
            iindex = i % x_N;
            jindex = floor(i/x_N);
            bodyHeart[i][0] = (dx / 2) + iindex * dx;
            bodyHeart[i][1] = (dy / 2) + jindex * dy;
        };
        std::cout << std::endl;
        for(int i = 0; i < x_N * y_N; ++i)
        {
            std::cout << "(" << bodyHeart[i][0] << "," << bodyHeart[i][1] << ")" << std::endl;
        };
    }

    void GaussSeidel(vector_2d &a,vector_1d &b,vector_1d &x, const double eps)
    {
        int flag = 0,
            iter = 0,
            n = a.size();
        double xold = 0.0,
                sum = 0.0,
                error = 0.0;
        do
        {
            flag = 0;
            iter ++;
            std::cout << "第" << iter << "次迭代" << std::endl;
            for(int i = 0; i < n; ++i)
            {
                xold = x[i];
                sum = 0.0;
                for(int j = 0; j < n; ++j)
                {
                    if(j != i) sum += a[i][j] * x[j];  
                };
                x[i] = (b[i] - sum) / a[i][i];
                error = fabs(xold - x[i]) / x[i];
                if(error >= eps){ flag = 1; }
            };
        }while(flag == 1);
        std::cout << "迭代完成" << std::endl;
    }

    void RelaxationGaussSeidel(vector_2d &a,vector_1d &b, vector_1d &x,const double eps,const double w)
    {
        int flag,
                iter = 0;
        double error,
                sum;
        int n = a.size();
        vector_1d xold(n,0);
        do
        {
            flag = 0;
            iter++;
            std::cout << "\n第" << iter << "次迭代" << std::endl;
            for(int i = 0; i < n; ++i)
            {
                sum = 0.0;
                xold[i] = x[i];
                for(int j = 0; j < n; ++j)
                {
                    if(j != i) {sum += a[i][j] * x[j];}
                };
                x[i] = (b[i] - sum) / a[i][i];
            };
            for(int i = 0; i < n; ++i)
            {
                x[i] = w * x[i] + (1 - w) * xold[i];
            };
            for(int i = 0; i < n; ++i)
            {
                error = fabs((xold[i] - x[i]) / x[i]);
                (error >= eps) ? (flag = 1) : (flag = 0);
            };

        }while(flag == 1);
    }
}