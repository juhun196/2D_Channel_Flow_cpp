#include "OperatorOverload.h"
#include <iomanip>
#include <cmath>
using namespace std;

Matrix::Matrix()
:value(0.0),  Nx(200), Ny(40), H(0.01), L(0.1), dx(0.0005), dy(0.00025)
{
}

Matrix::~Matrix()
{
}

Matrix::Mat operator+(const Matrix::Mat& lhs, const Matrix::Mat& rhs)
{
    Matrix::Mat result(lhs);
    for(unsigned int i = 0; i < lhs.size(); i++)
    {
        for(unsigned int j = 0; j < lhs[i].size(); j++)
        {
            result[i][j].value += rhs[i][j].value;
        }
    }
    
    return result;
}

Matrix::Mat operator%(const double dblvalue, const Matrix::Mat& rhs)
{
    Matrix::Mat result(rhs);
    for(unsigned int i = 0; i < rhs.size(); i++)
    {
        for(unsigned int j = 0; j < rhs[i].size(); j++)
        {
            result[i][j].value += dblvalue;
        }
    }
    
    return result;
}

Matrix::Mat operator-(const Matrix::Mat& lhs, const Matrix::Mat& rhs)
{
    Matrix::Mat result(lhs);
    for(unsigned int i = 0; i < lhs.size(); i++)
    {
        for(unsigned int j = 0; j < lhs[i].size(); j++)
        {
            result[i][j].value -= rhs[i][j].value;
        }
    }
    
    return result;
}

Matrix::Mat operator&&(const Matrix::Mat& lhs, const Matrix::Mat& rhs)
{
    Matrix::Mat result(lhs);
    for(unsigned int i = 0; i < lhs.size(); i++)
    {
        for(unsigned int j = 0; j < lhs[i].size(); j++)
        {
            result[i][j].value *= rhs[i][j].value;
        }
    }
    
    return result;
}

Matrix::Mat operator*(const double dblvalue, const Matrix::Mat& rhs)
{
    Matrix::Mat result(rhs);
    for(unsigned int i = 0; i < rhs.size(); i++)
    {
        for(unsigned int j = 0; j < rhs[i].size(); j++)
        {
            result[i][j].value *= dblvalue;
        }
    }
    
    return result;
}

Matrix::Mat operator/(const double dblvalue, const Matrix::Mat& rhs)
{
    Matrix::Mat result(rhs);
    for(unsigned int i = 0; i < rhs.size(); i++)
    {
        for(unsigned int j = 0; j < rhs[i].size(); j++)
        {
            result[i][j].value = dblvalue / rhs[i][j].value;
        }
    }
    
    return result;
}

Matrix::Mat1d operator+(const Matrix::Mat1d& lhs, const Matrix::Mat1d& rhs)
{
    Matrix::Mat1d result(lhs);
    for(unsigned int i = 0; i < lhs.size(); i++)
    {
        result[i].value += rhs[i].value;
    }
    
    return result;
}

Matrix::Mat1d operator%(const double dblvalue, const Matrix::Mat1d& rhs)
{
    Matrix::Mat1d result(rhs);
    for(unsigned int i = 0; i < rhs.size(); i++)
    {
        result[i].value += dblvalue;
    }
    
    return result;
}


Matrix::Mat1d operator-(const Matrix::Mat1d& lhs, const Matrix::Mat1d& rhs)
{
    Matrix::Mat1d result(lhs);
    for(unsigned int i = 0; i < lhs.size(); i++)
    {
        result[i].value -= rhs[i].value;
    }
    
    return result;
}


Matrix::Mat1d operator&&(const Matrix::Mat1d& lhs, const Matrix::Mat1d& rhs)
{
    Matrix::Mat1d result(lhs);
    for(unsigned int i = 0; i < lhs.size(); i++)
    {
        result[i].value *= rhs[i].value;
    }
    
    return result;
}


Matrix::Mat1d operator*(const double dblvalue, const Matrix::Mat1d& rhs)
{
    Matrix::Mat1d result(rhs);
    for(unsigned int i = 0; i < rhs.size(); i++)
    {
        result[i].value *= dblvalue;
    }
    
    return result;
}

double sumMatrixAbs(Matrix::Mat& vec)
{
    double result = 0;

    for(unsigned int i = 0; i < vec.size(); i++)
    {
        for(unsigned int j = 0; j < vec[0].size(); j++)
        {
            result += abs(vec[i][j].value);
        }
    }
    
    return result;
}



void Matrix::print2dmat(Mat& vec)
{
    cout << "size : " << vec.size() << "x" <<  vec[0].size() << endl;
    for(unsigned int i = 0; i < vec.size(); i++)
    {
        for(unsigned int j = 0; j < vec[i].size(); j++)
        {
            std::cout << std::setprecision(3)<< vec[i][j].value << ' ';
        }
        cout << endl;
    }
}

void Matrix::print1dmat(Mat1d& vec)     // +
{   cout << "size : " << vec.size() << endl;
    for(unsigned int i = 0; i < vec.size(); i++)
    {
        std::cout << std::setprecision(3) << vec[i].value << ' ';
    }
    cout << endl;
}

void Matrix::setX(Mat1d& vecX)
{
    for(int i =1; i < vecX.size(); i++)
    {
        vecX[i].value = vecX[i-1].value + dx;
    }
    vecX[vecX.size()-1] = vecX[vecX.size()-2]; 
}

void Matrix::setY(Mat1d& vecY)
{
    for(int i =1; i < vecY.size(); i++)
    {
        vecY[i].value = vecY[i-1].value + dy;
    }
    vecY[vecY.size()-1] = vecY[vecY.size()-2];   
}

void Matrix::setXC(Mat1d& vecX, Mat1d& vecXC)
{
    for(int i =1; i < vecX.size(); i++)
    {
        vecXC[i].value = (vecX[i].value + vecX[i-1].value)*0.5;
    }
    vecXC[vecX.size()-1] = vecX[vecX.size()-2];
}

void Matrix::setYC(Mat1d& vecY, Mat1d& vecYC)
{
    for(int i =1; i < vecY.size(); i++)
    {
        vecYC[i].value = (vecY[i].value + vecY[i-1].value)*0.5;
    }
    vecYC[vecY.size()-1] = vecY[vecY.size()-2];
}

void Matrix::inlet(Mat& vec, const double U)
{
    for(int i = 0; i < vec.size(); i++)
    {
        for(int j = 1; j < vec[i].size()-1; j++)
        {
            vec[i][j].value = U;
        }
    }
}

void Matrix::pressuredrop(Mat& vec, Mat1d& vecX)
{
    for(int i = 1; i < vec.size()-2; i++)
    {
        for(int j = 1; j < vec[0].size()-1; j++)
        {
            vec[i][j].value = vecX[i].value;
        }
    }
}

void Matrix::setpx(Mat1d& vecX, const double p1)
{
    vecX[1].value = p1;
    for(int i = 2; i < vecX.size(); i++)
    {
        vecX[i].value = vecX[i-1].value - p1/(Nx-1);
    }
}

void Matrix::setT(Mat& vec, const double Ti, const double Tw)
{
    for(int i = 0; i < vec.size(); i++)
    {   
        for(int j = 0; j< vec[i].size(); j++)
        {   
            if(j == 0)
                vec[i][j].value = Tw;
            else if(j == vec[i].size()-1)
                vec[i][j].value = Tw;
            else
                vec[i][j].value = Ti;
        }
    }
}

void Matrix::replace(Mat& vecNew, Mat& vecOld)
{
    for(int i = 0; i < vecOld.size(); i++)
    {
        for(int j = 0; j < vecOld[i].size(); j++)
        {
            vecNew[i][j].value = vecOld[i][j].value;
        }
    }
}

void Matrix::replaceValue(Mat& vec, const double value)
{
    for(int i = 0; i < vec.size(); i++)
    {
        for(int j = 0; j < vec[i].size(); j++)
        {
            vec[i][j].value = value;
        }
    }
}

void Matrix::replaceSection(Matrix::Mat& vecOld, Matrix::Mat& vecNew, Matrix::Mat1d& vecX, Matrix::Mat1d& vecY)
{
    for(int i = 0; i < vecNew.size(); i++)
    {
        for(int j = 0; j < vecNew[i].size(); j++)
        {
            vecOld[vecX[i].value][vecY[j].value].value = vecNew[i][j].value;
        }
    }
}

void Matrix::replaceSameSection(Matrix::Mat& vecOld, Matrix::Mat& vecNew, Matrix::Mat1d& vecX, Matrix::Mat1d& vecY)
{
    for(int i = vecX[0].value; i <= vecX[vecX.size()-1].value; i++)
    {
        for(int j = vecY[0].value; j <= vecY[vecY.size()-1].value; j++)
        {
            vecOld[i][j].value = vecNew[i][j].value;
        }
    }
}

// Matrix::Mat replaceSection(Matrix::Mat& vecOld, Matrix::Mat& vecNew, Matrix::Mat1d& vecX, Matrix::Mat1d& vecY)
// {
//     Matrix::Mat result(vecOld);
//     for(int i = 0; i < vecNew.size(); i++)
//     {
//         for(int j = 0; j < vecNew[i].size(); j++)
//         {
//             vecOld[vecX[i].value][vecY[j].value].value = vecNew[i][j].value;
//         }
//     }

//     return result;
// }

void Matrix::replaceline(Matrix::Mat& vec, const double value, const int line, const int rowOrcol)
{  
    if(rowOrcol == 1)   // row
    {
        for(int j = 0; j < vec[0].size(); j++)
        {
            vec[line][j].value = value;
        }
    }

    else               // column
    {
        for(int i = 0; i < vec.size(); i++)
        {
            vec[i][line].value = value;
        }
    }
    
}

Matrix::Mat extract(Matrix::Mat& vec, Matrix::Mat1d& vecX, Matrix::Mat1d& vecY)
{
    Matrix::Mat result(vecX.size(), vecY);
    for(int i = 0; i < vecX.size(); i++)
    {
        for(int j = 0; j < vecY.size(); j++)
        {
            result[i][j].value = vec[vecX[i].value][vecY[j].value].value;
        }
    }

    return result;
}

void Matrix::nodenumber(Mat1d& vecX)
{
    for(int i = 0; i < vecX.size(); i++)
    {   
        vecX[i].value = i;
    }
}

void Matrix::nodenumber2(Mat1d& vecX)
{
    for(int i = 0; i < vecX.size(); i++)
    {   
        vecX[i].value = i+1;
    }
}


// void Matrix::setY(Mat1d& vec)
// {

// }