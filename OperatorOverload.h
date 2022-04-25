#include <vector>
#include <iostream>
//#include "forAllOperations.h"
//#include "Fields.h"
using std::vector;

// int Nx = 10, Ny = 5; // Nx = column, Ny = row
// int H = 0.01;
// int L = 10*H;      // L = 0.1
// double dx = L/Nx;
// double dy = H/dy;
// double dz = 0.01;

class Matrix
{
public:
    int Nx, Ny; // Nx = column, Ny = row
    double H, L;
    double dx, dy;

    // double dz = 0.01;

    Matrix();
    virtual ~Matrix();

    typedef vector<vector<Matrix>> Mat;         // finiteMat == two dimensional vector of the class "FiniteMatrix"
    typedef vector<Matrix> Mat1d; // +

    double value;

    void print1dmat(Mat1d&);    // +
    void print2dmat(Mat&);

    void setX(Mat1d&);
    void setY(Mat1d&);
    void setXC(Mat1d&, Mat1d&);
    void setYC(Mat1d&, Mat1d&);

    void inlet(Mat&, const double);
    void pressuredrop(Mat&, Mat1d&);
    void setpx(Mat1d&, const double);
    void setT(Mat&, const double, const double);
    void replace(Mat&, Mat&);
    void replaceSection(Mat&, Mat&, Mat1d&, Mat1d&);
    void replaceSameSection(Mat&, Mat&, Mat1d&, Mat1d&);
    void replaceline(Matrix::Mat&, const double, const int, const int);
    void replaceValue(Mat&, const double);

    double sumMatrix(Matrix::Mat&);
    void nodenumber(Mat1d&);
    void nodenumber2(Mat1d&);

    Matrix::Mat extract(Mat&, Mat1d&, Mat1d&);

    friend Matrix::Mat operator+(const Matrix::Mat&, const Matrix::Mat&);
    friend Matrix::Mat operator%(const double, const Matrix::Mat&);
    friend Matrix::Mat operator-(const Matrix::Mat&, const Matrix::Mat&);
    friend Matrix::Mat operator&&(const Matrix::Mat&, const Matrix::Mat&);
    friend Matrix::Mat operator*(const double, const Matrix::Mat&);
    friend Matrix::Mat operator/(const double, const Matrix::Mat& );

    friend Matrix::Mat1d operator+(const Matrix::Mat1d&, const Matrix::Mat1d&);
    friend Matrix::Mat1d operator%(const double, const Matrix::Mat1d&);
    friend Matrix::Mat1d operator-(const Matrix::Mat1d&, const Matrix::Mat1d&);
    friend Matrix::Mat1d operator&&(const Matrix::Mat1d&, const Matrix::Mat1d&);
    friend Matrix::Mat1d operator*(const double, const Matrix::Mat1d&);
    
};