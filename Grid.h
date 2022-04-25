#pragma once
#include <vector>

int Nx = 200, Ny = 40; // Nx = column, Ny = row
double H = 0.01;
double L = 10*H;      // L = 0.1
double dx = L/Nx;     // 0.0005
double dy = H/Ny;     // 0.00025
double dz = 0.01;


// linspace
// x  = dx/2:dx:L-dx/2;    % x-locations of main grid points (m)
// xu = 0:dx:L;            % x-locations of u-velocities (m)
// y  = dy/2:dy:H-dy/2;    % y-locations of main grid points (m)
// yv = 0:dy:H;            % y-locations of v-velocities (m)

// iu = 1:Nx+1; Ju = 2:Ny+1;   % Interior node numbers for u
// Iv = 2:Nx+1; jv = 1:Ny+1;   % Interior node numbers for v
// Ip = 2:Nx+1; Jp = 2:Ny+1;   % Interior node numbers for p
// iF = 2:Nx;   jF = 2:Ny;     % Node numbers for face flow (or advection)



    // Properties (air at STP)
double rho   = 1.2;       // Density (kg/m^3)
double mu    = 1.8e-5;    // Absolute viscosity (N-s/m^2)
double nu    = mu/rho;    // Kinematic viscosity (m^2/s)
double kt    = 0.025;     // Thermal conductivity (W/m-k)
double cp    = 1006;      // specific heat (J/kg-K)
double alpha = kt/(rho*cp); // Theramal diffusivity (m^2/s)
double Pr    = nu/alpha;     // Prandtl number

// Boundary conditions
int     Re   = 100;          // Reynolds number
double  U    = Re*nu/(2*H);  // Average velocity (m/s)
int     Ti   = 20;           // Inlet temperature (deg. C)
int     Tw   = 100;          // Wall temperature (deg. C)
int     qw   = 100;          // Wall heat flux (W/m^2)

int     BC_N = 1;            // BC_N = 0 for Tw, BC_N = 1 for qw
int     BC_S = 1;            // BC_S = 0 for Tw, BC_S = 1 for symmetry

// Solution controls
double   alphaU = 0.3;        // Velocity relaxation (under)
double   alphaP = 0.2;        // Pressure relaxation (under)
int     NmaxSIM = 1e+4;       // Iteration max for SIMPLE algorithm (-)
int     NmaxGSI = 1e+1;       // Iteration max for numerical method (-)
double    err   = 1e-5;       // Convergence criteria (-)
int       div   = 1e+1;       // Divergence criteria (-)

//Initialize to linear pressure drop for fully developed flow
double p1 = 12*mu*U*L/((2*H)*(2*H));  

double Dx = (mu/dx)*dy*dz;
double Dy = (mu/dy)*dx*dz;

double res_sum;
double res2_sum;
double res3_sum;
double res4_sum;
double cTest;
double Tres;

Matrix finiteObj;

Matrix::Mat1d X(Nx);
Matrix::Mat1d XC(Nx+2);
Matrix::Mat1d Y(Ny);
Matrix::Mat1d YC(Ny+2);
Matrix::Mat1d xu(Nx+2);
Matrix::Mat1d yv(Ny+2);

Matrix::Mat1d px(Nx);

Matrix::Mat1d oneToNx_1(Nx-1);
Matrix::Mat1d oneToNy_1(Ny-1);
Matrix::Mat1d oneToNx_2(Nx);
Matrix::Mat1d oneToNy_2(Ny);

Matrix::Mat1d iF(Nx-1);
Matrix::Mat1d jF(Ny-1);
Matrix::Mat1d Ip(Nx);
Matrix::Mat1d Jp(Ny);
Matrix::Mat1d Iv(Nx);
Matrix::Mat1d jv(Ny+1);
Matrix::Mat1d iu(Nx+1);
Matrix::Mat1d Ju(Ny);


// Initialize u, v, p, T, F, a, and residual matirces
Matrix::Mat u(Nx+1, vector<Matrix>(Ny+2));
Matrix::Mat uOld(Nx+1, vector<Matrix>(Ny+2));
Matrix::Mat v(Nx+2, vector<Matrix>(Ny+1));
Matrix::Mat vOld(Nx+2, vector<Matrix>(Ny+1));
Matrix::Mat uStar(Nx+1, vector<Matrix>(Ny+2));
Matrix::Mat vStar(Nx+2, vector<Matrix>(Ny+1));
Matrix::Mat uPrime(Nx+1, vector<Matrix>(Ny+2));
Matrix::Mat vPrime(Nx+2, vector<Matrix>(Ny+1));
Matrix::Mat dU(Nx+1, vector<Matrix>(Ny+2));
Matrix::Mat dV(Nx+2, vector<Matrix>(Ny+1));

Matrix::Mat T(Nx+2, vector<Matrix>(Ny+2));
Matrix::Mat T_ones(Nx+2, vector<Matrix>(Ny+2));
Matrix::Mat p(Nx+2, vector<Matrix>(Ny+2));
Matrix::Mat pStar(Nx+2, vector<Matrix>(Ny+2));
Matrix::Mat pPrime(Nx+2, vector<Matrix>(Ny+2));

Matrix::Mat Fe(Nx+1, vector<Matrix>(Ny+1));
Matrix::Mat Fw(Nx+1, vector<Matrix>(Ny+1));
Matrix::Mat Fn(Nx+1, vector<Matrix>(Ny+1));
Matrix::Mat Fs(Nx+1, vector<Matrix>(Ny+1));
Matrix::Mat DF(Nx+1, vector<Matrix>(Ny+1));
Matrix::Mat aE(Nx+1, vector<Matrix>(Ny+1));
Matrix::Mat aW(Nx+1, vector<Matrix>(Ny+1));
Matrix::Mat aN(Nx+1, vector<Matrix>(Ny+1));
Matrix::Mat aS(Nx+1, vector<Matrix>(Ny+1));
Matrix::Mat aP(Nx+1, vector<Matrix>(Ny+1));
Matrix::Mat bP(Nx+1, vector<Matrix>(Ny+1));



Matrix::Mat ures(NmaxSIM, vector<Matrix>(1));
Matrix::Mat vres(NmaxSIM, vector<Matrix>(1));
Matrix::Mat pres(NmaxSIM, vector<Matrix>(1));






// void setX(vector<double> vecX);

// Matrix finiteObj;
// Initialize u,v, p, T, F, a, and residual matrices


// Matrix::Mat AE(Ny, vector<Matrix>(Nx));
// Matrix::Mat AE2(Ny, vector<Matrix>(Nx));
// Matrix::Mat U_result(9.0*AE2);

// finiteObj.print2dmat(U_result);


