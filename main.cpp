#include "OperatorOverload.cpp"
#include <iostream>
#include <ctime>
using namespace std;


int main()
{
    clock_t start, finish;
    double duration;

    start = clock();
    #include "Grid.h"    

    //finiteObj.print2dmat(u);
    //Matrix::Mat U_result(1.0 % u);
    //finiteObj.print2dmat(U_result);

    // finiteObj.setX(xu);
    // finiteObj.print1dmat(xu);
    // finiteObj.setY(yv);
    // finiteObj.print1dmat(yv);
    

    // Inlet velocity for uniform flow at inlet
    finiteObj.inlet(u, U);
    finiteObj.setpx(px, p1);
    finiteObj.pressuredrop(p, px);
    finiteObj.setT(T, Ti, Tw);
    finiteObj.nodenumber2(iF);
    finiteObj.nodenumber2(Ju);
    finiteObj.nodenumber2(oneToNx_1);
    finiteObj.nodenumber2(oneToNy_1);

    finiteObj.nodenumber2(oneToNx_2);
    finiteObj.nodenumber2(oneToNy_2);


    //finiteObj.print1dmat(oneToNx_1);

    // /////////////////////// 나중에 지워야 함!!!!!!/////////////////////////////////////////
    // finiteObj.replace(u, uOld);
    // finiteObj.replace(v, vOld);
    // finiteObj.replace(p, pStar);
    // /////////////////////// 나중에 지워야 함!!!!!!/////////////////////////////////////////


    // /////////////////////////  FVM_u  ///////////////////////////////////////////////////
    // Fw
    

    // NmaxSIM
    for(int n = 0; n < NmaxSIM; n++)
    {
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // STEP 1a : solve x-momentum as uStar

        // Setup coefficients
        finiteObj.replace(uOld, u);
        finiteObj.replace(vOld, v);
        finiteObj.replace(pStar, p);

        Matrix::Mat u1 = extract(u, iF, Ju);
        Matrix::Mat1d iF1(1.0 % iF);
        Matrix::Mat1d iF2(-1.0 % iF);
        Matrix::Mat u2 = extract(u, iF2, Ju);

        Matrix::Mat u3(rho*dy*dz*0.5*(u1+u2));
        finiteObj.replaceSection(Fw, u3, iF, Ju);
        // Fe
        u1 = extract(u, iF, Ju);
        u2 = extract(u, iF1, Ju);
        u3 = (rho*dy*dz*0.5*(u1+u2));
        finiteObj.replaceSection(Fe, u3, iF, Ju);

        // Fs
        Matrix::Mat1d Ju2(-1.0 % Ju);
        Matrix::Mat v1 = extract(v, iF, Ju2);
        Matrix::Mat v2 = extract(v, iF1, Ju2);
        Matrix::Mat v3 = (rho*dy*dz*0.5*(v1+v2));
        finiteObj.replaceSection(Fs, v3, iF, Ju);

        //Fn
        v1 = extract(v, iF, Ju);
        v2 = extract(v, iF1, Ju);
        v3 = (rho*dy*dz*0.5*(v1+v2));
        finiteObj.replaceSection(Fn, v3, iF, Ju);
        
        for(int i = iF[0].value; i <= iF[iF.size()-1].value; i++)
        {
            for(int j = Ju[0].value; j <= Ju[Ju.size()-1].value; j++)
            {
                aW[i][j].value = max(0.0, max(Fw[i][j].value, Dx+0.5*Fw[i][j].value));
                aE[i][j].value = max(0.0, max(-Fe[i][j].value, Dx-0.5*Fe[i][j].value));
                aS[i][j].value = max(0.0, max(Fs[i][j].value, Dy+0.5*Fs[i][j].value));
                aN[i][j].value = max(0.0, max(-Fn[i][j].value, Dy-0.5*Fn[i][j].value));
            }
        }

        if(BC_S == 0)   // Wall or symmetry at bottom
        {
            finiteObj.replaceline(aS, 2*Dy, 1, 0);
        }
        else
        {
            finiteObj.replaceline(aS, 0.0, 1, 0);
        }


        finiteObj.replaceline(aN, 2*Dy, Ny, 0);     // Wall at top
        finiteObj.replaceline(aE, 1.0, Nx, 1);     // Right boundary (outlet, du/dx = 0)

        Matrix::Mat DF_prev = Fe - Fw + Fn - Fs;
        finiteObj.replaceSection(DF, DF_prev, iF, Ju);

        Matrix::Mat aP_prev = aE + aW + aN + aS + DF;
        finiteObj.replaceSameSection(aP, aP_prev, iF, Ju);

        Matrix::Mat pStar1 = extract(pStar, iF, Ju);
        Matrix::Mat pStar2 = extract(pStar, iF1, Ju);
        Matrix::Mat bP_prev = dy*(dz*(pStar1 - pStar2));
        finiteObj.replaceSection(bP, bP_prev, iF, Ju);


        Matrix::Mat dU_prev = (dy*dz)/((1/alphaU)*aP);
        finiteObj.replaceSameSection(dU, dU_prev, iF, Ju);


        // Use previous calculation as initial guess in numerical method

        Matrix::Mat res(Nx-1, vector<Matrix>(Ny));
        finiteObj.replaceValue(res, 0.0);

        for(int k = 0; k < NmaxGSI; k++)
        {
            for(int j = 1; j < Ny+1; j++)         // Ny = 5
            {
                for(int i = 1; i < Nx; i++)     // Nx = 10
                {
                    uOld[i][j].value = (aW[i][j].value*uOld[i-1][j].value + aE[i][j].value*uOld[i+1][j].value + aS[i][j].value*uOld[i][j-1].value +
                                        aN[i][j].value*uOld[i][j+1].value + bP[i][j].value)/(aP[i][j].value/alphaU) + (1-alphaU)*uOld[i][j].value;
                }
                
            }
            
            for(int j = 1; j < Ny+1; j++)
            {
                for(int i = 1; i < Nx; i++)
                {
                    res[i-1][j-1].value = (aW[i][j].value*uOld[i-1][j].value + aE[i][j].value*uOld[i+1][j].value + aS[i][j].value*uOld[i][j-1].value + 
                                            aN[i][j].value*uOld[i][j+1].value + bP[i][j].value) - aP[i][j].value*uOld[i][j].value;
                }
            }            

            Matrix::Mat aP_ex = extract(aP, oneToNx_1, oneToNy_2)&&extract(uOld, oneToNx_1, oneToNy_2);
            res_sum = sumMatrixAbs(res) / sumMatrixAbs(aP_ex);
            
            if(res_sum < err)
            {
                break;
            }
        }

        finiteObj.replace(uStar, uOld);
        ures[n][0].value = res_sum;

        //cout << ures[0].size() << endl;

        // STEP 1b: solve y-momentum as vStar
        // Set up coefficients, FVM_v
        finiteObj.nodenumber2(Iv);
        finiteObj.nodenumber2(jF);
        
        Matrix::Mat1d Iv1(1.0 % Iv);
        Matrix::Mat1d jF1(1.0 % jF);
        Matrix::Mat1d Iv2(-1.0 % Iv);
        Matrix::Mat1d jF2(-1.0 % jF);

        u1 = extract(u, Iv2, jF);
        u2 = extract(u, Iv2, jF1);
        u3 = extract(u, Iv, jF);
        Matrix::Mat u4 = extract(u, Iv, jF1);

        v1 = extract(v, Iv, jF);
        v2 = extract(v, Iv, jF2);
        v3 = extract(v, Iv, jF1);

        Matrix::Mat Fw_prev = rho*dy*dz*0.5*(u1 + u2);
        Matrix::Mat Fe_prev = rho*dy*dz*0.5*(u3 + u4);
        Matrix::Mat Fs_prev = rho*dy*dz*0.5*(v1 + v2);
        Matrix::Mat Fn_prev = rho*dy*dz*0.5*(v1 + v3);

        finiteObj.replaceSection(Fw, Fw_prev, Iv, jF);
        finiteObj.replaceSection(Fe, Fe_prev, Iv, jF);
        finiteObj.replaceSection(Fs, Fs_prev, Iv, jF);
        finiteObj.replaceSection(Fn, Fn_prev, Iv, jF);
     
        // Convective flux using hybrid upwinding
        for(int i = Iv[0].value; i <= Iv[Iv.size()-1].value; i++)
        {
            for(int j = jF[0].value; j <= jF[jF.size()-1].value; j++)
            {
                aW[i][j].value = max(0.0, max(Fw[i][j].value, Dx+0.5*Fw[i][j].value));
                aE[i][j].value = max(0.0, max(-Fe[i][j].value, Dx-0.5*Fe[i][j].value));
                aS[i][j].value = max(0.0, max(Fs[i][j].value, Dy+0.5*Fs[i][j].value));
                aN[i][j].value = max(0.0, max(-Fn[i][j].value, Dy-0.5*Fn[i][j].value));
            }
        }

        finiteObj.replaceline(aE, 0.0, Nx, 1);
                
        DF_prev = Fe - Fw + Fn - Fs;
        finiteObj.replaceSameSection(DF, DF_prev, Iv, jF);

        aP_prev = aE + aW + aN + aS + DF;
        finiteObj.replaceSameSection(aP, aP_prev, Iv, jF);

        pStar1 = extract(pStar, Iv, jF);
        pStar2 = extract(pStar, Iv, iF1);

        bP_prev = dx*dz*(pStar1 - pStar2);
        finiteObj.replaceSection(bP, bP_prev, Iv, jF);

        Matrix::Mat dV_prev = dx*dz/((1/alphaU)*aP);
        finiteObj.replaceSameSection(dV, dV_prev, Iv, jF);

        //Use previous calculation as initial guess in numerical method
        
        Matrix::Mat res2(Nx, vector<Matrix>(Ny-1));
        finiteObj.replaceValue(res2, 0.0);
        for(int k = 0; k < NmaxGSI; k++)
        {
            for(int j = 1; j < Ny; j++)         // Ny = 5
            {
                for(int i = 1; i < Nx+1; i++)     // Nx = 10
                {
                    vOld[i][j].value = (aW[i][j].value*vOld[i-1][j].value + aE[i][j].value*vOld[i+1][j].value + aS[i][j].value*vOld[i][j-1].value +
                                        aN[i][j].value*vOld[i][j+1].value + bP[i][j].value)/(aP[i][j].value/alphaU) + (1-alphaU)*vOld[i][j].value;
                }
            }
            
            for(int j = 1; j < Ny; j++)
            {
                for(int i = 1; i < Nx+1; i++)
                {
                    res2[i-1][j-1].value = (aW[i][j].value*vOld[i-1][j].value + aE[i][j].value*vOld[i+1][j].value + aS[i][j].value*vOld[i][j-1].value + 
                                            aN[i][j].value*vOld[i][j+1].value + bP[i][j].value) - aP[i][j].value*vOld[i][j].value;
                }
            }


            Matrix::Mat aP_ex = extract(aP, oneToNx_2, oneToNy_1)&&extract(vOld, oneToNx_2, oneToNy_1);
            res2_sum = sumMatrixAbs(res2) / sumMatrixAbs(aP_ex);
            
            if(res2_sum < err)
            {
                break;
            }
        }

        finiteObj.replace(vStar, vOld);
        vres[n][0].value = res2_sum;

        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // STEP 2 : Solve pressure correction equation (PCE)
        // Setup coefficients
        // FVM_pcorr

        // Setup coefficients
        finiteObj.nodenumber2(Ip);
        finiteObj.nodenumber2(Jp);


        Matrix::Mat1d Ip2 = -1.0 % Ip;
        Matrix::Mat1d Jp2 = -1.0 % Jp;

        Matrix::Mat aE_prev = rho*dy*dz * extract(dU, Ip, Jp);
        Matrix::Mat aW_prev = rho*dy*dz * extract(dU, Ip2, Jp);
        Matrix::Mat aN_prev = rho*dx*dz * extract(dV, Ip, Jp);
        Matrix::Mat aS_prev = rho*dx*dz * extract(dV, Ip, Jp2);

        finiteObj.replaceSection(aE, aE_prev, Ip, Jp);
        finiteObj.replaceSection(aW, aW_prev, Ip, Jp);
        finiteObj.replaceSection(aN, aN_prev, Ip, Jp);
        finiteObj.replaceSection(aS, aS_prev, Ip, Jp);        

        aP_prev = aE + aW + aN + aS;

        finiteObj.replaceSameSection(aP, aP_prev, Ip, Jp);


        bP_prev = dy*dz*rho*(extract(uStar, Ip2, Jp) - extract(uStar, Ip, Jp)) + dx*dz*rho*(extract(vStar, Ip, Jp2) - extract(vStar, Ip, Jp));
        

        finiteObj.replaceSection(bP, bP_prev, Ip, Jp);

        //Matrix::Mat test = extract(uStar, Ip2, Jp);
        
        // Fix pressure at outlet to zero
        for (int i = 1; i <= Jp.size(); i++)
        {
            aE[Nx][i].value = 0.0;
            aW[Nx][i].value = 0.0;
            aN[Nx][i].value = 0.0;
            aS[Nx][i].value = 0.0;
            aP[Nx][i].value = 1.0;
            bP[Nx][i].value = 0.0;
        }

        finiteObj.replaceValue(pPrime, 0.0);

        // Use numerical method to calculate pressure correction
        // [pPrime, pres(n)] = FVM_GS_ext_mesh(Nx+1, Ny+1, 1, NmaxGSI, err, pPrime);
        Matrix::Mat res3(Nx, vector<Matrix>(Ny));
        finiteObj.replaceValue(res3, 0.0);
        finiteObj.replaceValue(pPrime, 0.0);

        double alpha2 = 1.0;
       

        for(int k = 0; k < NmaxGSI; k++)
        {
            for(int j = 1; j < Ny+1; j++)         // Ny = 5
            {
                for(int i = 1; i < Nx+1; i++)     // Nx = 10
                {
                    pPrime[i][j].value = (aW[i][j].value*pPrime[i-1][j].value + aE[i][j].value*pPrime[i+1][j].value + aS[i][j].value*pPrime[i][j-1].value +
                                        aN[i][j].value*pPrime[i][j+1].value + bP[i][j].value)/(aP[i][j].value/alpha2) + (1-alpha2)*pPrime[i][j].value;                   
                }
                
            }

            
            for(int j = 1; j < Ny+1; j++)
            {
                for(int i = 1; i < Nx+1; i++)
                {
                    res3[i-1][j-1].value = (aW[i][j].value*pPrime[i-1][j].value + aE[i][j].value*pPrime[i+1][j].value + aS[i][j].value*pPrime[i][j-1].value + 
                                            aN[i][j].value*pPrime[i][j+1].value + bP[i][j].value) - aP[i][j].value*pPrime[i][j].value;
                }
            }

            //finiteObj.print2dmat(res3);

            Matrix::Mat aP_ex = extract(aP, oneToNx_2, oneToNy_2)&&extract(pPrime, oneToNx_2, oneToNy_2);
            res3_sum = sumMatrixAbs(res3) / sumMatrixAbs(aP_ex);
            
            if(res3_sum < err)
            {
                break;
            }
        }
       
        pres[n][0].value = res3_sum;



        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // STEP 3 : calculate corrected pressure and velocity

        // p corrections with under -relaxation
        Matrix::Mat p_prev = extract(pStar, Ip, Jp)+alphaP*extract(pPrime, Ip,Jp);
        finiteObj.replaceSection(p, p_prev, Ip, Jp);
        
        // u corrections
        Matrix::Mat uPrime_prev = extract(dU, iF, Ju)&&(extract(pPrime, iF, Ju) - extract(pPrime, iF1, Ju));
        finiteObj.replaceSection(uPrime, uPrime_prev, iF, Ju);

        Matrix::Mat u_prev = extract(uStar, iF, Ju) + extract(uPrime, iF, Ju);
        finiteObj.replaceSection(u, u_prev, iF, Ju);

        // v corrections
        Matrix::Mat vPrime_prev = extract(dV, Iv, jF)&&(extract(pPrime, Iv, jF) - extract(pPrime, Iv, jF1));
        finiteObj.replaceSection(vPrime, vPrime_prev, Iv, jF);

        Matrix::Mat v_prev = extract(vStar, Iv, jF) + extract(vPrime, Iv, jF);
        finiteObj.replaceSection(v, v_prev, Iv, jF);

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // STEP 4 : Check for convergence or divergence
        if(n >= 0 )
        {
            cout << "n = " << n << "\t\t" << "u = " << ures[n][0].value << "\t\t" << "v = " << vres[n][0].value << "\t\t" << "p = " << pres[n][0].value << endl;
            cTest = max(ures[n][0].value, vres[n][0].value);
            if (cTest < err)
                break;
            else if(cTest > 10 || isnan(cTest))
            {
                cout << "Residuals are too high" << endl;
            }
        }

        // Apply right boundary condition (outlet, du/dx = dv/dx = 0)
        for(int j = 0; j < u[0].size(); j++)
        {
            u[Nx][j] = u[Nx-1][j];
        }
        for(int j = 0; j <v[0].size(); j++)
        {
            v[Nx+1][j] = v[Nx][j];
        }

    
    }
    
    finish = clock();

    duration = (double)(finish - start) / CLOCKS_PER_SEC;

    cout << "time : " << duration << "sec" << endl;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // STEP 5 : solve for temperature distribution
    // Setup coefficients, FVM_phi

    // Coefficients for advection-diffusion equation for phi
    // Setup diffusion coefficients, FVM_phi(Nx,Ny,dx,dy,dz,rho,kt/cp,qw/cp,Ip,Jp,u,v,BC_S,BC_N)
    double gamma = kt/cp;
    double dphi = qw/cp;

    double Dx2 = (gamma/dx)*(dy*dz);
    double Dy2 = (gamma/dy)*(dx*dz);


    // 4/21일 수정 필요 ***************************************************************************************

    Matrix::Mat1d Ip2 = -1.0 % Ip;
    Matrix::Mat1d Jp2 = -1.0 % Jp;
    // Setup flow (or advection) coefficients
    Matrix::Mat Fe_prev = rho*dy*dz*extract(u, Ip, Jp);
    Matrix::Mat Fw_prev = rho*dy*dz*extract(u, Ip2, Jp);
    Matrix::Mat Fn_prev = rho*dx*dz*extract(v, Ip, Jp);
    Matrix::Mat Fs_prev = rho*dx*dz*extract(v, Ip, Jp2);

    finiteObj.replaceSection(Fe, Fe_prev, Ip, Jp);
    finiteObj.replaceSection(Fw, Fw_prev, Ip, Jp);
    finiteObj.replaceSection(Fn, Fn_prev, Ip, Jp);
    finiteObj.replaceSection(Fs, Fs_prev, Ip, Jp);

    // Convective flux using hybrid upwinding
    for(int i = Ip[0].value; i <= Ip[Ip.size()-1].value; i++)
    {
        for(int j = Jp[0].value; j <= Jp[Jp.size()-1].value; j++)
        {
            aE[i][j].value =  max(0.0, max(-Fe[i][j].value, Dx2-0.5*Fe[i][j].value));
            aW[i][j].value =  max(0.0, max(Fw[i][j].value, Dx2+0.5*Fw[i][j].value));
            aN[i][j].value =  max(0.0, max(-Fn[i][j].value, Dy2-0.5*Fn[i][j].value));
            aS[i][j].value =  max(0.0, max(Fs[i][j].value, Dy2+0.5*Fs[i][j].value));

            bP[i][j].value = 0.0;
        }
    }

    // cout << "---------------------------------aE---------------------------------------"  << endl;
    // finiteObj.print2dmat(aE);
    // cout << "---------------------------------aW---------------------------------------"  << endl;
    // finiteObj.print2dmat(aW);
    // cout << "---------------------------------aN---------------------------------------"  << endl;
    // finiteObj.print2dmat(aN);
    // cout << "---------------------------------aS---------------------------------------"  << endl;
    // finiteObj.print2dmat(aS);

    if(BC_N == 0)       // Top : wall (no-slip) at Tw or qw
        finiteObj.replaceline(aN, 2*Dy2, Ny, 0);
    else
    {
        finiteObj.replaceline(aN, 0.0, Ny, 0);
        finiteObj.replaceline(bP, dphi*dx*dz, Ny, 0);
    }
    
    if(BC_S == 0)       // Bottom : wall (no-slip) at Tw or symmetry
        finiteObj.replaceline(aS, 2*Dy2, 1, 0);
    else
        finiteObj.replaceline(aS, 0.0, 1, 0);

    finiteObj.replaceline(aE, 0.0, Nx, 1);      // Right boundary (outlet, dT/dx = 0)

    Matrix::Mat DF_prev = Fe - Fw + Fn - Fs;
    finiteObj.replaceSameSection(DF, DF_prev, Ip, Jp);

    Matrix::Mat aP_prev = aE + aW + aN + aS + DF;
    finiteObj.replaceSameSection(aP, aP_prev, Ip, Jp);
    

    // cout << "---------------------------------aE---------------------------------------"  << endl;
    // finiteObj.print2dmat(aE);
    // cout << "---------------------------------aW---------------------------------------"  << endl;
    // finiteObj.print2dmat(aW);
    // cout << "---------------------------------aN---------------------------------------"  << endl;
    // finiteObj.print2dmat(aN);
    // cout << "---------------------------------aS---------------------------------------"  << endl;
    // finiteObj.print2dmat(aS);

    // cout << "Tres1 : " << Tres << endl;

    // Use numerical method to calculate temperature
    Matrix::Mat res4(Nx, vector<Matrix>(Ny));
    finiteObj.replaceValue(res4, 0.0);

    double alpha2 = 1.0;    

    for(int k = 0; k < 1e4; k++)
    {
        for(int j = 1; j < Ny+1; j++)         // Ny = 5
        {
            for(int i = 1; i < Nx+1; i++)     // Nx = 10
            {
                T[i][j].value = (aW[i][j].value*T[i-1][j].value + aE[i][j].value*T[i+1][j].value + aS[i][j].value*T[i][j-1].value +
                                    aN[i][j].value*T[i][j+1].value + bP[i][j].value)/(aP[i][j].value/alpha2) + (1-alpha2)*T[i][j].value;                   
            }
            
        }

        
        for(int j = 1; j < Ny+1; j++)
        {
            for(int i = 1; i < Nx+1; i++)
            {
                res4[i-1][j-1].value = (aW[i][j].value*T[i-1][j].value + aE[i][j].value*T[i+1][j].value + aS[i][j].value*T[i][j-1].value + 
                                        aN[i][j].value*T[i][j+1].value + bP[i][j].value) - aP[i][j].value*T[i][j].value;
            }
        }


        Matrix::Mat aP_ex = extract(aP, oneToNx_2, oneToNy_2)&&extract(T, oneToNx_2, oneToNy_2);
        res4_sum = sumMatrixAbs(res4) / sumMatrixAbs(aP_ex);
        
        if(res4_sum < 10e-8)
        {
            break;
        }
    }
    
    Tres = res4_sum;

    cout << "Tres : " << Tres << endl;

    finiteObj.print2dmat(T);

    return 0;
}
