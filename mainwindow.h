/* c++ */
#ifndef MAINWINDOW_H_INCLUDED
#define MAINWINDOW_H_INCLUDED

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <string>

#include "percol2d.h"

//---------------------------------------------------------------------
// Define main window 

struct X_of_T
{
    double T;
    double G;
};

struct MainWindow
{
    MainWindow();
    ~MainWindow(){ clear(); }
    void setModel();

    void computeModel();
    void clear();
    void selectSigma(int);
    void computeAreaE();
    double AreaE(double E);
    void computeRT(std::vector<X_of_T> & data);
    void computeReffT(std::vector<X_of_T> & data);
    void computeRU(std::vector<X_of_T> & data);
    void computeReffU(std::vector<X_of_T> & data);
    void computeEF_TU();
    void computeEFT();
    void computeEFU();

    void randomizeSigma_0();
    void randomizeSigma_1();
    void randomizeSigma_2();

    inline double singleSigma(double r, double rEx);
    inline double singleSigmaT0(double E, double V, double r, double EFc);
    inline double sedlo(double E, double Ey, double Ex, double V);
    inline double Vbarrier(double r);
    inline double computeSum(int NE, double dE, double Vd, double EFTU);
    inline double Vdot(void);
    double effective_medium(double y);
    inline double average(double y);

    void randRcr();
    Percol2D *model;
    int typeResistor;

    double T, U, Tmin, Tmax, dT, Umin, Umax, dU, Ex, Ey, EF,EFT;
    int cols, rows, seed; 
    double sigmaU, sigmaMin, r_c;
    double kappa, CUTOFF_SIGMA;
    std::vector<double> EFTarray;
    std::vector<double> EFUarray;
    std::vector<double> AreaEf;
//    std::vector<double> CondDist;
//    std::vector<double> NumJDist;
//    std::vector<int> index_for_sorted_CondDist() const;
//    std::vector<int> index_for_sorted_NumJDist() const;
    double gTun, gOv;
    double dispCond; // Display of condition number of last computation
    double dispCapac; // Display of capacity value of last computation
    double dispConduct; // Display of conductivity value of last computation
    double dispFerr; // Display of ferr
    double dispBerr; // Display of berr
    double dispDeltaI; // Display of deltaI
    int typeCond;
};

#endif /*MAINWINDOW_H_INCLUDED*/
