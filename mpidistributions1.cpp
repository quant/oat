#define _CRT_SECURE_NO_WARNINGS 1

#include <mpi.h>
#include "mkl.h"
#include "percol2d.h"
#include "mainwindow.h"
#include <algorithm>
#include <vector>
#include <valarray>
#ifndef _WIN32
#include <unistd.h>
#endif

#define MYARRAY  std::valarray
#define MYVECTOR  std::vector
//#define IAMHERE printf("I (%i) am at %s:%i\n",myrank,__FILE__,__LINE__)
#define IAMHERE (void)0

double global_Tmin = 0.1, global_Tmax = 5.3, global_dT = 0.325;//1.3;//2.6;//0.2;
double global_Umin = 165., global_Umax = 240.0, global_dU = 5;
double global_U = 230.;//235.;//190;
double global_T = 0.13;//235.;//190;

int global_rows=-1;//250;
int global_cols=-1;//253;
int global_seed=1;
int npV=201;
int npI=201;
int npJ=201;
int nRollsPerNode=-1;
double global_kappa=30;
double global_CUTOFF_SIGMA=1.e-13;
//0 -- kappa, 1 -- 0 or 1, 2 -- experiment
int global_typeResistor=-1;//2;

static void parse_args(int argc, char *argv[])
{
    for (int i=0; i<argc; ++i)
    {
        if (1==sscanf(argv[i],"global_rows=%i",&global_rows)) continue;
        if (1==sscanf(argv[i],"global_cols=%i",&global_cols)) continue;
        if (1==sscanf(argv[i],"nRollsPerNode=%i",&nRollsPerNode)) continue;
        if (1==sscanf(argv[i],"global_typeResistor=%i",&global_typeResistor)) continue;
    }
    if (global_rows==-1) goto err;
    else printf("global_rows=%i\n",global_rows);
    if (global_cols==-1) goto err;
    else printf("global_cols=%i\n",global_cols);
    if (nRollsPerNode==-1) goto err;
    else printf("nRollsPerNode=%i\n",nRollsPerNode);
    if (global_typeResistor==-1) goto err;
    else printf("global_typeResistor=%i\n",global_typeResistor);
    return;
err:
    printf("Error in parse_args\n");
    exit(1);
}

static void print_mkl_version(void);

static int main1(int argc, char **argv);
int main(int argc, char **argv)
{
    int res;
    printf("Starting\n"); fflush(0);
    res = MPI_Init( &argc, &argv );
    printf("MPI_Init returned %i\n",res); fflush(0);
    res = main1( argc, argv );
    printf("main1 returned %i\n",res); fflush(0);
    res = MPI_Finalize();
    printf("MPI_Finalize returned %i\n",res); fflush(0);
    return res;
}

static int main1(int argc, char **argv)
{
    MYARRAY<double> CondDist;
    MYARRAY<double> GeffDist;

    parse_args(argc,argv);

    MainWindow mainWindow;

    int mysize, myrank;

    MPI_Comm_size(MPI_COMM_WORLD,&mysize);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    printf("I am %i of %i\n",myrank,mysize);
    mainWindow.U=global_U;
    mainWindow.T=global_T;
    mainWindow.rows=global_rows;
    mainWindow.cols=global_cols;
    mainWindow.CUTOFF_SIGMA=global_CUTOFF_SIGMA;
    mainWindow.kappa=global_kappa;
    mainWindow.typeResistor=global_typeResistor;
    CondDist.resize(nRollsPerNode,0);
    GeffDist.resize(nRollsPerNode,0);


    // Each process will do nRollsPerNode realizations, so 'mysize'
    // processors will do nR realizations

    int iseed = global_seed + nRollsPerNode * myrank * 10;
    // 10 means we assume giving 10 seeds per roll is enough
    double t_final=60*20500;//calculation time in seconds
    double dVmin=log(1.e-20);
    double dVmax=log(2.);
    double dV=(dVmax-dVmin)/(npV+1);
    MYVECTOR<double> pdV( npV+2,0);
    double imin=log(1.e-25);
    double imax=-1.;
    double dI=(imax-imin)/(npI+1);
    MYVECTOR<double> pI( npI+3 );
    double Jmin=log(1.e-25);
    double Jmax=-1.;
    double dJ=(Jmax-Jmin)/(npJ+1);
    MYVECTOR<double> pJ( npJ+3,0 );
    double cond_min=1e-12;
    //double logGmin=log(cond_min);
    double logmin=log(mainWindow.CUTOFF_SIGMA);
    int e=0;
    double sumI=0.; 
    double y_old=0;
    double sigma_mid=0; 
    // FIRST IT IS NECESSARY TO CALCULATE THE FERMI ENERGY
    if(global_typeResistor==2) mainWindow.computeEF_TU();
    if (myrank==0)
    {
        printf( "rows=%i cols=%i seed=%i nRollsPerNode=%i\n", 
            mainWindow.rows, mainWindow.cols, iseed, nRollsPerNode);
	if(global_typeResistor==2)
	    printf( "U=%lg T=%lg EFT=%lg\n", 
                    mainWindow.U, mainWindow.T, mainWindow.EFT);

	if(global_typeResistor==0) printf("kappa=%lg\n", global_kappa);
        printf( "CUTOFF_SIGMA=%lg\n",global_CUTOFF_SIGMA);
   
        print_mkl_version();
    }
    double t_begin = MPI_Wtime();
    double t_current,t_av;
IAMHERE;
    int j;
    double yeff=0.;
    for (j = 0; j < nRollsPerNode; ++j)
    {
//printf("j=%i\n",j);
        int ni;
        double yy;
        for (;;) // in fact we have only 10 seeds to try in the worst case
        {
            mainWindow.seed = iseed;
IAMHERE;
            mainWindow.computeModel();
IAMHERE;
            iseed++;
            ni = mainWindow.model->nI();
IAMHERE;
            yy = mainWindow.model->conductivity;
            //yn = mainWindow.model->capacity;

            // If conductivity is ok, proceed with computation (break
            // from this loop), otherwise try another realization
            if (yy > cond_min) break;
        }

        CondDist[j] = yy;
 //       double yeff;
        if(global_typeResistor!=1) 
        {
           yeff = mainWindow.effective_medium(y_old);
 //       printf( "myrank=%i j=%i G=%lg Geff=%lg\n",myrank,j,yy,yeff);
           y_old=yeff;
           GeffDist[j] = yeff;//yn;
        }
IAMHERE;
        if(global_typeResistor!=1) sigma_mid=sigma_mid+yeff;
        for (int i = 0; i < ni; ++i)
        {
            double Vi = log(fabs(mainWindow.model->difV[i]));
            double Ji = log(fabs(mainWindow.model->IdifV[i]));
            double Ii = log(fabs(mainWindow.model->I[i]));
            if(mainWindow.model->Sigma[i]<mainWindow.sigmaU)
            {
                sumI=sumI+1; 
                //current distribution          
                if(Ii<imin) 
                {
                  pI[0]++;
                    //              pI[myrank,0]=pI[myrank,0]+1;
                }
                else
                {
                    e=1+int((Ii-imin)/dI);
                    if(e>(npI+2)) pI[npI+2]=pI[npI+2]+1;
                    if(e==(npI+2)) pI[e-1]=pI[e-1]+1;
                    if(e<(npI+2)) pI[e]=pI[e]+1;
                }
                //Joule heat distribution          
                if(Ji<Jmin) 
                {
                    pJ[0]=pJ[0]+1;
                }
                else
                {
                    e=1+int((Ji-Jmin)/dJ);
                    if(e>(npJ+2)) pJ[npJ+2]=pJ[npJ+2]+1;
                    if(e==(npJ+2)) pJ[e-1]=pJ[e-1]+1;
                    if(e<(npJ+2)) pJ[e]=pJ[e]+1;
                }
                //voltage dufference distribution          
                if(Vi<dVmin) 
                {
                    pdV[0]=pdV[0]+1;
                }
                else
                {
                    e=1+int((Vi-dVmin)/dV);
#if 0
                    if(e>=(npV+2)) pdV[e-1]=pdV[e-1]+1; //XXX:e-1 out of range
                    if(e<(npV+2)) pdV[e]=pdV[e]+1;
#else
                    if (e >= pdV.size()) e = pdV.size()-1;
                    pdV[e] += 1;
#endif
                }
            }
        }
    t_current = MPI_Wtime()-t_begin;
    t_av=t_current/(j+1);
    if(t_current>(t_final-2*t_av)) break;

    }
IAMHERE;
    int Jmyrank=j+1;
    if(j==nRollsPerNode) Jmyrank=j;   
        printf( "myrank=%i Jmyrank=%i t_av=%lg t_culc=%lg\n",myrank,j,t_av,t_current);
    //--------------------------------------
    //---------conductance distribution pG
    //---------the 0.1 heat fraction distribution pNumJ
#if 0
    if (0) {
        int k;
        FILE *f = fopen("cd1.txt","wt");
        for (k=0; k<CondDist.size(); ++k)
            fprintf(f,"CondDist[%i]=%lg\n",k,CondDist[k]);
        for (k=0; k<idxS.size(); ++k)
            fprintf(f,"sorted.CondDist[%i]=%lg\n",idxS[k],CondDist[idxS[k]]);
        fclose(f);
    }
#endif
    //-----------------------------
//IAMHERE;
    double t_mid = MPI_Wtime();
//IAMHERE;
    //-----------------------------
    // sending and receiving data
    //CondDist[0:nRollsPerNode-1], GeffDist[0:nRollsPerNode-1], pdV[0:npV+1]  pI[0:npI+2] pJ[0:npJ+2] sigma_mid
    // Reducing number of calculated conductance values: sum of Jmyrank
     int Jtotal;
    {
        int res;
        MPI_Reduce(&Jmyrank, &res, 1,
            MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        Jtotal = res;
    }
IAMHERE;
    MYVECTOR<double> totalG, totalGeff;
    {
        int NITEMS = CondDist.size();
        if(myrank==0) totalG.resize(mysize*NITEMS);

        MPI_Gather(&CondDist[0], NITEMS, MPI_DOUBLE,
            &totalG[0], NITEMS, MPI_DOUBLE,
            0, MPI_COMM_WORLD);
        std::sort(totalG.begin(),totalG.end());
//      int k;
 
//     if(myrank==0){
//      printf("NITEMS=%i totalG.size()=%i\n",NITEMS,totalG.size());
//     for (k=0; k<totalG.size(); ++k)
//           printf("totalG[%i]=%lg\n",k,totalG[k]);
//}
 
       if(global_typeResistor!=1)
         {
        if(myrank==0) totalGeff.resize(mysize*NITEMS);
        MPI_Gather(&GeffDist[0], NITEMS, MPI_DOUBLE,
            &totalGeff[0], NITEMS, MPI_DOUBLE,
            0, MPI_COMM_WORLD);
        std::sort(totalGeff.begin(),totalGeff.end());
         }
    }
//IAMHERE;
    // Reducing pdV
    {
        MYVECTOR<double> res( pdV.size() );
        MPI_Reduce(&pdV[0], &res[0], pdV.size(),
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        pdV = res;
    }
//IAMHERE;

    // Reducing pI
    {
        MYVECTOR<double> res( pI.size() );
        MPI_Reduce(&pI[0], &res[0], pI.size(),
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        pI = res;
    }
//IAMHERE;

    // Reducing pJ
    {
        MYVECTOR<double> res( pJ.size() );
        MPI_Reduce(&pJ[0], &res[0], pJ.size(),
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        pJ = res;
    }
//IAMHERE;
/*
    // Reducing pS1
    {
        MYVECTOR<double> res( pS1.size() );
        int s = MPI_Reduce(&pS1[0], &res[0], pS1.size(),
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        pS1 = res;
    }

    // Reducing pNumJ
    {
        MYVECTOR<double> res( pNumJ.size() );
        int s = MPI_Reduce(&pNumJ[0], &res[0], pNumJ.size(),
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        pNumJ = res;
    }
*/
    // Reducing sigma_mid
    {
        double res;
        MPI_Reduce(&sigma_mid, &res, 1,
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        sigma_mid = res;
    }
//IAMHERE;

    //--------------------------------------
    // final processing
    if (myrank != 0) 
    {
    	MPI_Finalize();
	exit(0);
    }
    int nG=51;
    int nGeff=51;
    MYVECTOR<double> pG( nG+1,0 );// from 0-box to Ns-box
    MYVECTOR<double> pGeff( nGeff+1,0 );// from 0-box to Ns-box
    double Gmin, Geffmin;
    for (int i = 0; i < totalG.size(); ++i)
    {   double y=totalG[i];
        if(y>0) { 
           Gmin=y;
           break;}
    }
    if(global_typeResistor!=1)
    {
    for (int i = 0; i < totalGeff.size(); ++i)
    {   double y=totalGeff[i];
        if(y>0) { 
           Geffmin=y;
           break;}
    }
    }  
//    double Gmin = totalG[0];
    double Gmax = totalG[totalG.size()-1];
//    Gmin = ( Gmin>0 ? log(Gmin): logGmin);
    Gmin = log(Gmin);
    Gmax = log(Gmax);
    double dG = (Gmax-Gmin)/(nG+1);
    double dGeff=0;
//    double Geffmin=totalGeff[0];
    if(global_typeResistor!=1)
    {
    double Geffmax=totalGeff[totalGeff.size()-1];
//    Geffmin = (Geffmin>0 ? log(Geffmin) : logGmin);
    Geffmin = log(Geffmin);
    Geffmax = log(Geffmax);
    dGeff=(Geffmax-Geffmin)/(nGeff+1);
    }
//IAMHERE;
    int jG=0;
    int jGeff=0;
    for (int i = 0; i < totalG.size(); ++i)
    { 
        double y=totalG[i];
        if(y>0) 
         {y=log(y);
          e=int((y-Gmin)/dG);
          if(e==(nG+1)) pG[e-1]=pG[e-1]+1;
          else pG[e]=pG[e]+1;
          jG=jG+1;}
        if(global_typeResistor!=1)
        {
        double yn=totalGeff[i];
        if(yn>0) 
         {yn=log(yn);
          int e1=int((yn-Geffmin)/dGeff);
          if(e1==(nGeff+1)) pGeff[e1-1]=pGeff[e1-1]+1;
          else pGeff[e1]=pGeff[e1]+1;
          jGeff=jGeff+1;}
        }
    }

//IAMHERE;
    sigma_mid /= Jtotal;
//    sigma_mid /= nRollsPerNode*mysize;
    printf("Effective sigma_mid=%lg Jtotal=%i jG=%i jGeff=%i\n",sigma_mid,Jtotal,jG,jGeff);

    double p;
    double x0=double(mainWindow.cols-3);
    x0=log(x0);
    printf("log(sigma_mid)/x0=%lg\n",log(sigma_mid)/x0);

    double sumV = Jtotal * dV;
//    double sumV = nRollsPerNode*mysize * dV;

    p = (pdV[0] > 0 ? log(pdV[0]/sumV)/x0 : logmin);

    printf( "Distribution dV -----------------------\n");
    printf( "%20s %20s\n", "x/x0", "log(pdV/sumV)/x0");
    printf( "%20.7lg %20.7lg\n", dVmin/x0, p);
    for (int i = 1; i < npV+2; ++i)
    {   
        double x = dVmin + i*dV - 0.5*dV;
        p = (pdV[i] > 0 ? log(pdV[i]/sumV)/x0 : logmin);
        printf( "%20.7lg %20.7lg\n", x/x0, p);
    }
    //current
    double sumCur = Jtotal*dI;
//    double sumCur = nRollsPerNode*mysize*dI;
    p = (pI[0] > 0 ? log(pI[0]/sumCur)/x0 : logmin);

    printf( "Distribution dI -----------------------\n");
    printf( "%20s %20s\n", "x/x0", "log(pI/sumCur)/x0");
    printf( "%20.7lg %20.7lg\n", imin/x0, p);
    for (int i = 1; i < npI+2; ++i)
    {   
        double x = imin + i*dI - 0.5*dI;
        p = (pI[i] > 0 ? log(pI[i]/sumCur)/x0 : logmin);
        printf( "%20.7lg %20.7lg\n", x/x0, p);
    }
    p = (pI[npI+2] > 0 ? log(pI[npI+2]/sumCur)/x0 : logmin);
    printf( "%20.7lg %20.7lg\n", imax/x0, p);

    //---------Joule Heat
    double sumJ = Jtotal*dJ;
//    double sumJ = nRollsPerNode*mysize*dJ;
    p = (pJ[0] > 0 ? log(pJ[0]/sumJ)/x0 : logmin);
    printf( "Distribution dJ -----------------------\n");
    printf( "%20s %20s\n", "x/x0", "log(pJ/sumJ)/x0");
    printf( "%20.7lg %20.7lg\n", Jmin/x0, p);
    for (int i = 1; i < npJ+2; ++i)
    {   
        double x=Jmin+i*dJ-0.5*dJ;
        p = (pJ[i] > 0 ? log(pJ[i]/sumJ)/x0 : logmin);
        printf( "%20.7lg %20.7lg\n", x/x0, p);
    }
    p = (pJ[npJ+2] > 0 ? log(pJ[npJ+2]/sumJ)/x0 : logmin);
    printf( "%20.7lg %20.7lg\n", Jmax/x0, p);
    // conductance distribution
    double sumG=Jtotal*dG;  
//    double sumG=nRollsPerNode*mysize*dG;  
    printf( "Distribution pG ----------------------\n");
    printf( "%20s %20s\n", "x/x0", "log(pG/sumG)/x0");
    for (int i = 0; i < (nG+1); ++i)
    {   
        double x=Gmin+i*dG+0.5*dG;
        p = pG[i]*exp(x);
        p = (p > 0 ? log(p/sumG)/x0 : logmin);
        printf( "%20.7lg %20.7lg\n", x/x0, p);
    }
    // conductance distribution
    if(global_typeResistor!=1)
    {
    double sumGeff=Jtotal*dGeff;  
//    double sumGeff=nRollsPerNode*mysize*dGeff;  
    printf( "Distribution of effective conductivity-\n");
    printf( "%20s %20s\n", "x/x0", "log(pGeff/sumGeff)/x0");
    for (int i = 0; i < (nGeff+1); ++i)
    {   
        double x=Geffmin+(i+0.5)*dGeff;
        p = pGeff[i]*exp(x);
        p = (p > 0 ? log(p/sumGeff)/x0 : logmin);
        printf( "%20.7lg %20.7lg\n", x/x0, p);
    }
    }
    //        double t_mid = MPI_Wtime();
    double t_end = MPI_Wtime();


    printf( "t_mid-t_begin=%lg\n",t_mid-t_begin);
    printf( "time_exec=%lg\n",t_end-t_begin);

    return 0;
}

static void print_mkl_version(void)
{
    MKLVersion v;
    MKLGetVersion(&v);
    printf("Using Intel MKL %i.%i.%i %s %s for %s/%s\n",
        v.MajorVersion, v.MinorVersion, v.UpdateVersion,
        v.ProductStatus, v.Build, v.Processor, v.Platform);
}
