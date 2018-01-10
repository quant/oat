#include "mainwindow.h"
#include "mkl.h"
#include <math.h>
#include <algorithm>
#include <vector>
#define MYVECTOR std::vector
#include "PercolRectSides.h"

//const int NCUT = 1.e4;
//const int nJ=5000;
const double E0=560.;//meV
const double Vg0=350.;
const double Delta_r=22.;
const double Cg0=-0.09;//12;//-0.05;//-0.04;
//double sigma_m=0.1;
const double G_ser=1.;//0.9;//1.;//1.;//1.;//3.;
double a_barrier=180.;
double EF0=35;//42.;
void MainWindow::clear(void)
{
}

void MainWindow::setModel()
{
    if (model) delete model;
    model = new PercolRectSides(this->rows,this->cols);
}

MainWindow::MainWindow()
: typeCond(3),sigmaU(100.0),
T(0.1), Tmin(0.1),Tmax(5.301), dT(0.1),
U(460.), Umin(460.), Umax(600.), dU(10.),
r_c(0.0),Ex(22.),Ey(6.), EF(22.),EFT(22.),kappa(10.),
//r_c(0.0),Ex(22.),Ey(6.), EF(22.),EFT(22.),kappa(10.),
rows(90), cols(153), seed(1), model(0)
{
    this->setModel();
    this->typeResistor = 2;//0 sigma_j=exp(-kappa*r_j), 1-- sigma_j is 0 or 1, 2 -- two diiferent exponent
    this->selectSigma( this->typeResistor );
}

void MainWindow::computeModel()
{
    int r_type = this->typeResistor;
    this->selectSigma(r_type);
    model->compute();
}

void MainWindow::computeRT(MYVECTOR<X_of_T> & data)
{
    int G_type = typeCond;
    if (G_type != 4)
    {
        computeEFT();
    }
    int NT=1 +int((this->Tmax-this->Tmin)/this->dT);
    for(int j=0; j<NT; j++)
    {
        double x=this->Tmax-this->dT*j;
        this->T=x;
        if(G_type==4) this->EFT=EF0;
        else
        {
            double EFTU=this->EFTarray[j];
            this->EFT=EFTU;
        }
        this->computeModel();
        double y=model->conductivity;
        y=(this->cols-3)*y/this->rows;
        y=6.28*y;
        printf("T=%lg G=%lg\n",x,y);
        X_of_T  xy;
        xy.T = x;
        xy.G = y;
        data.push_back(xy);
    }
}
void MainWindow::computeEF_TU()
{
    this->EFT=36;//30;//27;//25.965;//27.;
}
/*void MainWindow::computeEF_TU()
{
    double E,EFT0,EFT1,EFT2 ;
    double dE=0.1;
    double sum, sum1, Area, sum10, sum11, sum12;
    double Ucur=this->U;//!!!!!!!!!!!
    this->U=Vg0;
    double Vd0=Vdot();
    this->U=Ucur;
    double aa=(1000.*sqrt(1.25)-350.)/Delta_r;
    aa=aa*aa;
    aa=4*this->U/(1+aa);
//    aa=0;
    if(this->T==0)
    {
        EFT1=aa+EF0+Vdot()-Vd0+Cg0*(Ucur-Vg0);//!!!!!!!!!!
    }
    else
    {
        //T!=0
        double Vd=Vdot();
        this->EF=aa+EF0+Vd-Vd0+Cg0*(Ucur-Vg0);
        int NE=int((this->EF+40*this->T-Vd)/dE);
        this->AreaEf.resize(NE,0.0);
        sum=0;
        EFT1=this->EF-1;//!!!!!!!!!!!!!!!!!!!!!!!!
        for (int i=0; i< NE; ++i)
        {
            E=dE*(i+1)+Vd;
            Area=AreaE(E)/10000;
            this->AreaEf[i]=Area;
            if(E<=this->EF) sum=sum+Area;
        }
        //        this->density=sum;
        if(sum!=0)
        {
            EFT0=EFT1;//this->EF-1.;//!!!!!!!!!!!!!!!!!
            sum1=computeSum(NE, dE, Vd, EFT0);
            while(sum1>sum)
            {
                EFT0=EFT0-1;
                sum1=computeSum(NE, dE, Vd, EFT0);
            }
            sum10=sum1;
            EFT1=EFT0+1;
            sum11=computeSum(NE, dE, Vd, EFT1);
            //        int j=0;
            while(fabs(sum11-sum)>0.001*sum)
            {
                EFT2=EFT1-(sum11-sum)*(EFT1-EFT0)/(sum11-sum10);
                sum12=computeSum(NE, dE, Vd, EFT2);
                //            j++;
                if(sum12>sum&&sum11<sum||sum12<sum&& sum11>sum)
                {
                    sum10=sum11;
                    EFT0=EFT1;
                }
                sum11=sum12;
                EFT1=EFT2;
            }
        }
    }
    this->EFT=EFT1;
}
*/
void MainWindow::computeEFT()
{
    double E,EFT0,EFT1,EFT2 ;
    int NT=1+int( (this->Tmax-this->Tmin)/this->dT );
    this->EFTarray.resize(NT,0.0);
    //    printf("EFTarray has %i points",NT);
    double dE=0.1;
    double Ucur=this->U;
    this->U=Vg0;
    double Vd0=Vdot();
    this->U=Ucur;
    double Vd=Vdot();
    double aa1=(1000.*sqrt(1.25)-350.)/Delta_r;
    aa1=aa1*aa1;
    aa1=4./(1+aa1);
//    aa1 = 0;
    double aa=aa1*this->U;
    this->EF=aa+EF0+Vd-Vd0+Cg0*(this->U-Vg0);
    int NE = 1+int( (this->EF+40*this->Tmax-Vdot())/dE );
    if(NE<0) return;
    printf("AreaEf has %i points",NE);
    this->AreaEf.resize(NE,0.0);
    double sum, sum1, Area, sum10, sum11, sum12;
    sum=0;
    for (int i=0; i< NE; ++i)
    {
        E=dE*(i+1)+Vd;
        Area=AreaE(E)/10000;
        this->AreaEf[i]=Area;
        if(E<=this->EF) sum=sum+Area;
    }
    if(sum==0) return;
    else
    {
        EFT1=this->EF-1.;
        for(int j=0; j<NT; j++)
        {
            this->T=this->Tmax-this->dT*j;
            EFT0=EFT1;
            sum1=computeSum(NE, dE, Vd, EFT0);
            while(sum1>sum)
            {
                EFT0=EFT0-1;
                sum1=computeSum(NE, dE, Vd, EFT0);
            }
            sum10=sum1;
            EFT1=EFT0+1;
            sum11=computeSum(NE, dE, Vd, EFT1);
            while(fabs(sum11-sum)>0.001*sum)
            {
                EFT2=EFT1-(sum11-sum)*(EFT1-EFT0)/(sum11-sum10);
                sum12=computeSum(NE, dE, Vd, EFT2);
                if(sum12>sum&&sum11<sum||sum12<sum&& sum11>sum)
                {
                    sum10=sum11;
                    EFT0=EFT1;
                }
                sum11=sum12;
                EFT1=EFT2;
            }
            this->EFTarray[j]=EFT1;
        }
    }
}

double MainWindow::computeSum(int NE, double dE, double Vd, double EFT)
{
    double sum=0;
    for (int i=0; i< NE; ++i)
    {
        double E=dE*(i+1)+Vd;
        sum=sum+this->AreaEf[i]/(1+exp((E-EFT)/this->T));
    }
    return sum;
}

void MainWindow::computeReffT(MYVECTOR<X_of_T> & data)
{
    int r_type = this->typeResistor;
    int G_type = typeCond;
    if(G_type!=4)
    {
        computeEFT();
    }
    int NT=1 +int((this->Tmax-this->Tmin)/this->dT);
    double y_old=0;
    for(int j=0; j<NT; j++)
    {
        double x=this->Tmax-this->dT*j;
        this->T=x;
        if(G_type==4) this->EFT=EF0;
        else
        {
            double EFTU=this->EFTarray[j];
            this->EFT=EFTU;
        }
        this->selectSigma(r_type);
        double y=effective_medium(y_old);
        y_old=y;
        double y1=6.28*y;
        X_of_T  xy;
        xy.T = x;
        xy.G = y1;
        data.push_back(xy);
    }
}

void MainWindow::computeRU(MYVECTOR<X_of_T> & data)
{
    int G_type = typeCond;
    if(G_type!=4)  computeEFU();
    int NU=1 +int((this->Umax-this->Umin)/this->dU);
    for(int j=0; j<NU; j++)
    {
        double x=this->Umin+this->dU*j;
        this->U=x;
        if(G_type==4) this->EFT=EF0;
        else
        {
            double EFTU=this->EFUarray[j];
            this->EFT=EFTU;
            if(EFTU<=0) break;

        }
        this->computeModel();
        double y=model->conductivity;
        y=(this->cols-3)*y/this->rows;
        y=6.28*y;
        printf("T=%lg G=%lg\n",x,y);
        X_of_T  xy;
        xy.T = x;
        xy.G = y;
        data.push_back(xy);
    }
}

void MainWindow::computeReffU(MYVECTOR<X_of_T> & data)
{
    int r_type = this->typeResistor;
    int G_type = typeCond;
    if(G_type!=4) computeEFU();
    int NU=1 +int((this->Umax-this->Umin)/this->dU);
    double y_old=0;
    for(int j=0; j<NU; j++)
    {
        double x=this->Umin+this->dU*j;
        this->U=x;
        if(G_type==4) this->EFT=EF0;
        else
        {
            double EFTU=this->EFUarray[j];
            this->EFT=EFTU;
            if(EFTU<=0) break;

        }
        this->selectSigma(r_type);
        double y=effective_medium(y_old);
        y_old=y;
        double y1=6.28*y;
        printf("T=%lg G=%lg\n",x,y);
        X_of_T  xy;
        xy.T = x;
        xy.G = y1;
        data.push_back(xy);
    }
}

double MainWindow::average(double y)
{
    double sum=0;
    for (int i = 0; i < model->nI(); ++i)
    {
        double sigma_ik=fabs(model->Sigma[i]);
        if (sigma_ik!=this->sigmaU)
        {
            sum = sum+(sigma_ik-y)/(sigma_ik+fabs(y));
        }
    }
    return sum;
}

double MainWindow::effective_medium(double y_old)
{
    double y0=y_old;
    double sum10=0;
    double sum11=1;
    double Totsum=0;
    double y1=0;
    for (int i = 0; i < model->nI(); ++i)
    {
        double sigma_ik=fabs(model->Sigma[i]);
        if (sigma_ik!=this->sigmaU)
        {
            Totsum=Totsum+1;
            sum10 = sum10+(sigma_ik-y0)/(sigma_ik+fabs(y0));
        }
    }
    if(y0==0)
    {
        sum10=1;
        double dy=0.1;
        double sumold;
        while(sum11>0)
        {
            y1=y1+dy;
            sumold=sum11;
            sum11=average(y1)/Totsum;
        }
        if(sumold!=1)
        {
            sum10=sumold;
            y0=y1-dy;
        }
    }
    else
    {
        sum10=average(y0)/Totsum;
        if(sum10 < 0) y1=0.5*y0;
        else y1=1.5*y0;
        sum11=average(y1)/Totsum;
    }
    double y2=0;
    double y11=y0;
    int j=0;
    while(fabs(y1-y11)>0.001*y1)
//    while(fabs(y1-y11)>0.0001*y1)
    {
        y2=y1-sum11*(y1-y0)/(sum11-sum10);
        double sum12=average(y2)/Totsum;
        j++;
        if(sum12>0&&sum11<0||sum12<0&& sum11>0)
        {
            sum10=sum11;
            y0=y1;
        }
        sum11=sum12;
        y11=y1;
        y1=y2;
    }
    return y1;
}

double MainWindow::AreaE(double E)
{
    double Uxy,x1,x2, y1, y2, r1, r2, r3, r4;
    double Area=0;
    for (double x =0.5; x <= 499.5; x += 1)
    {
        for (double y =0.5; y <= 499.5; y += 1)
        {
            if((x-500)*(x-500)+(y-500)*(y-500)>=122500)
            {
                x1=500+x;
                x2=500-x;
                y1=500+y;
                y2=500-y;
                r1=sqrt(x1*x1+y2*y2)-350;
                r2=sqrt(x2*x2+y2*y2)-350;
                r3=sqrt(x1*x1+y1*y1)-350;
                r4=sqrt(x2*x2+y1*y1)-350;
                if(r1>0&&r2>0&&r3>0&&r4>0)
                {
                    r1=r1/Delta_r;
                    r2=r2/Delta_r;
                    r3=r3/Delta_r;
                    r4=r4/Delta_r;
                    Uxy=this->U*(1/(1+r1*r1)+1/(1+r2*r2)+1/(1+r3*r3)+1/(1+r4*r4));
                    if(Uxy<=E) Area=Area+1;
                }
            }
        }
    }
    return Area;
}

void MainWindow::computeEFU()
{
    double E,EFT0,EFT1,EFT2 ;
    double dE=0.1;
    double sum, sum1, Area, sum10, sum11, sum12;
    int NU=1+int( (this->Umax-this->Umin)/this->dU );
    this->EFUarray.resize(NU,0.0);
    this->U=Vg0;
    double Vd0=Vdot();
    double aa1=(1000.*sqrt(1.25)-350.)/Delta_r;
    aa1=aa1*aa1;
    aa1=4./(1+aa1);
    if(this->T==0) {
        for(int j=0; j<NU; j++)
        {
            double x=this->Umin+this->dU*j;
            this->U=x;
            double aa=aa1*this->U;
//            aa = 0;
            EFT1=aa+EF0+Vdot()-Vd0+Cg0*(this->U-Vg0);
            this->EFUarray[j]=EFT1;
            this->EFT=EFT1;
        }
    }
    for(int j=0; j<NU; j++)
    {
        double x=this->Umin+this->dU*j;
        this->U=x;
        double Vd=Vdot();
        double aa=aa1*this->U;
//        aa = 0;
        this->EF=aa+EF0+Vd-Vd0+Cg0*(this->U-Vg0);
        int NE = int( (this->EF+40*this->T-Vd)/dE );
        this->AreaEf.resize(NE,0.0);
        sum=0;
        EFT1=this->EF-1;
        for (int i=0; i< NE; ++i)
        {
            E=dE*(i+1)+Vd;
            Area=AreaE(E)/10000;
            this->AreaEf[i]=Area;
            if(E<=this->EF) sum=sum+Area;
        }
        if(sum>0)
        {
            EFT0=EFT1;
            sum1=computeSum(NE, dE, Vd, EFT0);
            while(sum1>sum)
            {
                EFT0=EFT0-1;
                sum1=computeSum(NE, dE, Vd, EFT0);
            }
            sum10=sum1;
            EFT1=EFT0+1;
            sum11=computeSum(NE, dE, Vd, EFT1);
            while(fabs(sum11-sum)>0.001*sum)
            {
                EFT2=EFT1-(sum11-sum)*(EFT1-EFT0)/(sum11-sum10);
                sum12=computeSum(NE, dE, Vd, EFT2);
                if(sum12>sum&&sum11<sum||sum12<sum&& sum11>sum)
                {
                    sum10=sum11;
                    EFT0=EFT1;
                }
                sum11=sum12;
                EFT1=EFT2;
            }
            this->EFUarray[j]=EFT1;
            this->EFT=EFT1;
        }
    }
}

double MainWindow::sedlo(double E, double Ey, double Ex, double V)
{
    double  alpha,G0,g,exp0,EE, Ep,Uc;
    Uc=Vdot();
    double a1=a_barrier;
    double deltaV = 10;//11;
//    double deltaV = 11;
    double Va=V-0*deltaV;//meV
    V=(V+(this->U-364.3)/120.*this->T+0.75*deltaV)*(1+sqrt(0.1/this->T)); //Ef=36; U=453 a
    double a0=10.;
//    double a0=2/Ex*sqrt(E0*(V-Va));
    Ex=2./a0*sqrt(E0*(V-Va));
    double a2=a0/a1;
    double U00=V-Va;
    double U01=Va/(1-a2*a2);
    this->gTun=0;
    this->gOv=0;
    G0=0;
    EE=E-0.5*Ey-V;
    Ep=E-0.5*Ey;
    while(Ep>Uc)
    {
        if(Ep>Va)
        {
            alpha=-6.2832*EE/Ex;
            exp0=exp(alpha);
            g=1./(1+exp0);
        }
        else
        {
            double b0=sqrt(U00/(V-Ep));
            double b1=a2*sqrt(U01/(U01-Ep));
            double asinb0=asin(b0);
            double asinb1=asin(b1);
            double b2=1/(b0*b0);
            double Z=-2*a0*sqrt(U00/E0)*(sqrt(b2-1)+asinb0*b2);
            b2=1/(b1*b1);
            Z=Z-2*a0*a2*sqrt(U01/E0)*((3.14159/2-asinb1)*b2-sqrt(b2-1));
            g=exp(Z);
        }
        if(g<0.5) this->gTun+=g;
        else this->gOv+=g;
        G0+=g;
        EE-=Ey;
        Ep-=Ey;
    }
    return G0;
}

double MainWindow::Vbarrier(double r)
{
    double BB,rr;
    const double RMIN=190;//nm used for report
    const double RMAX=410;//nm
    rr=0.5*(RMIN+r*(RMAX-RMIN))/Delta_r;
    rr=rr*rr;
    BB=(1+rr);
    BB=2./BB;
    double aa1=(1000.*sqrt(1.25)-350)/Delta_r;
    aa1=aa1*aa1;
    aa1=4/(1+aa1);

//    double coeff = 0.01532*(this->T-2.2098)*(this->T-2.2098) + 0.4775;
    double coeff = 0.00974*(this->T-1.976)*(this->T-1.976) + 0.491;
//    coeff = 0.53;
    return (aa1+BB)*this->U*(1.+coeff*(this->U-Vg0)/170.);
//    return (aa1+BB)*this->U*(1.+0.53*(this->U-Vg0)/170.);
//    return (aa1+BB)*this->U*(1.+0.5*(this->U-Vg0)/170.);
//    return (aa1+BB)*this->U*(1.+0.7*(this->U-Vg0)/170.);
//    return (aa1+BB)*this->U*(1.+0.5*(this->U-Vg0)/(this->Umax-Vg0));
//    return (aa1+BB)*this->U*(1.+(this->U-Vg0)/(this->Umax-Vg0));
}

double MainWindow::Vdot(void)
{
    double x1,y1,r1;
    x1=500;
    y1=500;
    r1=sqrt(x1*x1+y1*y1)-350;
    r1=r1/Delta_r;
    return 4*this->U/(1+r1*r1);
}

double MainWindow::singleSigma(double r, double rEx)
{
    double Gtot, E,dE, Emin,csh,sum,alpha,Ec,Uc,V;
    V=Vbarrier(r);
    Uc=Vdot();
    double kT=this->T;
    dE=0.1;
    if(dE>=0.6) dE=0.5;
    Ec=this->EFT;
//    double g0=sedlo(Ec, this->Ey, rEx, V);
    if(kT==0)
    {
        if(Ec>Uc)  Gtot=sedlo(Ec, this->Ey,rEx, V);
        else Gtot=this->CUTOFF_SIGMA;
    }
    else
    {
        Gtot=0;
        double GTunnel=0.;
        double GOver=0.;
        int G_type = this->typeCond;//->currentIndex();
        Emin=Ec-40*kT;
        double sumt=0.;
        double aa=0.25*dE/kT;
        if(Emin<Uc) Emin=Uc+dE;
        for(E=Emin; E<=Ec+40*kT; E+=dE)
        {
            alpha=0.5*(E-Ec)/kT;
            csh=1./cosh(alpha);
            sum=aa*csh*csh;
            sumt=sumt+sum;
            double g=sedlo(E, this->Ey, rEx, V);
            GTunnel+= this->gTun*sum;
            GOver+= this->gOv*sum;
            Gtot+=g*sum;
        }
        if(G_type==1) Gtot=GTunnel;
        if(G_type==2) Gtot=GOver;
        double eps=0;//0.0075*this->U-1.02;
        if(G_type==3) Gtot=GOver+GTunnel*exp(-eps/kT);
    }
    Gtot=Gtot*G_ser/(Gtot+G_ser);

    if(Gtot<CUTOFF_SIGMA)
        return CUTOFF_SIGMA;
    else
        return Gtot;
}
void MainWindow::randomizeSigma_2()
{
    double x1=0, x2;
    std::valarray<double> t( model->nI() );
    VSLStreamStatePtr stream;
    vslNewStream( &stream, VSL_BRNG_MCG31, this->seed );
    vdRngUniform( 0, stream, model->nI(), &model->Sigma[0], 0.0, 1.0 );
    vslDeleteStream( &stream );
    for (int i=0; i < model->nI(); ++i)
    {
        std::pair<int,int> ends = model->ends(i);
        int from = ends.first;
        int to   = ends.second;
        std::pair<double,double> xy0 = model->xy(from);
        std::pair<double,double> xy1 = model->xy(to);

#if 0
        if (xy0.first==0 && xy1.first==0
            || xy0.first==0 && xy1.first==1
            || xy0.first==1 && xy1.first==0
            || xy0.first==1 && xy1.first==1
            || xy0.first==model->xmax() && xy1.first==model->xmax()
            || xy0.first==model->xmax()-1 && xy1.first==model->xmax()
            || xy0.first==model->xmax()   && xy1.first==model->xmax()-1
            || xy0.first==model->xmax()-1 && xy1.first==model->xmax()-1
            )
        {
            model->Sigma[i]=this->sigmaU;
        }
        else
#endif
        {
            x1=model->Sigma[i];
            x2=this->Ex;
            if (x1 < this->r_c) model->Sigma[i] = CUTOFF_SIGMA;
            else model->Sigma[i] = singleSigma(x1,x2);
        }
    }
}

void MainWindow::randomizeSigma_0()
{
    double x1=0;
	    VSLStreamStatePtr stream;
    vslNewStream( &stream, VSL_BRNG_MCG31, this->seed );
    vdRngUniform( 0, stream, model->nI(), &model->Sigma[0], 0.0, 1.0 );
    vslDeleteStream( &stream );
    for (int i=0; i < model->nI(); ++i)
    {
        std::pair<int,int> ends = model->ends(i);
        int from = ends.first;
        int to   = ends.second;
        std::pair<double,double> xy0 = model->xy(from);
        std::pair<double,double> xy1 = model->xy(to);

#if 0
        if (xy0.first==0 && xy1.first==0
            || xy0.first==0 && xy1.first==1
            || xy0.first==1 && xy1.first==0
            || xy0.first==1 && xy1.first==1
            || xy0.first==model->xmax() && xy1.first==model->xmax()
            || xy0.first==model->xmax()-1 && xy1.first==model->xmax()
            || xy0.first==model->xmax()   && xy1.first==model->xmax()-1
            || xy0.first==model->xmax()-1 && xy1.first==model->xmax()-1
            )
        {
            model->Sigma[i]=this->sigmaU;
        }
        else
#endif
        {
            x1=model->Sigma[i];
            if (x1 > this->r_c)
                model->Sigma[i] = this->CUTOFF_SIGMA;
            else
                model->Sigma[i] = exp(-this->kappa*x1);
            if (model->Sigma[i]<this->CUTOFF_SIGMA)
                model->Sigma[i] = this->CUTOFF_SIGMA;
        }
    }
}

void MainWindow::randomizeSigma_1()
{
    double x1=0;
    VSLStreamStatePtr stream;
    vslNewStream( &stream, VSL_BRNG_MCG31, this->seed );
    vdRngUniform( 0, stream, model->nI(), &model->Sigma[0], 0.0, 1.0 );
    vslDeleteStream( &stream );
    for (int i=0; i < model->nI(); ++i)
    {
        std::pair<int,int> ends = model->ends(i);
        int from = ends.first;
        int to   = ends.second;
        std::pair<double,double> xy0 = model->xy(from);
        std::pair<double,double> xy1 = model->xy(to);

#if 0
        if (xy0.first==0 && xy1.first==0
            || xy0.first==0 && xy1.first==1
            || xy0.first==1 && xy1.first==0
            || xy0.first==1 && xy1.first==1
            || xy0.first==model->xmax() && xy1.first==model->xmax()
            || xy0.first==model->xmax()-1 && xy1.first==model->xmax()
            || xy0.first==model->xmax()   && xy1.first==model->xmax()-1
            || xy0.first==model->xmax()-1 && xy1.first==model->xmax()-1
            )
        {
            model->Sigma[i]=this->sigmaU;
        }
        else
#endif
        {
            x1=model->Sigma[i];
            if (x1 < this->r_c)
                model->Sigma[i] = CUTOFF_SIGMA;
            else
                model->Sigma[i] =1;
        }
    }
}

void MainWindow::selectSigma(int i)
{
    this->setModel();
    switch(i)
    {
    case 0: /* random Sigma_j=exp(-kappa*r_j) */
        this->randomizeSigma_0();
        break;
    case 1: /* uniform Sigma=1 or random Sigma_j=[0,1]*/
        this->randomizeSigma_1();
        break;
    case 2: /* adjusted to experiment */
        this->randomizeSigma_2();
        break;

    default:
        break;
    }
}

//MYARRAY<int> MainWindow::index_for_sorted_CondDist() const
//{
//    MYARRAY<int> res(this->nJ);
//    for (int i=0; i<res.size(); ++i)
//    {
//        res[i] = i;
//    }
//
//    struct MyLessThan
//    {
//        const MYARRAY<double>& w;
//        MyLessThan(const  MYARRAY<double>& w_) : w(w_) {}
//        bool operator()(int a, int b) const { return w[a] < w[b]; }
//    };
//
//    std::sort(res.begin(),res.end(),MyLessThan(CondDist));
//
//    return res;
//}

//MYARRAY<int> MainWindow::index_for_sorted_NumJDist() const
//{
//    MYARRAY<int> res(this->nJ);
//    for (int i=0; i<res.size(); ++i)
//    {
//        res[i] = i;
//    }
//
//    struct MyLessThan
//    {
//        const MYARRAY<double>& w;
//        MyLessThan(const MYARRAY<double>& w_) : w(w_) {}
//        bool operator()(int a, int b) const { return w[a] < w[b]; }
//    };
//
//    std::sort(res.begin(),res.end(),MyLessThan(NumJDist));
//
//    return res;
//}

