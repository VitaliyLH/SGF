#include "Probabilities.h"

Probabilities::Probabilities()
{

}



Probabilities::~Probabilities()
{

}



Probabilities::Probabilities(double T, unsigned int L)
{
    alpha = 0.95;
    beta=1.0/T;
    l=L;
    GT=0.0;
    TG=0.0;
    fL=0.0;
    fR=0.0;
    qL=0.0;
    qR=0.0;
    Renormalization=0.0;
}



void Probabilities::UpdateWeightL(double wGT)
{
    GT=wGT;
}



void Probabilities::UpdateWeightR(double wTG)
{
    TG=wTG;
}



double Probabilities::ExponentialHelper(double r,double L)
{
    double u=r*(1-exp(L));
    double s=u/(2-u);
    double ss=s*s;
    double R=((((((ss*Lg7+Lg6)*ss+Lg5)*ss+Lg4)*ss+Lg3)*ss+Lg2)*ss+Lg1)*ss;
    double a= (L==0) ? exp(L/2) : 2*exp(L/2)*sinh(L/2)/L;
    return r*a*(2+R)/(2-u);
}



// Returns a number in the interval [0,1] with distribution exp(A x)
double Probabilities::Exponential(MTRand& mt, double L)
{
    // Do not dare to touch this before you read the documentation above.
    const double T=0.346629;
    unsigned long n=1;
    unsigned long i=0;

    while(fabs(L)>=T) {
        L/=2;
        i<<=1;
        double r=mt.randDblExc();
        if(exp(-L)<=r/(1-r))
            i+=1;
        n<<=1;
    }
    double r=mt.randDblExc();
    if(L>0)
        return (i+ExponentialHelper(r,L))/n;
    else
        return (i+1-ExponentialHelper(1-r,-L))/n;
}



// Returns a number in the interval [0,1] with distribution exp(A x)
double Probabilities::NZExponential(MTRand &mt, double A)
{
    double result=0;
    do { result=Exponential(mt,A); }
    while(result==0);
    return result;
}



double Probabilities::Time(MTRand &mt, double dV, double dt)
{
    return dt*NZExponential(mt, beta*dV*dt);
}



void Probabilities::EvaluateRunningPars(double dV, double dt, bool is_empty)
{
    double expL=0;
    double expR=0;

    if (!is_empty)
    {
        expL=(fabs(beta*dt*dV)==0 ? 1.0/(beta*dt) : dV/(1-exp(-beta*dt*dV)));
        expR=(fabs(beta*dt*dV)==0 ? 1.0/(beta*dt) : -dV/(1-exp(beta*dt*dV)));
    }

    fL=(GT+expL);
    fR=(TG+expR);

    qL=fL*(1-alpha*fmin(1,fR/fL));
    qR=fR*(1-alpha*fmin(1,fL/fR));

    Renormalization = qL+qR;
}



bool Probabilities::ChooseDirectionRight(MTRand& mt)
{
    return qR > mt.randDblExc()*Renormalization;
}



bool Probabilities::ChooseCreateRight(MTRand& mt)
{
    return TG > mt.randDblExc()*fR;
}



bool Probabilities::ChooseCreateLeft(MTRand& mt)
{
    return GT > mt.randDblExc()*fL;
}



bool Probabilities::KeepGoingRight(MTRand& mt)
{
    return alpha*fmin(1,fR/fL) > mt.randDblExc();
}



bool Probabilities::KeepGoingLeft(MTRand& mt)
{
    return alpha*fmin(1,fL/fR) > mt.randDblExc();
}
