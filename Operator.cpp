#include "Operator.h"

Operator::Operator()
{

}

Operator::~Operator()
{

}

Operator::Operator(int type, int poff, double soffs, double sofft, int ssite, int tsite)
{
    ty=type;
    pof=poff;
    sofs=soffs;
    soft=sofft;
    ss=ssite;
    ts=tsite;
}



int Operator::MatrixElement(pair<double,double> sstate, pair<double,double> tstate)
{
    int me=0;

    bool pcontp1 = (sstate.first==0 && tstate.first==1);
    bool scontp1= (sstate.second==-0.5 && tstate.second==0);

    bool pcontp2 = (sstate.first==0 && tstate.first==1);
    bool scontp2= (sstate.second==0.5 && tstate.second==0);

    bool pcontpn1 = (sstate.first==0 && tstate.first==0);
    bool scontpn1 = (sstate.second==-0.5 && tstate.second==0.5);

    bool pcontpn2 = (sstate.first==0 && tstate.first==0);
    bool scontpn2 = (sstate.second==0.5 && tstate.second==-0.5);

    bool pcontpn = (sstate.first==-1 && tstate.first==1);
    bool scontpn = (sstate.second==0 && tstate.second==0);

    bool pcontn1 = (sstate.first==-1 && tstate.first==0);
    bool scontn1 = (sstate.second==0 && tstate.second==0.5);

    bool pcontn2 = (sstate.first==-1 && tstate.first==0);
    bool scontn2 = (sstate.second==0 && tstate.second==-0.5);

    if (ty==0 && (pcontp1 && scontp1)) me=1;
    if (ty==1 && (pcontp2 && scontp2)) me=1;

    if (ty==2 && ((pcontpn1 && scontpn1) || (pcontpn && scontpn))) me=1;
    if (ty==3 && ((pcontpn2 && scontpn2) || (pcontpn && scontpn))) me=1;

    if (ty==4 && (pcontn1 && scontn1)) me=1;
    if (ty==5 && (pcontn2 && scontn2)) me=1;

    if (ty==6 && pcontpn && scontpn) me=1;

    return me;
}
