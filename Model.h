#ifndef MODEL_H
#define MODEL_H
#include "Operator.h"
#include "State.h"
#include "Site.h"

#include "MTRand.h"

class Model
{
public:

    static void SetDelta(double Delta) { delta=Delta; }
    static void SetV(double V) { v=V; }
    static void SetJ(double J) { j=J; }
    static void SetTp(double Tp) { tp=Tp; }
    static void SetTpn(double Tpn) { tpn=Tpn; }
    static void SetTn(double Tn) { tn=Tn; }
    static void SetT2(double T2) { t2=T2; }

    static double SiteEnergy(double State)
    {
        return delta*State*State;
    }

    static double JBondEnergy(double Sstate , double Tstate)
    {
        return j*Sstate*Tstate;
    }

    static double VBondEnergy(double Sstate , double Tstate)
    {
        return v*Sstate*Tstate;
    }

    static double UpdateEnergy(vector<Operator> &op, State &stO, State &stN)
    {
        double en1 = 0.0;
        double en2 = 0.0;

        en1 += (SiteEnergy(stN.SiteState(op.at(0).GetSsite()).first)-SiteEnergy(stO.SiteState(op.at(0).GetSsite()).first));
        en1 += (SiteEnergy(stN.SiteState(op.at(0).GetTsite()).first)-SiteEnergy(stO.SiteState(op.at(0).GetTsite()).first));

        en1 += (SiteEnergy(stN.SiteState(op.at(1).GetSsite()).first)-SiteEnergy(stO.SiteState(op.at(1).GetSsite()).first));
        en1 += (SiteEnergy(stN.SiteState(op.at(1).GetTsite()).first)-SiteEnergy(stO.SiteState(op.at(1).GetTsite()).first));

        en1*=0.5;

        for (unsigned int i=0;i < op.size();i++) en2 += (
                    VBondEnergy(stN.SiteState(op.at(i).GetSsite()).first, stN.SiteState(op.at(i).GetTsite()).first)
                    +
                    JBondEnergy(stN.SiteState(op.at(i).GetSsite()).second, stN.SiteState(op.at(i).GetTsite()).second)
                    -
                    VBondEnergy(stO.SiteState(op.at(i).GetSsite()).first, stO.SiteState(op.at(i).GetTsite()).first)
                    -
                    JBondEnergy(stO.SiteState(op.at(i).GetSsite()).second, stO.SiteState(op.at(i).GetTsite()).second)
                    );

        en2*=0.5;

        return en1+en2;
    }

    static double Energy(State &st)
    {
        double en = 0.0;

        for (int i=0;i < Site::L()*Site::L();i++)  en += SiteEnergy(st.SiteState(i).first);

        for (unsigned int y=0;y < Site::L();y++)
            for (unsigned int x=0;x < Site::L();x++) {

                en += JBondEnergy(st.SiteState(Site::GetSite(x,y)).second, st.SiteState(Site::GetSite(x+1,y)).second);
                en += JBondEnergy(st.SiteState(Site::GetSite(x,y)).second, st.SiteState(Site::GetSite(x,y+1)).second);

                en += VBondEnergy(st.SiteState(Site::GetSite(x,y)).first, st.SiteState(Site::GetSite(x+1,y)).first);
                en += VBondEnergy(st.SiteState(Site::GetSite(x,y)).first, st.SiteState(Site::GetSite(x,y+1)).first);
            }

        return en;
    }

    static double DeltaEnergy(State &stL, State &stR)
    {
        return Energy(stL)-Energy(stR);
    }

private:

    static double delta;
    static double v;
    static double j;
    static double tp;
    static double tpn;
    static double tn;
    static double t2;
};

#endif // MODEL_H
