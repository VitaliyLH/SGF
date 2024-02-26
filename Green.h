#ifndef GREEN_H
#define GREEN_H

#include "CircularTime.h"
#include "State.h"
#include "Site.h"
#include <math.h>
#include <vector>
#include <iterator>

class Green
{
public:

    static CircularTime Time() { return ti; }

    static void Time(CircularTime t) { ti = t; }

    static void SetCutOff(int cut) { cutoff=cut; }

    static double Me(int n)
    {
        double result=1.0;
        if(n < cutoff && n!=0) result=1.0/(Site::L()*Site::L());
        if(n==cutoff) result=1.0/(Site::L()*Site::L()*Site::L()*Site::L());
        if(n > cutoff) result=pow(1.0/(Site::L()*Site::L()),1.0*(n-1));
        return result;

//        double result=1.0;
//        if(n<=cutoff && n!=0) result=1.0/(Site::L()*Site::L());
//        if(n>cutoff) result=pow(1.0/(Site::L()*Site::L()),1.0*n);
//        return result;
    }


    static double Offset(State &stL, State &stR)
    {
        double offset = 0;

        for (unsigned int i=0;i < stL.GetState().size();i++)
        {
            offset += fabs(stL.SiteState(i).first - stR.SiteState(i).first);
            offset += fabs(stL.SiteState(i).second - stR.SiteState(i).second);
        }

        return offset;
    }



    static double UpdateOffset(Operator &op, State &stL, State &stR)
    {
        double offset = 0;

        offset += fabs(stL.SiteState(op.GetSsite()).first-stR.SiteState(op.GetSsite()).first);

        offset += fabs(stL.SiteState(op.GetTsite()).first-stR.SiteState(op.GetTsite()).first);

        offset += fabs(stL.SiteState(op.GetSsite()).second-stR.SiteState(op.GetSsite()).second);

        offset += fabs(stL.SiteState(op.GetTsite()).second-stR.SiteState(op.GetTsite()).second);

        return offset;
    }



    static double TermOffset(Operator &op, State &stL, State &stR)
    {
        double toffset = 0;

        toffset += fabs(
                    stL.SiteState(op.GetSsite()).first
                    -
                    stR.SiteState(op.GetSsite()).first
                    -
                    op.GetPOffset()
                    );

        toffset += fabs(
                    stL.SiteState(op.GetTsite()).first
                    -
                    stR.SiteState(op.GetTsite()).first
                    +
                    op.GetPOffset()
                    );

        toffset -= fabs(
                    stL.SiteState(op.GetSsite()).first
                    -
                    stR.SiteState(op.GetSsite()).first
                    );

        toffset -= fabs(
                    stL.SiteState(op.GetTsite()).first
                    -
                    stR.SiteState(op.GetTsite()).first
                    );


        toffset += fabs(
                    stL.SiteState(op.GetSsite()).second
                    -
                    stR.SiteState(op.GetSsite()).second
                    -
                    op.GetSOffsetS()
                    );

        toffset += fabs(
                    stL.SiteState(op.GetTsite()).second
                    -
                    stR.SiteState(op.GetTsite()).second
                    +
                    op.GetSOffsetT()
                    );

        toffset -= fabs(
                    stL.SiteState(op.GetSsite()).second
                    -
                    stR.SiteState(op.GetSsite()).second
                    );

        toffset -= fabs(
                    stL.SiteState(op.GetTsite()).second
                    -
                    stR.SiteState(op.GetTsite()).second
                    );

        return toffset;
    }

private:

    static int cutoff;
    static CircularTime ti;
};

#endif // GREEN_H
