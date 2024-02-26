#ifndef PROBABILITIES_H
#define PROBABILITIES_H

#include "MTRand.h"

class Probabilities
{
public:
    Probabilities();
    ~Probabilities();

    Probabilities(double T, unsigned int L);

    void UpdateWeightL(double wGT);
    void UpdateWeightR(double wTG);

    double GetBeta() { return beta; }

    double Time(MTRand& mt, double dV, double dt);

    bool ChooseDirectionRight(MTRand& mt);

    bool ChooseCreateRight(MTRand& mt);

    bool ChooseCreateLeft(MTRand& mt);

    bool KeepGoingRight(MTRand& mt);

    bool KeepGoingLeft(MTRand& mt);

    void EvaluateRunningPars(double dV , double dt , bool is_empty);

    double GetRenormalization() { return Renormalization; }

private:

    double GT;
    double TG;

    double alpha;
    double beta;
    unsigned int l;

    double fL;
    double fR;
    double qL;
    double qR;

    double Renormalization;

    // This variables are used by Exponential(double)
    double Lg1 = 6.666666666666735130e-01;
    double Lg2 = 3.999999999940941908e-01;
    double Lg3 = 2.857142874366239149e-01;
    double Lg4 = 2.222219843214978396e-01;
    double Lg5 = 1.818357216161805012e-01;
    double Lg6 = 1.531383769920937332e-01;
    double Lg7 = 1.479819860511658591e-01;

    double ExponentialHelper(double r,double L);

    // Returns a number in the interval [0,1[ with distribution exp(A x)
    double Exponential(MTRand &mt, double L);

    // Returns a number in the interval ]0,1[ with distribution exp(A x)
    double NZExponential(MTRand& mt, double A);
};

#endif // PROBABILITIES_H
