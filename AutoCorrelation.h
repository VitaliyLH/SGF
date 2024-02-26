#ifndef AUTOCORRELATION_H
#define AUTOCORRELATION_H
#include <vector>
#include <math.h>

#include <QTextStream>
#include <QFile>

#include "Accumulator.h"

using namespace std;

class AutoCorrelation
{
public:
    AutoCorrelation();
    ~AutoCorrelation();

    AutoCorrelation(unsigned int NumberBins);

    void Add(unsigned int BinNumber, double V);

    double GetAverage() { return accumulator.Average(); }

    double GetError() { return accumulator.Sigma(); }

    void Print(QTextStream &stream);

    void Count(unsigned int BinNumber, double N);

private:

    Accumulator accumulator;

    vector<double> value;
};

#endif // AUTOCORRELATION_H
