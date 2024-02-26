#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include <cassert>

#include <QTextStream>
#include <QFile>

#include "AutoCorrelation.h"
#include "Probabilities.h"
#include "Procedure.h"

class Measurement
{
public:
    Measurement();
    ~Measurement();

    Measurement(unsigned int Size, double T, unsigned int BinsLenght, unsigned int NumberBins, double MaxError);

    bool isDone();
    bool isDoneKinetic();
    bool isDonePotential();
    void Coefficient(unsigned int &Times);
    void SetBinsLenght(unsigned int BinsLenght) { bins_lenght=BinsLenght; }

    void Make(Procedure& p);
    pair<double,double> GetWindingNumber(Procedure& p, unsigned int size);
    vector<pair<int,int>> GetTermDistribution(Procedure& p, unsigned int size);

    void PrintOutput(QFile &File);
    void Printsz(QFile &File);
    void PrintSz(QFile &File);
    void PrintP0(QFile &File);
    void PrintLocalsz(QFile &File);
    void PrintLocalSz(QFile &File);
    void PrintStiffness(QFile &File);
    void PrintStructureFactorsz(QFile &File);
    void PrintStructureFactorSz(QFile &File);

private:

    double norm;
    double max_error;
    unsigned int count;
    double temperature;
    unsigned int count_bad;
    unsigned int bin_number;
    unsigned int bins_lenght;
    unsigned int number_sites;
    double coefficient_potential;
    double coefficient_kinetic;

    //    void GetStates(Procedure& p);

    //    vector<State> state;
    //    vector<double> deltatau;

    vector<double> normalization;

    vector<AutoCorrelation> auto_mean;
    vector<AutoCorrelation> local_sz;
    vector<AutoCorrelation> local_Sz;
};

#endif // MEASUREMENT_H
