#include "AutoCorrelation.h"
#include <iostream>

AutoCorrelation::AutoCorrelation()
{

}



AutoCorrelation::~AutoCorrelation()
{

}



AutoCorrelation::AutoCorrelation(unsigned int NumberBins)
{
  value.resize(NumberBins);
}



void AutoCorrelation::Print(QTextStream &stream)
{
    stream.setFieldWidth(5);
    stream.setFieldAlignment(QTextStream::AlignLeft);
    stream.setRealNumberPrecision(10);

    stream << accumulator.Average() << "\t" << accumulator.Sigma();
}



void AutoCorrelation::Add(unsigned int BinNumber, double V)
{
    value[BinNumber]+=V;
}



void AutoCorrelation::Count(unsigned int BinNumber, double N)
{
    accumulator.Push(value.at(BinNumber)/N);
}
