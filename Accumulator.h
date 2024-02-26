#ifndef ACCUMULATOR_H
#define ACCUMULATOR_H

#include <math.h>

using namespace std;

class Accumulator {
public:
  Accumulator() : n(0), mean(0), moment(0) {}
  void Push(const double& x) {
    ++n;
    double delta = x - mean;
    double R = delta / n;
    mean = mean + R;
    moment = moment + delta * (x - mean);
  }
  double Sigma() {
    return sqrt(moment / (n - 1));
  }
  double Average() {
    return mean;
  }

private:
  int n;
  double mean;
  double moment;

};

#endif // ACCUMULATOR_H
