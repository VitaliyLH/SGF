#ifndef CIRCULARTIME_H
#define CIRCULARTIME_H

#include <cmath>

using namespace std;

class CircularTime {
public:
  explicit CircularTime(double t=1.0) : ti(t+1.0-ceil(t)) {}
  double  Time() const { return ti; }
  CircularTime operator+(const CircularTime &a) {return CircularTime(ti+a.Time());}
  CircularTime operator-(const CircularTime &a) {return CircularTime(ti-a.Time());}
private:
 double ti;
};



#endif // CIRCULARTIME_H
