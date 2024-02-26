#ifndef OPERATOR_H
#define OPERATOR_H

#include <vector>

using namespace std;

class Operator
{
public:
    Operator();
    ~Operator();

    Operator(int type, int poff, double soffs, double sofft, int ssite, int tsite);

    bool operator==(const Operator& op) { return (op.pof == pof && op.sofs == sofs && op.soft == soft && op.ty == ty && op.ss == ss && op.ts == ts); }

    int& GetType() { return ty; }
    int& GetPOffset() { return pof ;}
    double& GetSOffsetS() { return sofs;}
    double& GetSOffsetT() { return soft;}
    int& GetSsite() { return ss; }
    int& GetTsite() { return ts; }

    int MatrixElement(pair<double, double> sstate, pair<double, double> tstate);

private: 

   int ty;
   int pof;
   double sofs;
   double soft;
   int ss;
   int ts;

};

#endif // OPERATOR_H
