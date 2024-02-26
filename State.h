#ifndef STATE_H
#define STATE_H

#include <vector>
#include "MTRand.h"
#include "Operator.h"

class State
{
public:
    State();
    ~State();

    State(int L, int x);

    vector<pair<double,double>> GetState() { return st; }

    void SetState( vector<pair<double,double>> &state) { st=state; }

    void ChangeState(Operator &op);

    void Conjugate();

    pair<double,double> SiteState(int site) {return st[site]; }

private:

    vector<pair<double,double>> st;

};

#endif // STATE_H
