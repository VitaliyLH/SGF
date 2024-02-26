#ifndef UPDATEOPERATORS_H
#define UPDATEOPERATORS_H

#include "Operators.h"
#include "MTRand.h"
#include "State.h"
#include "Green.h"
#include "OffSetOperator.h"

#include <map>
#include <cassert>

class UpdateOperators
{
public:
    UpdateOperators();
    ~UpdateOperators();

    UpdateOperators(const State &stL, const State &stR);

    double GetWeightL(int Offset, int type);
    double GetWeightR(int Offset, int type);

    void UpdateL(int oper, State &stL, State &stR);
    void UpdateR(int oper, State &stL, State &stR);

    int ChooseOperatorL(int offset, int type, MTRand& mt);
    int ChooseOperatorR(int offset, int type, MTRand& mt);

private:

    map<int, OffSetOperator> Initialize(OffSetOperator ofp0, OffSetOperator ofp1, int offset);

    void ActiveHoppings(State stL, State stR);

    void RemoveL(int oper);
    void RemoveR(int oper);

    vector<map<int,OffSetOperator>> l;
    vector<map<int,OffSetOperator>> r;

    unsigned int size_one;

    vector<int> opsL;
    vector<int> opsR;
};

#endif // UPDATEOPERATORS_H
