#ifndef OPERATORS_H
#define OPERATORS_H

#include "Operator.h"
#include "Site.h"
#include <vector>

class Operators
{
public:

    static void SetOperators(vector<Operator>& op) { operators=op; }
    static void SetIntNeighbors(vector<vector<int>>& op) { int_neighbors=op; }
    static void SetOpNeighbors(vector<vector<Operator>>& op) { op_neighbors=op; }

    static vector<Operator>& GetOperators() { return operators; }
    static vector<int>& GetIntNeighbors(int oper) { return int_neighbors.at(oper); }
    static vector<Operator>& GetOpNeighbors(int oper) { return op_neighbors.at(oper); }

    static int Conjugate(int op)
    {
        return op+(1-2*(op%2));
    }

    static int Transformation(Operator op)
    {
        int i=op.GetType()*operators.size()/7;
        while (!((operators.at(i))==op)) { i++; }

        return i;
    }

    static Operator& Transformation(int op)
    {
        return operators.at(op);
    }

private:

    static vector<Operator> operators;
    static vector<vector<int>> int_neighbors;
    static vector<vector<Operator>> op_neighbors;

};

#endif // OPERATORS_H
