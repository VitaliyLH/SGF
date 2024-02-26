#ifndef PROCEDURE_H
#define PROCEDURE_H
#include <vector>
#include <algorithm>
#include "Green.h"
#include "Model.h"
#include "Site.h"
#include "Operator.h"
#include "Operators.h"
#include "CircList.h"
#include "MTRand.h"
#include "Probabilities.h"
#include "UpdateOperators.h"
#include "rapidjson/document.h"
#include <QTextStream>
#include <QFile>


class Procedure
{
public:
    Procedure(const rapidjson::Value&parameters);
    Procedure(const rapidjson::Value &parameters, const rapidjson::Value &configuration);
    ~Procedure();

    Probabilities &GetProbabilities() { return p; }

    State &GetState() { return stL; }

    CircList &GetKinks() { return c; }

    void DirectUpdate();

    vector<double>& GetWindingNumber() { return winding_number; }

    void WriteConfiguration(QFile &File);

private:

    vector<vector<int>> InitializeIntNeighbors(vector<Operator> &operators);
    vector<vector<Operator>> InitializeOpNeighbors(vector<Operator> &operators);

    void FinalizeIntNeighbors(vector<vector<int>> &int_neighbors);

    vector<Operator> Initialize(int type, int poff, double soffs, double sofft);

    void ActiveHoppings();

    int ChooseL();
    int ChooseR();

    void UpdateWeights();

    Probabilities p;

    CircList c;

    MTRand mt;

    State stL;
    State stR;

    State stLO;
    State stRO;

    int Uop2L;
    int Uop2R;

    double GT;
    double TG;

    double green_offset;

    UpdateOperators u;

    vector<double> me;

    vector<double> winding_number;
};

#endif // PROCEDURE_H
