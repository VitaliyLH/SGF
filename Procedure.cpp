#include "Procedure.h"

Procedure::Procedure(const rapidjson::Value& parameters)
{
    unsigned int L=parameters["L"].GetUint();
    double T=parameters["T"].GetDouble();
    int x=parameters["x"].GetInt();
    double Delta=parameters["Delta"].GetDouble();
    double V=parameters["V"].GetDouble();
    double J=parameters["J"].GetDouble();
    double tp=parameters["tp"].GetDouble();
    double tpn=parameters["tpn"].GetDouble();
    double tn=parameters["tn"].GetDouble();
    double t2=parameters["t2"].GetDouble();
    unsigned int Seed=parameters["Seed"].GetUint();
    int CutOff=parameters["Green"].GetInt();

    vector<vector<int>> site;

    site.resize(L);

    for (int x = 0; x < L; x++)
        site[x].resize(L, L * L);

    int count = 0;

    for (int y = 0; y < L; y++)
        for (int x = 0; x < L; x++)
            site[x][y] = count++;

    Site::SetL(L);
    Site::SetSite(site);

    Green::SetCutOff(CutOff);

    Model::SetDelta(Delta);
    Model::SetV(V);
    Model::SetJ(J);
    Model::SetTp(tp);
    Model::SetTpn(tpn);
    Model::SetTn(tn);
    Model::SetT2(t2);

    vector<Operator> op;
    vector<Operator> op0;
    vector<Operator> op1;
    vector<Operator> op2;
    vector<Operator> op3;
    vector<Operator> op4;
    vector<Operator> op5;
    vector<Operator> op6;

    op0=Initialize(0,1,0.5,0.5);
    op1=Initialize(1,1,-0.5,-0.5);
    op2=Initialize(2,1,0.5,0.5);
    op3=Initialize(3,1,-0.5,-0.5);
    op4=Initialize(4,1,0.5,0.5);
    op5=Initialize(5,1,-0.5,-0.5);
    op6=Initialize(6,2,0,0);

    op.insert(op.begin(), op0.begin(), op0.end());
    op.insert(op.end(), op1.begin(), op1.end());
    op.insert(op.end(), op2.begin(), op2.end());
    op.insert(op.end(), op3.begin(), op3.end());
    op.insert(op.end(), op4.begin(), op4.end());
    op.insert(op.end(), op5.begin(), op5.end());
    op.insert(op.end(), op6.begin(), op6.end());

    vector<vector<int>> int_neighbors;
    vector<vector<Operator>> op_neighbors;

    Operators::SetOperators(op);

    int_neighbors=InitializeIntNeighbors(op);
    op_neighbors=InitializeOpNeighbors(op);

    FinalizeIntNeighbors(int_neighbors);

    Operators::SetIntNeighbors(int_neighbors);
    Operators::SetOpNeighbors(op_neighbors);

    stL=State(L,x);
    stR=stL;

    stLO=stL;
    stRO=stL;

    mt=MTRand(Seed);
    p=Probabilities(T,L);
    u=UpdateOperators(stL,stR);


    Uop2L=0;
    Uop2R=0;
    GT=0;
    TG=0;
    green_offset=0;
    me.resize(7);
    me.at(0)=2.0*tp;
    me.at(1)=2.0*tp;
    me.at(2)=2.0*tpn;
    me.at(3)=2.0*tpn;
    me.at(4)=2.0*tn;
    me.at(5)=2.0*tn;
    me.at(6)=4.0*t2;
    vector<double> wn(2,0);
    winding_number=wn;
}


Procedure::Procedure(const rapidjson::Value& parameters, const rapidjson::Value& configuration)
{
    unsigned int L=parameters["L"].GetUint();
    double T=parameters["T"].GetDouble();
    double Delta=parameters["Delta"].GetDouble();
    double V=parameters["V"].GetDouble();
    double J=parameters["J"].GetDouble();
    double tp=parameters["tp"].GetDouble();
    double tpn=parameters["tpn"].GetDouble();
    double tn=parameters["tn"].GetDouble();
    double t2=parameters["t2"].GetDouble();
    int CutOff=parameters["Green"].GetInt();

    double wnx=configuration["wnx"].GetDouble();
    double wny=configuration["wny"].GetDouble();

    vector<double> wn(2,0);
    winding_number=wn;

    winding_number[0]=wnx;
    winding_number[1]=wny;

    double gtime=configuration["greentime"].GetDouble();

    bool IsArray1=configuration["state"].IsArray();
    if (IsArray1)
    {
        vector<pair<double,double>> st;
        const rapidjson::Value& state=configuration["state"];
        for (unsigned long i=0; i < configuration["numbersites"].GetUint(); i++)
        {
            const rapidjson::Value& s = state[i];

            unsigned long l=0;
            unsigned long m=1;

            double ps=s[l].GetDouble();
            double sp=s[m].GetDouble();

            st.push_back(make_pair(ps,sp));
        }
        stL.SetState(st);
        stR=stL;

        stLO=stL;
        stRO=stL;
    }

    MTRand::uint32 randstate[ mt.SAVE ];

    bool IsArray2=configuration["randstate"].IsArray();
    if (IsArray2)
    {
        const rapidjson::Value& rand=configuration["randstate"];
        for (unsigned long i=0; i < mt.SAVE; i++)
        {
            const rapidjson::Value& rs = rand[i];
            randstate[i]=rs.GetUint();
        }
    }

    if (configuration["numberkinks"].GetUint()!=0)
    {
        bool IsArray3=configuration["kinks"].IsArray();
        if (IsArray3)
        {
            const rapidjson::Value& ki=configuration["kinks"];
            for (unsigned long i=0; i < configuration["numberkinks"].GetUint(); i++)
            {
                const rapidjson::Value& k = ki[i];

                unsigned long l=0;
                unsigned long m=1;
                double time=k[l].GetDouble();
                int op=k[m].GetInt();

                c.Push(-1,Kink{CircularTime(time),op});
            }
        }
    }

    vector<vector<int>> site;

    site.resize(L);

    for (int x = 0; x < L; x++)
        site[x].resize(L, L * L);

    int count = 0;

    for (int y = 0; y < L; y++)
        for (int x = 0; x < L; x++)
            site[x][y] = count++;

    Site::SetL(L);
    Site::SetSite(site);

    Green::SetCutOff(CutOff);
    Green::Time(CircularTime(gtime));

    Model::SetDelta(Delta);
    Model::SetV(V);
    Model::SetJ(J);
    Model::SetTp(tp);
    Model::SetTpn(tpn);
    Model::SetTn(tn);
    Model::SetT2(t2);

    vector<Operator> op;
    vector<Operator> op0;
    vector<Operator> op1;
    vector<Operator> op2;
    vector<Operator> op3;
    vector<Operator> op4;
    vector<Operator> op5;
    vector<Operator> op6;

    op0=Initialize(0,1,0.5,0.5);
    op1=Initialize(1,1,-0.5,-0.5);
    op2=Initialize(2,1,0.5,0.5);
    op3=Initialize(3,1,-0.5,-0.5);
    op4=Initialize(4,1,0.5,0.5);
    op5=Initialize(5,1,-0.5,-0.5);
    op6=Initialize(6,2,0,0);

    op.insert(op.begin(), op0.begin(), op0.end());
    op.insert(op.end(), op1.begin(), op1.end());
    op.insert(op.end(), op2.begin(), op2.end());
    op.insert(op.end(), op3.begin(), op3.end());
    op.insert(op.end(), op4.begin(), op4.end());
    op.insert(op.end(), op5.begin(), op5.end());
    op.insert(op.end(), op6.begin(), op6.end());

    vector<vector<int>> int_neighbors;
    vector<vector<Operator>> op_neighbors;

    Operators::SetOperators(op);

    int_neighbors=InitializeIntNeighbors(op);
    op_neighbors=InitializeOpNeighbors(op);

    FinalizeIntNeighbors(int_neighbors);

    Operators::SetIntNeighbors(int_neighbors);
    Operators::SetOpNeighbors(op_neighbors);

    p=Probabilities(T,L);
    mt.load(randstate);
    u=UpdateOperators(stL,stR);

    Uop2L=0;
    Uop2R=0;
    GT=0;
    TG=0;
    green_offset=0;
    me.resize(7);
    me.at(0)=2.0*tp;
    me.at(1)=2.0*tp;
    me.at(2)=2.0*tpn;
    me.at(3)=2.0*tpn;
    me.at(4)=2.0*tn;
    me.at(5)=2.0*tn;
    me.at(6)=4.0*t2;
}



Procedure::~Procedure()
{

}



vector<Operator> Procedure::Initialize(int type, int poff, double soffs, double sofft)
{
    vector<Operator> operators;

    for (int y=0; y < Site::L();y++) {
        for (int x=0; x < Site::L();x++) {
            operators.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(x,y), Site::GetSite(x+1,y)));
            operators.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(x+1,y), Site::GetSite(x,y)));

            operators.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(x,y), Site::GetSite(x,y+1)));
            operators.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(x,y+1), Site::GetSite(x,y)));
        }
    }

    return operators;
}


vector<vector<int>> Procedure::InitializeIntNeighbors(vector<Operator> &operators)
{
    vector<vector<int>> int_neighbors;

    for(unsigned int i=0;i<operators.size();i+=2)
    {
        vector<int> int_ops;

        int xi=operators.at(i).GetSsite()%Site::L();
        int yi=operators.at(i).GetSsite()/Site::L();

        int xj=operators.at(i).GetTsite()%Site::L();
        int yj=operators.at(i).GetTsite()/Site::L();

        int type=operators.at(i).GetType();
        int poff=operators.at(i).GetPOffset();
        double soffs=operators.at(i).GetSOffsetS();
        double sofft=operators.at(i).GetSOffsetT();

        int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xj,yj))));
        int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xi,yi))));

        if (yi==yj) {

            //(xj,yj)
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj+1,yj))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj+1,yj), Site::GetSite(xj,yj))));

            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj,yj+1))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj+1), Site::GetSite(xj,yj))));

            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj,yj-1))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj-1), Site::GetSite(xj,yj))));

            //(xi,yi)
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xi-1,yi))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi-1,yi), Site::GetSite(xi,yi))));

            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xi,yi+1))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi+1), Site::GetSite(xi,yi))));

            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xi,yi-1))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi-1), Site::GetSite(xi,yi))));

        } else {

            //(xj,yj)
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj,yj+1))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj+1), Site::GetSite(xj,yj))));

            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj+1,yj))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj+1,yj), Site::GetSite(xj,yj))));

            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj-1,yj))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xj-1,yj), Site::GetSite(xj,yj))));

            //(xi,yi)
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xi,yi-1))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi-1), Site::GetSite(xi,yi))));

            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xi+1,yi))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi+1,yi), Site::GetSite(xi,yi))));

            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xi-1,yi))));
            int_ops.push_back(Operators::Transformation(Operator(type, poff, soffs, sofft, Site::GetSite(xi-1,yi), Site::GetSite(xi,yi))));
        }

        int_neighbors.push_back(int_ops);
        int_neighbors.push_back(int_ops);
    }

    return int_neighbors;
}



vector<vector<Operator>> Procedure::InitializeOpNeighbors(vector<Operator> &operators)
{
    vector<vector<Operator>> op_neighbors;

    for(unsigned int i=0;i<operators.size();i+=2)
    {
        vector<int> int_ops;
        vector<Operator> op_ops;

        int xi=operators.at(i).GetSsite()%Site::L();
        int yi=operators.at(i).GetSsite()/Site::L();

        int xj=operators.at(i).GetTsite()%Site::L();
        int yj=operators.at(i).GetTsite()/Site::L();

        int type=operators.at(i).GetType();
        int poff=operators.at(i).GetPOffset();
        double soffs=operators.at(i).GetSOffsetS();
        double sofft=operators.at(i).GetSOffsetT();

        op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xj,yj)));
        op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xi,yi)));

        if (yi==yj) {

            //(xj,yj)
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj+1,yj)));
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj+1,yj), Site::GetSite(xj,yj)));

            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj,yj+1)));
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj+1), Site::GetSite(xj,yj)));

            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj,yj-1)));
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj-1), Site::GetSite(xj,yj)));

            //(xi,yi)
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xi-1,yi)));
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xi-1,yi), Site::GetSite(xi,yi)));

            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xi,yi+1)));
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi+1), Site::GetSite(xi,yi)));

            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xi,yi-1)));
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi-1), Site::GetSite(xi,yi)));

        } else {

            //(xj,yj)
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj,yj+1)));
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj+1), Site::GetSite(xj,yj)));

            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj+1,yj)));
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj+1,yj), Site::GetSite(xj,yj)));

            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj,yj), Site::GetSite(xj-1,yj)));
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xj-1,yj), Site::GetSite(xj,yj)));

            //(xi,yi)
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi), Site::GetSite(xi,yi-1)));
            op_ops.push_back(Operator(type, poff, soffs, sofft, Site::GetSite(xi,yi-1), Site::GetSite(xi,yi)));

            op_ops.push_back(Operator(type, poff, soffs, sofft,  Site::GetSite(xi,yi), Site::GetSite(xi+1,yi)));
            op_ops.push_back(Operator(type, poff, soffs, sofft,  Site::GetSite(xi+1,yi), Site::GetSite(xi,yi)));

            op_ops.push_back(Operator(type, poff, soffs, sofft,  Site::GetSite(xi,yi), Site::GetSite(xi-1,yi)));
            op_ops.push_back(Operator(type, poff, soffs, sofft,  Site::GetSite(xi-1,yi), Site::GetSite(xi,yi)));

        }

        op_neighbors.push_back(op_ops);
        op_neighbors.push_back(op_ops);
    }

    return op_neighbors;
}



void Procedure::FinalizeIntNeighbors(vector<vector<int>> &int_neighbors)
{
    unsigned int size=int_neighbors.size()/7;

    for (unsigned int i=0; i < size; i++)
    {
        vector<int> int_ops0=int_neighbors.at(i);
        vector<int> int_ops1=int_neighbors.at(i+size);
        vector<int> int_ops2=int_neighbors.at(i+2*size);
        vector<int> int_ops3=int_neighbors.at(i+3*size);
        vector<int> int_ops4=int_neighbors.at(i+4*size);
        vector<int> int_ops5=int_neighbors.at(i+5*size);
        vector<int> int_ops6=int_neighbors.at(i+6*size);

        //type=0
        for (unsigned int j=0; j < int_ops0.size(); j++) int_neighbors.at(i).push_back(int_ops1.at(j));
        for (unsigned int j=0; j < int_ops0.size(); j++) int_neighbors.at(i).push_back(int_ops2.at(j));
        for (unsigned int j=0; j < int_ops0.size(); j++) int_neighbors.at(i).push_back(int_ops3.at(j));
        for (unsigned int j=0; j < int_ops0.size(); j++) int_neighbors.at(i).push_back(int_ops4.at(j));
        for (unsigned int j=0; j < int_ops0.size(); j++) int_neighbors.at(i).push_back(int_ops5.at(j));
        for (unsigned int j=0; j < int_ops0.size(); j++) int_neighbors.at(i).push_back(int_ops6.at(j));

        //type=1
        for (unsigned int j=0; j < int_ops1.size(); j++) int_neighbors.at(i+size).push_back(int_ops0.at(j));
        for (unsigned int j=0; j < int_ops1.size(); j++) int_neighbors.at(i+size).push_back(int_ops2.at(j));
        for (unsigned int j=0; j < int_ops1.size(); j++) int_neighbors.at(i+size).push_back(int_ops3.at(j));
        for (unsigned int j=0; j < int_ops1.size(); j++) int_neighbors.at(i+size).push_back(int_ops4.at(j));
        for (unsigned int j=0; j < int_ops1.size(); j++) int_neighbors.at(i+size).push_back(int_ops5.at(j));
        for (unsigned int j=0; j < int_ops1.size(); j++) int_neighbors.at(i+size).push_back(int_ops6.at(j));

        //type=2
        for (unsigned int j=0; j < int_ops2.size(); j++) int_neighbors.at(i+2*size).push_back(int_ops0.at(j));
        for (unsigned int j=0; j < int_ops2.size(); j++) int_neighbors.at(i+2*size).push_back(int_ops1.at(j));
        for (unsigned int j=0; j < int_ops2.size(); j++) int_neighbors.at(i+2*size).push_back(int_ops3.at(j));
        for (unsigned int j=0; j < int_ops2.size(); j++) int_neighbors.at(i+2*size).push_back(int_ops4.at(j));
        for (unsigned int j=0; j < int_ops2.size(); j++) int_neighbors.at(i+2*size).push_back(int_ops5.at(j));
        for (unsigned int j=0; j < int_ops2.size(); j++) int_neighbors.at(i+2*size).push_back(int_ops6.at(j));

        //type=3
        for (unsigned int j=0; j < int_ops3.size(); j++) int_neighbors.at(i+3*size).push_back(int_ops0.at(j));
        for (unsigned int j=0; j < int_ops3.size(); j++) int_neighbors.at(i+3*size).push_back(int_ops1.at(j));
        for (unsigned int j=0; j < int_ops3.size(); j++) int_neighbors.at(i+3*size).push_back(int_ops2.at(j));
        for (unsigned int j=0; j < int_ops3.size(); j++) int_neighbors.at(i+3*size).push_back(int_ops4.at(j));
        for (unsigned int j=0; j < int_ops3.size(); j++) int_neighbors.at(i+3*size).push_back(int_ops5.at(j));
        for (unsigned int j=0; j < int_ops3.size(); j++) int_neighbors.at(i+3*size).push_back(int_ops6.at(j));

        //type=4
        for (unsigned int j=0; j < int_ops4.size(); j++) int_neighbors.at(i+4*size).push_back(int_ops0.at(j));
        for (unsigned int j=0; j < int_ops4.size(); j++) int_neighbors.at(i+4*size).push_back(int_ops1.at(j));
        for (unsigned int j=0; j < int_ops4.size(); j++) int_neighbors.at(i+4*size).push_back(int_ops2.at(j));
        for (unsigned int j=0; j < int_ops4.size(); j++) int_neighbors.at(i+4*size).push_back(int_ops3.at(j));
        for (unsigned int j=0; j < int_ops4.size(); j++) int_neighbors.at(i+4*size).push_back(int_ops5.at(j));
        for (unsigned int j=0; j < int_ops4.size(); j++) int_neighbors.at(i+4*size).push_back(int_ops6.at(j));

        //type=5
        for (unsigned int j=0; j < int_ops5.size(); j++) int_neighbors.at(i+5*size).push_back(int_ops0.at(j));
        for (unsigned int j=0; j < int_ops5.size(); j++) int_neighbors.at(i+5*size).push_back(int_ops1.at(j));
        for (unsigned int j=0; j < int_ops5.size(); j++) int_neighbors.at(i+5*size).push_back(int_ops2.at(j));
        for (unsigned int j=0; j < int_ops5.size(); j++) int_neighbors.at(i+5*size).push_back(int_ops3.at(j));
        for (unsigned int j=0; j < int_ops5.size(); j++) int_neighbors.at(i+5*size).push_back(int_ops4.at(j));
        for (unsigned int j=0; j < int_ops5.size(); j++) int_neighbors.at(i+5*size).push_back(int_ops6.at(j));

        //type=6
        for (unsigned int j=0; j < int_ops6.size(); j++) int_neighbors.at(i+6*size).push_back(int_ops0.at(j));
        for (unsigned int j=0; j < int_ops6.size(); j++) int_neighbors.at(i+6*size).push_back(int_ops1.at(j));
        for (unsigned int j=0; j < int_ops6.size(); j++) int_neighbors.at(i+6*size).push_back(int_ops2.at(j));
        for (unsigned int j=0; j < int_ops6.size(); j++) int_neighbors.at(i+6*size).push_back(int_ops3.at(j));
        for (unsigned int j=0; j < int_ops6.size(); j++) int_neighbors.at(i+6*size).push_back(int_ops4.at(j));
        for (unsigned int j=0; j < int_ops6.size(); j++) int_neighbors.at(i+6*size).push_back(int_ops5.at(j));
    }
}



void Procedure::WriteConfiguration(QFile &File)
{
    QTextStream stream( &File );

    stream << "{" << "\n";
    stream << '\"' << "configuration" << '\"' << ":{" << "\n";
    stream << '\"' << "numbersites" << '\"' << ":" << stL.GetState().size() << "," << "\n";
    stream << '\"' << "numberkinks" << '\"' << ":" << c.Length() << "," << "\n";
    stream << '\"' << "wnx" << '\"' << ":" << winding_number.at(0) << "," << "\n";
    stream << '\"' << "wny" << '\"' << ":" << winding_number.at(1) << "," << "\n";
    stream << '\"' << "greentime" << '\"' << ":" << Green::Time().Time() << "," << "\n";

    stream << '\"' << "state" << '\"' << ":[" << "\n";
    for (unsigned int i=0; i < stL.GetState().size()-1; i++) stream << "[" << stL.SiteState(i).first << "," << stL.SiteState(i).second << "]" << "," << "\n";
    stream << "[" << stL.SiteState(stL.GetState().size()-1).first << "," << stL.SiteState(stL.GetState().size()-1).second << "]" << "\n";
    stream << "]" << "," << "\n";

    MTRand::uint32 randstate[ mt.SAVE ];
    mt.save(randstate);

    stream << '\"' << "randstate" << '\"' << ":[" << "\n";
    for (unsigned int i=0; i < mt.SAVE-1; i++) stream << randstate[i] << "," << "\n";
    stream << randstate[mt.SAVE-1] << "\n";

    if (c.Length()!=0) stream << "]" << "," << "\n";
    else stream << "]" << "\n";

    if (c.Length()!=0)
    {
        stream << '\"' << "kinks" << '\"' << ":[" << "\n";
        for (unsigned int i=0; i < c.Length()-1; i++) stream << "[" << c[i].ti.Time() << "," << c[i].op << "]" << "," << "\n";
        stream << "[" << c[c.Length()-1].ti.Time() << "," << c[c.Length()-1].op << "]" << "\n";
        stream << "]" << "\n";
    }

    stream << "}" << "\n";
    stream << "}";
}



void Procedure::UpdateWeights()
{
    GT=0;
    TG=0;

    for (int Offset=-4;Offset<=4;Offset+=1)
    {

        GT+=Green::Me(green_offset+Offset)*
                (
                    me.at(0)*u.GetWeightL(Offset,0)+me.at(1)*u.GetWeightL(Offset,1)
                    +
                    me.at(2)*u.GetWeightL(Offset,2)+me.at(3)*u.GetWeightL(Offset,3)
                    +
                    me.at(4)*u.GetWeightL(Offset,4)+me.at(5)*u.GetWeightL(Offset,5)
                    +
                    me.at(6)*u.GetWeightL(Offset,6));

        TG+=Green::Me(green_offset+Offset)*
                (
                    me.at(0)*u.GetWeightR(Offset,0)+me.at(1)*u.GetWeightR(Offset,1)
                    +
                    me.at(2)*u.GetWeightR(Offset,2)+me.at(3)*u.GetWeightR(Offset,3)
                    +
                    me.at(4)*u.GetWeightR(Offset,4)+me.at(5)*u.GetWeightR(Offset,5)
                    +
                    me.at(6)*u.GetWeightR(Offset,6));
    }

    p.UpdateWeightL(GT/Green::Me(green_offset));
    p.UpdateWeightR(TG/Green::Me(green_offset));
}



int Procedure::ChooseL()
{
    double L=mt.randDblExc()*GT;

    int OffL=4;

    while ((L-= Green::Me(green_offset+OffL)*
            (
                me.at(0)*u.GetWeightL(OffL,0)+me.at(1)*u.GetWeightL(OffL,1)
                +
                me.at(2)*u.GetWeightL(OffL,2)+me.at(3)*u.GetWeightL(OffL,3)
                +
                me.at(4)*u.GetWeightL(OffL,4)+me.at(5)*u.GetWeightL(OffL,5)
                +
                me.at(6)*u.GetWeightL(OffL,6))) > 0) OffL-=1;

    int type=0;

    double LL=mt.randDblExc()*
            (
                me.at(0)*u.GetWeightL(OffL,0)+me.at(1)*u.GetWeightL(OffL,1)
                +
                me.at(2)*u.GetWeightL(OffL,2)+me.at(3)*u.GetWeightL(OffL,3)
                +
                me.at(4)*u.GetWeightL(OffL,4)+me.at(5)*u.GetWeightL(OffL,5)
                +
                me.at(6)*u.GetWeightL(OffL,6));

    while ((LL-= me.at(type)*u.GetWeightL(OffL,type)) > 0) type++;

    Uop2L=u.ChooseOperatorL(OffL,type,mt);

    return OffL;
}



int Procedure::ChooseR()
{
    double R=mt.randDblExc()*TG;

    int OffR=4;

    while ((R-= Green::Me(green_offset+OffR)*
            (
                me.at(0)*u.GetWeightR(OffR,0)+me.at(1)*u.GetWeightR(OffR,1)
                +
                me.at(2)*u.GetWeightR(OffR,2)+me.at(3)*u.GetWeightR(OffR,3)
                +
                me.at(4)*u.GetWeightR(OffR,4)+me.at(5)*u.GetWeightR(OffR,5)
                +
                me.at(6)*u.GetWeightR(OffR,6))) > 0) OffR-=1;


    int type=0;

    double RR=mt.randDblExc()*
            (
                me.at(0)*u.GetWeightR(OffR,0)+me.at(1)*u.GetWeightR(OffR,1)
                +
                me.at(2)*u.GetWeightR(OffR,2)+me.at(3)*u.GetWeightR(OffR,3)
                +
                me.at(4)*u.GetWeightR(OffR,4)+me.at(5)*u.GetWeightR(OffR,5)
                +
                me.at(6)*u.GetWeightR(OffR,6));

    while ((RR-= me.at(type)*u.GetWeightR(OffR,type)) > 0) type++;

    Uop2R=u.ChooseOperatorR(OffR,type,mt);

    return OffR;
}




void Procedure::DirectUpdate() {

    int dir;

    double deltae=0;

    green_offset=0;

    CircularTime deltatau;

    CircularTime defaultdeltatau;

    if (c.Length() >= 1) deltatau=c.Top(-1).ti-c.Top(1).ti;
    else deltatau=defaultdeltatau;

    UpdateWeights();

    p.EvaluateRunningPars(0,deltatau.Time(),c.Empty());

    do {
        if(p.ChooseDirectionRight(mt))
        {
            dir=1;
            do {

                if (p.ChooseCreateRight(mt)) {

                    green_offset+=ChooseR();

                    winding_number[((Uop2R/2)%2)]+=1.0*(1-(Uop2R%2)*2)/Site::L();

                    CircularTime dt(p.Time(mt,deltae,deltatau.Time()));

                    stL.ChangeState(Operators::Transformation(Operators::Conjugate(Uop2R)));

                    deltae+=Model::UpdateEnergy(Operators::GetOpNeighbors(Uop2R),stLO,stL);

                    stLO.ChangeState(Operators::Transformation(Operators::Conjugate(Uop2R)));

                    CircularTime newtime;

                    if (c.Empty()) newtime=CircularTime();
                    else newtime=c.Top(1).ti+dt;

                    c.Push(1,Kink{newtime,Uop2R});

                    u.UpdateL(Uop2R,stL,stR);
                    u.UpdateR(Uop2R,stL,stR);
                }
                else {

                    Kink kR=c.Top(dir);

                    winding_number[((kR.op/2)%2)]+=1.0*(1-(Operators::Conjugate(kR.op)%2)*2)/Site::L();

                    unsigned int intkop=Operators::Conjugate(kR.op);

                    Operator opkop=Operators::Transformation(intkop);

                    green_offset-=Green::UpdateOffset(opkop, stL,stR);

                    stR.ChangeState(opkop);

                    deltae-=Model::UpdateEnergy(Operators::GetOpNeighbors(intkop),stRO,stR);

                    stRO.ChangeState(opkop);

                    green_offset+=Green::UpdateOffset(opkop, stL,stR);

                    c.Pop(dir);

                    Green::Time(kR.ti);

                    u.UpdateL(kR.op,stL,stR);
                    u.UpdateR(kR.op,stL,stR);
                }

                if (c.Length() >= 1) deltatau=c.Top(-1).ti-c.Top(1).ti;
                else deltatau=defaultdeltatau;

                UpdateWeights();

                p.EvaluateRunningPars(deltae,deltatau.Time(),c.Empty());

            }   while (p.KeepGoingRight(mt));
        }
        else
        {
            dir=-1;
            do {

                if (p.ChooseCreateLeft(mt)) {

                    green_offset+=ChooseL();

                    winding_number[((Uop2L/2)%2)]+=1.0*(1-(Uop2L%2)*2)/Site::L();

                    CircularTime dt(p.Time(mt,deltae,deltatau.Time()));

                    stR.ChangeState(Operators::Transformation(Uop2L));

                    deltae-=Model::UpdateEnergy(Operators::GetOpNeighbors(Uop2L),stRO,stR);

                    stRO.ChangeState(Operators::Transformation(Uop2L));

                    CircularTime newtime;

                    if (c.Empty()) newtime=CircularTime();
                    else newtime=c.Top(1).ti+dt;

                    c.Push(-1,Kink{newtime,Uop2L});

                    u.UpdateL(Uop2L,stL,stR);
                    u.UpdateR(Uop2L,stL,stR);
                }
                else {

                    Kink kL=c.Top(dir);

                    winding_number[((kL.op/2)%2)]+=1.0*(1-(Operators::Conjugate(kL.op)%2)*2)/Site::L();

                    Operator opkop=Operators::Transformation(kL.op);

                    green_offset-=Green::UpdateOffset(opkop, stL,stR);

                    stL.ChangeState(opkop);

                    deltae+=Model::UpdateEnergy(Operators::GetOpNeighbors(kL.op),stLO,stL);

                    stLO.ChangeState(opkop);

                    green_offset+=Green::UpdateOffset(opkop, stL,stR);

                    c.Pop(dir);

                    Green::Time(kL.ti);

                    u.UpdateL(kL.op,stL,stR);
                    u.UpdateR(kL.op,stL,stR);
                }

                if (c.Length() >= 1) deltatau=c.Top(-1).ti-c.Top(1).ti;
                else deltatau=defaultdeltatau;

                UpdateWeights();

                p.EvaluateRunningPars(deltae,deltatau.Time(),c.Empty());

            } while (p.KeepGoingLeft(mt));
        }

    } while (green_offset!=0);

    if (c.Length() >= 1) deltatau=c.Top(-1).ti-c.Top(1).ti;
    else deltatau=defaultdeltatau;

    UpdateWeights();

    p.EvaluateRunningPars(0,deltatau.Time(),c.Empty());

//    int index=0;
//    bool isequal=true;

//    for (unsigned int i=0; i < stL.GetState().size(); i++)
//    {
//        if (stL.GetState().at(i).first!=stR.GetState().at(i).first)
//        {
//            index=i;
//            isequal=false;
//        }
//        if (stL.GetState().at(i).second!=stR.GetState().at(i).second)
//        {
//            index=i;
//            isequal=false;
//        }
//    }

//    double off=Green::Offset(stL,stR);

//    assert ((isequal) && ("StateL is not equal to StateR"));

//    assert ((GT==TG) && ("Weights is not equal"));

    CircularTime greentime(p.Time(mt,0,deltatau.Time()));

    if (c.Empty()) Green::Time(CircularTime());
    else Green::Time(c.Top(1).ti+greentime);
}
