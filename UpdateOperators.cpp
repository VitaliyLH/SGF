#include "UpdateOperators.h"

UpdateOperators::UpdateOperators()
{

}



UpdateOperators::~UpdateOperators()
{

}



UpdateOperators::UpdateOperators(const State &stL, const State &stR)
{   
    size_one=Operators::GetOperators().size()/7;
    for (unsigned int i=0; i < 7*size_one; i++) opsL.push_back(0);
    for (unsigned int i=0; i < 7*size_one; i++) opsR.push_back(0);
    ActiveHoppings(stL,stR);
}


map<int,OffSetOperator> UpdateOperators::Initialize(OffSetOperator ofp, OffSetOperator ofptype, int offset)
{
    pair<int,OffSetOperator> p;
    map<int,OffSetOperator> mp;

    for (int i=-4; i <= 4; i++)
    {
        p=(i==offset ? make_pair(i,ofptype) : make_pair(i,ofp));
        mp.insert(p);
    }

    return mp;
}



void UpdateOperators::ActiveHoppings(State stL, State stR)
{

    OffSetOperator ofpL(size_one);
    OffSetOperator ofpR(size_one);

    OffSetOperator ofpL0(size_one);
    OffSetOperator ofpL1(size_one);
    OffSetOperator ofpL2(size_one);
    OffSetOperator ofpL3(size_one);
    OffSetOperator ofpL4(size_one);
    OffSetOperator ofpL5(size_one);
    OffSetOperator ofpL6(size_one);

    OffSetOperator ofpR0(size_one);
    OffSetOperator ofpR1(size_one);
    OffSetOperator ofpR2(size_one);
    OffSetOperator ofpR3(size_one);
    OffSetOperator ofpR4(size_one);
    OffSetOperator ofpR5(size_one);
    OffSetOperator ofpR6(size_one);

    for (auto i=0; i < Operators::GetOperators().size(); i++)
    {
        auto op=Operators::GetOperators().at(i);

        int meL=op.MatrixElement(stR.SiteState(op.GetSsite()),stR.SiteState(op.GetTsite()));

        opsL[Operators::Transformation(op)]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
        opsR[Operators::Transformation(op)]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));

        if (meL!=0)
        {
            int ops=Operators::Transformation(op)-op.GetType()*size_one;

            if (op.GetType()==0)
            {
                ofpL0[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpL0.insert_in_subset(ops);
            }

            if (op.GetType()==1)
            {
                ofpL1[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpL1.insert_in_subset(ops);
            }

            if (op.GetType()==2)
            {
                ofpL2[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpL2.insert_in_subset(ops);
            }

            if (op.GetType()==3)
            {
                ofpL3[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpL3.insert_in_subset(ops);
            }

            if (op.GetType()==4)
            {
                ofpL4[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpL4.insert_in_subset(ops);
            }

            if (op.GetType()==5)
            {
                ofpL5[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpL5.insert_in_subset(ops);
            }

            if (op.GetType()==6)
            {
                ofpL6[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpL6.insert_in_subset(ops);
            }
        }

        int meR=op.MatrixElement(stL.SiteState(op.GetSsite()),stL.SiteState(op.GetTsite()));

        if (meR!=0)
        {
            int ops=Operators::Conjugate(Operators::Transformation(op)-op.GetType()*size_one);

            if (op.GetType()==0)
            {
                ofpR0[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpR0.insert_in_subset(ops);
            }

            if (op.GetType()==1)
            {
                ofpR1[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpR1.insert_in_subset(ops);
            }

            if (op.GetType()==2)
            {
                ofpR2[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpR2.insert_in_subset(ops);
            }

            if (op.GetType()==3)
            {
                ofpR3[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpR3.insert_in_subset(ops);
            }

            if (op.GetType()==4)
            {
                ofpR4[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpR4.insert_in_subset(ops);
            }

            if (op.GetType()==5)
            {
                ofpR5[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpR5.insert_in_subset(ops);
            }

            if (op.GetType()==6)
            {
                ofpR6[ops]=2*(op.GetPOffset()+abs(op.GetSOffsetS()));
                ofpR6.insert_in_subset(ops);
            }
        }
    }

    l.push_back(Initialize(ofpL,ofpL0,3));
    l.push_back(Initialize(ofpL,ofpL1,3));
    l.push_back(Initialize(ofpL,ofpL2,3));
    l.push_back(Initialize(ofpL,ofpL3,3));
    l.push_back(Initialize(ofpL,ofpL4,3));
    l.push_back(Initialize(ofpL,ofpL5,3));
    l.push_back(Initialize(ofpL,ofpL6,4));

    r.push_back(Initialize(ofpR,ofpR0,3));
    r.push_back(Initialize(ofpR,ofpR1,3));
    r.push_back(Initialize(ofpR,ofpR2,3));
    r.push_back(Initialize(ofpR,ofpR3,3));
    r.push_back(Initialize(ofpR,ofpR4,3));
    r.push_back(Initialize(ofpR,ofpR5,3));
    r.push_back(Initialize(ofpR,ofpR6,4));
}



int UpdateOperators::ChooseOperatorL(int offset, int type, MTRand &mt)
{
    int sh=type*size_one;
    int op=l.at(type).at(offset).operator_at(mt.randInt(l.at(type).at(offset).subset_size()-1));
    return op+sh;
}



int UpdateOperators::ChooseOperatorR(int offset, int type, MTRand &mt)
{
    int sh=type*size_one;
    int op=r.at(type).at(offset).operator_at(mt.randInt(r.at(type).at(offset).subset_size()-1));
    return op+sh;
}



double UpdateOperators::GetWeightL(int Offset, int type)
{
    return 1.0*(l.at(type).at(Offset).subset_size());
}



double UpdateOperators::GetWeightR(int Offset, int type)
{
    return 1.0*(r.at(type).at(Offset).subset_size());
}



void UpdateOperators::RemoveL(int oper)
{
    Operator oops=Operators::Transformation(oper);
    int off=opsL[oper];

    l.at(oops.GetType()).at(off).erase_from_subset(oper-oops.GetType()*size_one);
}



void UpdateOperators::UpdateL(int oper, State& stL, State& stR)
{
    vector<int> op=Operators::GetIntNeighbors(oper);

    for (unsigned int i=0;i<op.size();i++) RemoveL(op.at(i));

    for (unsigned int i=0;i<op.size();i++)
    {
        Operator oops=Operators::Transformation(op.at(i));

        int me=oops.MatrixElement(stR.SiteState(oops.GetSsite()),stR.SiteState(oops.GetTsite()));

        int toff=Green::TermOffset(oops,stL,stR);

        opsL[op.at(i)]=toff;

        if (me!=0) l.at(oops.GetType()).at(toff).insert_in_subset(op.at(i)-oops.GetType()*size_one);
    }
}



void UpdateOperators::RemoveR(int oper)
{
    Operator oops=Operators::Transformation(oper);
    int off=opsR[Operators::Conjugate(oper)];

    r.at(oops.GetType()).at(off).erase_from_subset(Operators::Conjugate(oper-oops.GetType()*size_one));
}



void UpdateOperators::UpdateR(int oper, State& stL, State& stR)
{
    vector<int> op=Operators::GetIntNeighbors(oper);

    for (unsigned int i=0;i<op.size();i++) RemoveR(op.at(i));

    for (unsigned int i=0;i<op.size();i++)
    {
        Operator oops=Operators::Transformation(op.at(i));

        int me=oops.MatrixElement(stL.SiteState(oops.GetSsite()),stL.SiteState(oops.GetTsite()));

        int toff=Green::TermOffset(Operators::Transformation(Operators::Conjugate(op.at(i))),stL,stR);

        opsR[Operators::Conjugate(op.at(i))]=toff;

        if (me!=0) r.at(oops.GetType()).at(toff).insert_in_subset(Operators::Conjugate(op.at(i)-oops.GetType()*size_one));
    }
}
