#include "State.h"

State::State()
{

}

State::~State()
{

}


State::State(int L, int x)
{
    int nsites=L*L;

    double fill=-nsites;

    for (int i=0; i < nsites; i++) st.push_back(make_pair(-1.0,0));

    int i = 0;

    while (x!=fill) {
        st[i].first+=2.0;
        fill+=2.0;
        i++;
    }

//    for (int i=0; i < nsites; i++) st.push_back(make_pair(0,-0.5));

//    int i = 0;

//    while (x!=fill) {
//        st[i].second+=1.0;
//        fill+=2.0;
//        i++;
//    }
}



void State::ChangeState(Operator &op)
{
//    if (op.GetType()==0 || op.GetType()==2 || op.GetType()==4)
//    {
//        st[op.GetSsite()].first+=1;
//        st[op.GetSsite()].second+=0.5;

//        st[op.GetTsite()].first-=1;
//        st[op.GetTsite()].second-=0.5;
//    }

//    if (op.GetType()==1 || op.GetType()==3 || op.GetType()==5 )
//    {
//        st[op.GetSsite()].first+=1;
//        st[op.GetSsite()].second-=0.5;

//        st[op.GetTsite()].first-=1;
//        st[op.GetTsite()].second+=0.5;
//    }

//    if (op.GetType()==6)
//    {
//        st[op.GetSsite()].first+=2;
//        st[op.GetTsite()].first-=2;
//    }

    st[op.GetSsite()].first+=op.GetPOffset();
    st[op.GetSsite()].second+=op.GetSOffsetS();

    st[op.GetTsite()].first-=op.GetPOffset();
    st[op.GetTsite()].second-=op.GetSOffsetT();
}
