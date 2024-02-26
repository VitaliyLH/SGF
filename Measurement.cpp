#include "Measurement.h"

Measurement::Measurement()
{

}



Measurement::~Measurement()
{

}



Measurement::Measurement(unsigned int L, double T, unsigned int BinsLenght, unsigned int NumberBins, double MaxError) :
    norm(0),
    max_error(MaxError),
    count(0),
    temperature(T),
    count_bad(0),
    bin_number(0),
    bins_lenght(BinsLenght),
    number_sites(L*L),
    coefficient_potential(0),
    coefficient_kinetic(0),
    auto_mean(11,AutoCorrelation(NumberBins)),
    local_sz(number_sites,AutoCorrelation(NumberBins)),
    local_Sz(number_sites,AutoCorrelation(NumberBins))
{
    for (unsigned int i=0; i < NumberBins; i++) normalization.push_back(0);
}



void Measurement::PrintOutput(QFile &File)
{
    QTextStream stream( &File );

    stream.setFieldWidth(5);
    stream.setRealNumberPrecision(10);
    stream.setFieldAlignment(QTextStream::AlignLeft);

    stream << "\n";
    stream << temperature << "\t";

    for (unsigned int i=0; i <= 1; i++)
    {
        auto_mean.at(i).Print(stream);
        stream << "\t";
    }
}



void Measurement::Printsz(QFile &File)
{
    QTextStream stream( &File );

    stream.setFieldWidth(5);
    stream.setRealNumberPrecision(10);
    stream.setFieldAlignment(QTextStream::AlignLeft);

    stream << "\n";
    stream << temperature << "\t";

    for (unsigned int i=2; i <= 3; i++)
    {
        auto_mean.at(i).Print(stream);
        stream << "\t";
    }
}



void Measurement::PrintSz(QFile &File)
{
    QTextStream stream( &File );

    stream.setFieldWidth(5);
    stream.setRealNumberPrecision(10);
    stream.setFieldAlignment(QTextStream::AlignLeft);

    stream << "\n";
    stream << temperature << "\t";

    for (unsigned int i=4; i <= 5; i++)
    {
        auto_mean.at(i).Print(stream);
        stream << "\t";
    }
}



void Measurement::PrintLocalsz(QFile &File)
{
    QTextStream stream( &File );

    stream.setFieldWidth(5);
    stream.setRealNumberPrecision(10);
    stream.setFieldAlignment(QTextStream::AlignLeft);

    stream << "\n";
    stream << sqrt(local_sz.size()) << "\n";
    for (unsigned int i=0; i < local_sz.size(); i++)
    {
        stream << i;
        local_sz.at(i).Print(stream);
        stream << "\n";
    }

    stream << 1024;
}



void Measurement::PrintLocalSz(QFile &File)
{
    QTextStream stream( &File );

    stream.setFieldWidth(5);
    stream.setRealNumberPrecision(10);
    stream.setFieldAlignment(QTextStream::AlignLeft);

    stream << "\n";
    stream << sqrt(local_Sz.size()) << "\n";
    for (unsigned int i=0; i < local_Sz.size(); i++)
    {
        stream << i;
        local_Sz.at(i).Print(stream);
        stream << "\n";
    }

    stream << 1024;
}



void Measurement::PrintP0(QFile &File)
{
    QTextStream stream( &File );

    stream.setFieldWidth(5);
    stream.setRealNumberPrecision(10);
    stream.setFieldAlignment(QTextStream::AlignLeft);

    stream << "\n";
    stream << temperature;
    for (unsigned int i=6; i <= 7; i++)
    {
        auto_mean.at(i).Print(stream);
        stream << "\t";
    }
}



void Measurement::PrintStructureFactorsz(QFile &File)
{
    QTextStream stream( &File );

    stream.setFieldWidth(5);
    stream.setRealNumberPrecision(10);
    stream.setFieldAlignment(QTextStream::AlignLeft);

    stream << "\n";
    stream << temperature;
    auto_mean.at(8).Print(stream);
}


void Measurement::PrintStructureFactorSz(QFile &File)
{
    QTextStream stream( &File );

    stream.setFieldWidth(5);
    stream.setRealNumberPrecision(10);
    stream.setFieldAlignment(QTextStream::AlignLeft);

    stream << "\n";
    stream << temperature;
    auto_mean.at(9).Print(stream);
}



void Measurement::PrintStiffness(QFile &File)
{
    QTextStream stream( &File );

    stream.setFieldWidth(5);
    stream.setRealNumberPrecision(10);
    stream.setFieldAlignment(QTextStream::AlignLeft);

    stream << "\n";
    stream << temperature;
    auto_mean.at(10).Print(stream);
}



//void Measurement::GetStates(Procedure& p)
//{
//    state.clear();
//    deltatau.clear();

//    State stL=p.GetState();

//    if (p.GetKinks().Length()==0) {
//        state.push_back(stL);
//        deltatau.push_back(1.0);
//    }
//    else {

//        double sumdt=0.0;

//        bool issigned=false;

//        double dt=(p.GetKinks().Top(-1).Time()-p.GetKinks().Top(1).Time()).Time();

//        state.push_back(stL);
//        deltatau.push_back(dt);

//        sumdt+=dt;

//        if (p.GetKinks().Length()!=0) {
//            std::deque<Kink> k=p.GetKinks().GetQue();
//            for (unsigned int i=0;i<p.GetKinks().Length()-1;i++) {
//                if (k.at(i+1).Time().Time() < k.at(i).Time().Time() && !issigned) {
//                    issigned=true;
//                    dt=1-k.at(i).Time().Time()+k.at(i+1).Time().Time();
//                    sumdt+=dt;
//                    stL.ChangeState(p.GetOperators().Transformation(k.at(i).KinkOp()));
//                    state.push_back(stL);
//                    deltatau.push_back(dt);
//                } else {
//                    dt=fabs(k.at(i+1).Time().Time()-k.at(i).Time().Time());
//                    sumdt+=dt;
//                    stL.ChangeState(p.GetOperators().Transformation(k.at(i).KinkOp()));
//                    state.push_back(stL);
//                    deltatau.push_back(dt);
//                }
//            }
//            stL.ChangeState(p.GetOperators().Transformation(k.at(p.GetKinks().Length()-1).KinkOp()));
//        }
//    }
//}



bool Measurement::isDone()
{
    return ( ( isDonePotential() && isDoneKinetic() ) || count_bad > 10 );
}



bool Measurement::isDonePotential()
{
    return (
                fabs(auto_mean.at(0).GetAverage()) < 1.0
                ?
                    auto_mean.at(0).GetError() < max_error
                  :
                    auto_mean.at(0).GetError() < max_error*fabs(auto_mean.at(0).GetAverage())
                    );
}



bool Measurement::isDoneKinetic()
{
    return (
                fabs(auto_mean.at(1).GetAverage()) < 1.0
                ?
                    auto_mean.at(1).GetError() < max_error
                  :
                    auto_mean.at(1).GetError() < max_error*fabs(auto_mean.at(1).GetAverage())
                    );
}



void Measurement::Coefficient(unsigned int& Times)
{
    unsigned int prevcp=coefficient_potential;
    unsigned int prevck=coefficient_kinetic;

    unsigned int cp=
            (
                fabs(auto_mean.at(0).GetAverage()) < 1.0
                ?
                    auto_mean.at(0).GetError()/max_error
                  :
                    auto_mean.at(0).GetError()/(max_error*(fabs(auto_mean.at(0).GetAverage())))
                    );

    unsigned int ck=
            (
                fabs(auto_mean.at(1).GetAverage()) < 1.0
                ?
                    auto_mean.at(1).GetError()/max_error
                  :
                    auto_mean.at(1).GetError()/(max_error*(fabs(auto_mean.at(1).GetAverage())))
                    );

    unsigned int c=cp+ck;

    if (c > 2) Times=c;
    else Times=1;

    unsigned int nextcp=cp;
    unsigned int nextck=ck;

    if ( (nextcp >= prevcp) && (nextck >= prevck) && !isDonePotential() && !isDoneKinetic()) count_bad++;

    coefficient_potential=cp;
    coefficient_kinetic=ck;

    bin_number=0;
}



vector<pair<int,int>> Measurement::GetTermDistribution(Procedure& p, unsigned int size)
{
    vector<int> TermDistribution(size,0);
    vector<pair<int,int> > TermDistributionXY;

    for (unsigned int i = 0;i < p.GetKinks().GetQue().size();i++)
        TermDistribution[p.GetKinks().GetQue().at(i).op] += 1;

    for (unsigned int i = 0;i < TermDistribution.size();i+=4)
    {
        TermDistributionXY.push_back(make_pair(TermDistribution[i] , TermDistribution[i+2]));
        TermDistributionXY.push_back(make_pair(TermDistribution[i+1] , TermDistribution[i+3]));
    }

    return TermDistributionXY;
}



pair<double,double> Measurement::GetWindingNumber(Procedure& p, unsigned int size)
{
    pair<double,double> WindingNumber(0,0);
    vector<pair<int,int> > TermDistribution = GetTermDistribution(p,size);

    for (unsigned int i = 0;i < TermDistribution.size();i+=2)
    {
        WindingNumber.first += ( TermDistribution[i].first - TermDistribution[i+1].first );
        WindingNumber.second += ( TermDistribution[i].second - TermDistribution[i+1].second );
    }

    WindingNumber.first /=  Site::L();

    WindingNumber.second /= Site::L();

    return WindingNumber;
}



void Measurement::Make(Procedure& p)
{
    //    GetStates(p);

    double n=1/p.GetProbabilities().GetRenormalization();

    norm+=n;

    double szA=0.0;
    double szB=0.0;

    double SzA=0.0;
    double SzB=0.0;

    double P0A=0.0;
    double P0B=0.0;

    for (int y=0;y < Site::L();y++)
    {
        for (int x=0;x < Site::L();x++)
        {
            if ((x+y)%2==0) {
                SzA+=p.GetState().SiteState(Site::GetSite(x,y)).first;
                P0A+=(1-p.GetState().SiteState(Site::GetSite(x,y)).first*p.GetState().SiteState(Site::GetSite(x,y)).first);
                szA+=p.GetState().SiteState(Site::GetSite(x,y)).second;
            }
            if ((x+y)%2==1) {
                SzB+=p.GetState().SiteState(Site::GetSite(x,y)).first;
                P0B+=(1-p.GetState().SiteState(Site::GetSite(x,y)).first*p.GetState().SiteState(Site::GetSite(x,y)).first);
                szB+=p.GetState().SiteState(Site::GetSite(x,y)).second;
            }
        }
    }

    szA/=(0.5*Site::L()*Site::L());
    szB/=(0.5*Site::L()*Site::L());

    SzA/=(0.5*Site::L()*Site::L());
    SzB/=(0.5*Site::L()*Site::L());

    P0A/=(0.5*Site::L()*Site::L());
    P0B/=(0.5*Site::L()*Site::L());

    double stfactSz=(SzA-SzB)*(SzA-SzB)/4.0;

    double stfactsz=(szA-szB)*(szA-szB)/4.0;

    double windingnumber2=0.5*(p.GetWindingNumber().at(0)*p.GetWindingNumber().at(0)+p.GetWindingNumber().at(1)*p.GetWindingNumber().at(1))/p.GetProbabilities().GetBeta();

//    unsigned int Size=p.GetOperators().GetOperators().size();

//    pair<double,double> windingnumber=GetWindingNumber(p,Size);

//    double windingnumber2test=0.5*(windingnumber.first*windingnumber.first+windingnumber.second*windingnumber.second)/p.GetProbabilities().GetBeta();

//    assert ((windingnumber2==windingnumber2test) && ("Windings is not equal"));

    double potenergy=Model::Energy(p.GetState());

    //    double potenergy=0.0;

    //    for (unsigned int i=0;i < state.size();i++) { potenergy+=p.GetModel().Energy(state.at(i),p.GetSite())*deltatau.at(i); }

    double potentialenergy=potenergy/(Site::L()*Site::L());

    double kineticenergy=(-1.0*p.GetKinks().Length())/(p.GetProbabilities().GetBeta()*Site::L()*Site::L());

    for (unsigned int i=0; i < local_sz.size(); i++) local_sz.at(i).Add(bin_number,p.GetState().SiteState(i).second*n);
    for (unsigned int i=0; i < local_Sz.size(); i++) local_Sz.at(i).Add(bin_number,p.GetState().SiteState(i).first*n);

    auto_mean.at(0).Add(bin_number,potentialenergy*n);
    auto_mean.at(1).Add(bin_number,kineticenergy*n);
    auto_mean.at(2).Add(bin_number,szA*n);
    auto_mean.at(3).Add(bin_number,szB*n);
    auto_mean.at(4).Add(bin_number,SzA*n);
    auto_mean.at(5).Add(bin_number,SzB*n);
    auto_mean.at(6).Add(bin_number,P0A*n);
    auto_mean.at(7).Add(bin_number,P0B*n);
    auto_mean.at(8).Add(bin_number,stfactsz*n);
    auto_mean.at(9).Add(bin_number,stfactSz*n);
    auto_mean.at(10).Add(bin_number,0.5*windingnumber2*n);

    count++;

    if (count==bins_lenght)
    {
        normalization[bin_number]+=norm;
        for (unsigned int i=0; i < auto_mean.size(); i++) auto_mean.at(i).Count(bin_number,normalization.at(bin_number));
        for (unsigned int i=0; i < local_sz.size(); i++) local_sz.at(i).Count(bin_number,normalization.at(bin_number));
        for (unsigned int i=0; i < local_Sz.size(); i++) local_Sz.at(i).Count(bin_number,normalization.at(bin_number));
        bin_number++;
        count=0;
        norm=0;
    }
}
