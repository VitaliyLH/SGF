#include "Model.h"
#include "Procedure.h"
#include "Measurement.h"
#include "rapidjson/document.h"

#include <QTextStream>
#include <QFile>

#include "OffSetOperator.h"

#include <iomanip>
#include <iostream>

double Model::delta=0;
double Model::v=3;
double Model::j=0;
double Model::tp=0.1;
double Model::tpn=0.1;
double Model::tn=0.1;
double Model::t2=1;

int Green::cutoff=4;
CircularTime Green::ti=CircularTime();

int Site::l=8;
vector<vector<int>> Site::si;

vector<Operator> Operators::operators;
vector<vector<int>> Operators::int_neighbors;
vector<vector<Operator>> Operators::op_neighbors;

static string ReadInputFile( const char* path )
{
    FILE* file = fopen( path, "rb" );
    if ( !file )
        return string("");
    fseek( file, 0, SEEK_END );
    long size = ftell( file );
    fseek( file, 0, SEEK_SET );
    string text;
    char* buffer = new char[size + 1];
    buffer[size] = 0;
    if ( fread( buffer, 1, size, file ) == (unsigned long)size )
        text = buffer;
    fclose( file );
    delete[] buffer;
    return text;
}



QString GetFileName(QString number, QString name, QString extension1, QString extension2)
{
    number.remove(extension1);
    number.append(extension2);

    int index = number.indexOf(extension2);

    QString FileName(number);

    FileName.insert(index,name);

    return FileName;
}



int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        puts("Input File (*.json) not found!");
        return 1;
    }

    string input=ReadInputFile(argv[2]);

    string json=ReadInputFile(argv[1]);

    rapidjson::Document d;
    d.Parse<0>(json.c_str());

    assert(d["Parameters"].IsObject());
    const rapidjson::Value& parameters = d["Parameters"];

    assert(d["SGF"].IsObject());
    const rapidjson::Value& sgf = d["SGF"];

    unsigned int L=parameters["L"].GetUint();
    double T=parameters["T"].GetDouble();
    double Error=sgf["Error"].GetDouble();
    unsigned int NumBins=sgf["Bins"].GetUint();
    unsigned int ThermalizationSteps=sgf["WarmSteps"].GetUint();
    unsigned int MinMeasureSteps=sgf["MinMeasSteps"].GetUint();
    unsigned int InMeasureSteps=sgf["InMeasSteps"].GetUint();

    if (input == "")
    {
        Procedure p(parameters);

        for (unsigned int Steps=1; Steps <= ThermalizationSteps; Steps++)
        {
            p.DirectUpdate();
        }

        unsigned int Times=1;

        unsigned int BinsLenght=MinMeasureSteps/NumBins;

        unsigned int MinBinsLenght=InMeasureSteps/NumBins;

        Measurement me(L,T,BinsLenght,NumBins,Error);

        for (unsigned int Steps=1; Steps <= MinMeasureSteps; Steps++)
        {
            p.DirectUpdate();
            me.Make(p);
        }

        //            me.Coefficient(Times);

        //            do
        //            {
        //                BinsLenght=Times*MinBinsLenght;

        //                me.SetBinsLenght(BinsLenght);

        //                for (unsigned int Steps=1; Steps <= Times*InMeasureSteps; Steps++)
        //                {
        //                    p.DirectUpdate();
        //                    me.Make(p);
        //                }

        //                me.Coefficient(Times);

        //            } while (!me.isDone());

        QFile FileConfiguration(GetFileName(argv[1],"Configuration",".json",".yaml"));
        FileConfiguration.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileSpin(GetFileName(argv[1],"Spin",".json",".dat"));
        FileSpin.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FilePseudo(GetFileName(argv[1],"Pseudo",".json",".dat"));
        FilePseudo.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileP0(GetFileName(argv[1],"P0",".json",".dat"));
        FileP0.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileOutput(GetFileName(argv[1],"Output",".json",".dat"));
        FileOutput.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileLocalSpin(GetFileName(argv[1],"LocalSpin",".json",".dat"));
        FileLocalSpin.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileLocalPseudo(GetFileName(argv[1],"LocalPseudo",".json",".dat"));
        FileLocalPseudo.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileStiffness(GetFileName(argv[1],"Stiffness",".json",".dat"));
        FileStiffness.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileStructureFactorAFM(GetFileName(argv[1],"StructureFactorAFM",".json",".dat"));
        FileStructureFactorAFM.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileStructureFactorCO(GetFileName(argv[1],"StructureFactorCO",".json",".dat"));
        FileStructureFactorCO.open(QIODevice::WriteOnly | QIODevice::Text);

        p.WriteConfiguration(FileConfiguration);

        me.PrintOutput(FileOutput);
        me.Printsz(FileSpin);
        me.PrintSz(FilePseudo);
        me.PrintP0(FileP0);
        me.PrintLocalsz(FileLocalSpin);
        me.PrintLocalSz(FileLocalPseudo);
        me.PrintStiffness(FileStiffness);
        me.PrintStructureFactorsz(FileStructureFactorAFM);
        me.PrintStructureFactorSz(FileStructureFactorCO);

        FileSpin.close();
        FilePseudo.close();
        FileP0.close();
        FileOutput.close();
        FileLocalSpin.close();
        FileLocalPseudo.close();
        FileStiffness.close();
        FileStructureFactorAFM.close();
        FileStructureFactorCO.close();
    }
    else
    {
        rapidjson::Document c;
        c.Parse<0>(input.c_str());
        assert(c["configuration"].IsObject());
        const rapidjson::Value& configuration = c["configuration"];

        Procedure p(parameters,configuration);

        for (unsigned int Steps=1; Steps <= ThermalizationSteps; Steps++)
        {
            p.DirectUpdate();
        }

        unsigned int Times=1;

        unsigned int BinsLenght=MinMeasureSteps/NumBins;

        unsigned int MinBinsLenght=InMeasureSteps/NumBins;

        Measurement me(L,T,BinsLenght,NumBins,Error);

        for (unsigned int Steps=1; Steps <= MinMeasureSteps; Steps++)
        {
            p.DirectUpdate();
            me.Make(p);
        }

        //            me.Coefficient(Times);

        //            do
        //            {
        //                BinsLenght=Times*MinBinsLenght;

        //                me.SetBinsLenght(BinsLenght);

        //                for (unsigned int Steps=1; Steps <= Times*InMeasureSteps; Steps++)
        //                {
        //                    p.DirectUpdate();
        //                    me.Make(p);
        //                }

        //                me.Coefficient(Times);

        //            } while (!me.isDone());

        QFile FileConfiguration(GetFileName(argv[1],"Configuration",".json",".yaml"));
        FileConfiguration.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileSpin(GetFileName(argv[1],"Spin",".json",".dat"));
        FileSpin.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FilePseudo(GetFileName(argv[1],"Pseudo",".json",".dat"));
        FilePseudo.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileP0(GetFileName(argv[1],"P0",".json",".dat"));
        FileP0.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileOutput(GetFileName(argv[1],"Output",".json",".dat"));
        FileOutput.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileLocalSpin(GetFileName(argv[1],"LocalSpin",".json",".dat"));
        FileLocalSpin.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileLocalPseudo(GetFileName(argv[1],"LocalPseudo",".json",".dat"));
        FileLocalPseudo.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileStiffness(GetFileName(argv[1],"Stiffness",".json",".dat"));
        FileStiffness.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileStructureFactorAFM(GetFileName(argv[1],"StructureFactorAFM",".json",".dat"));
        FileStructureFactorAFM.open(QIODevice::WriteOnly | QIODevice::Text);

        QFile FileStructureFactorCO(GetFileName(argv[1],"StructureFactorCO",".json",".dat"));
        FileStructureFactorCO.open(QIODevice::WriteOnly | QIODevice::Text);

        p.WriteConfiguration(FileConfiguration);

        me.PrintOutput(FileOutput);
        me.Printsz(FileSpin);
        me.PrintSz(FilePseudo);
        me.PrintP0(FileP0);
        me.PrintLocalsz(FileLocalSpin);
        me.PrintLocalSz(FileLocalPseudo);
        me.PrintStiffness(FileStiffness);
        me.PrintStructureFactorsz(FileStructureFactorAFM);
        me.PrintStructureFactorSz(FileStructureFactorCO);

        FileSpin.close();
        FilePseudo.close();
        FileP0.close();
        FileOutput.close();
        FileLocalSpin.close();
        FileLocalPseudo.close();
        FileStiffness.close();
        FileStructureFactorAFM.close();
        FileStructureFactorCO.close();
    }

    return 0;
}
