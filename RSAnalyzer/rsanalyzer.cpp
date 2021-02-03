#include "rsanalyzer.h"
#include <cmath>
using namespace std;

#include <QString>
#include <QByteArray>
#include <QFile>
#include <QVector>

RSAnalyzer::RSAnalyzer()
{
}
bool RSAnalyzer::Load_CSV(QString InputFile)
{
    // Načítací funkce programu
    inFileName = InputFile;
    QFile inFile(inFileName);
    QByteArray process_line;
    QByteArray znaky = "0123456789.";
    QByteArray Buffer_1;
    QByteArray Buffer_2;

    int lines = 0; // Počet řádků vstupního souboru
    int i, m;
    bool ok;

    if (!inFile.open(QIODevice::ReadOnly | QIODevice::Text)) // Načítací část programu
        return (false);

    while (!inFile.atEnd()){
        process_line = inFile.readLine();
        lines++;}
    inFile.close();

    Input.clear();
    Input.reserve(lines);
    pInput = Input.data();

    if (!inFile.open(QIODevice::ReadOnly | QIODevice::Text)) // Načítací část programu
        return (false);

    while ((!inFile.atEnd())) {
        process_line = inFile.readLine();
        m = process_line.size();
        for (i = 0; i < m; i++)
           if (znaky.indexOf(process_line[i]) != -1)
               Buffer_1.append(process_line[i]);
        Input.append(Buffer_1.toDouble(&ok));
        Buffer_1 = "";
        }
    inFile.close();
    return(true);
}

bool RSAnalyzer::Load_Wav(QString InputFile)
{
// Načítací funkce programu
int i;
long j;
inFileName = InputFile;
QFile inFile(inFileName);
QByteArray BA;

if (!inFile.open(QIODevice::ReadOnly)) // Načítací část programu
    return (false);

int size = inFile.size();
Input.clear();
Input.reserve((size/CountBytes) + 1);
BA.resize(size);
BA = inFile.readAll();

for (i = 0; i < (size); i++)
{
    if (i % CountBytes == 0)
    {
    j = BA[i];
    if (CountBytes > 1)
        j = j + (BA[i+1] << 8);
    if (CountBytes > 2)
        j = j +  (BA[i+2] << 16);
    if (CountBytes > 3)
        j = j + (BA[i+3] << 24);
    if (CountBytes > 4)
        j = j + (BA[i+4] << 32);
    if (CountBytes > 5)
        j = j + (BA[i+5] << 40);
    if (CountBytes > 6)
        j = j + (BA[i+6] << 48);
    if (CountBytes > 7)
        j = j + (BA[i+7] << 56);

    Input.append(double(j));
    }
}
inFile.close();
return(true);
}

bool RSAnalyzer::Load_Image(QString InputFile)
{
// Načítací funkce programu
Image.load(InputFile);
return(true);
}

bool RSAnalyzer::Save_File(QString outFileName, int format_x, int format_y)
{ // Ukládací funkce programu
    int c;
    int d = Output_x.size();
    QFile outFile(outFileName);
    QByteArray outString;
    QByteArray buffer1;
    QByteArray buffer2;
    QByteArray buffer3;
    if (!outFile.open(QIODevice::WriteOnly | QIODevice::Text))
      return (false);

pOutput_x = Output_x.data();
pOutput_RS = Output_RS.data();
pOutput_2 = Output_2.data();
outString.append("log2(N); log2(R/S); log2(box count) \n");

for (c = 0; c < d; c++)
{
    switch (format_x) {
    case 0:
        buffer1.setNum(pOutput_x[c]);
        break;
    case 1:
        buffer1.setNum((float)pOutput_x[c]);
        break;
    case 2:
        buffer1.setNum((int)pOutput_x[c]);
        break;}
    switch (format_y) {
    case 0:{
        buffer2.setNum(pOutput_RS[c]);
        buffer3.setNum(pOutput_2[c]);
        break;
    }
    case 1:{
        buffer2.setNum((float)pOutput_RS[c]);
        buffer3.setNum((float)pOutput_2[c]);

        break;}
    case 2:{
        buffer2.setNum((int)pOutput_RS[c]);
        buffer3.setNum((int)pOutput_2[c]);
        break;}}

    outString.append(buffer1);
    outString.append("; ");
    outString.append(buffer2);
    outString.append("; ");
    outString.append(buffer3);
    outString.append("\n");
    outFile.write(outString);
    outString ="";
}
outFile.close();
return(true);
}
void RSAnalyzer::Analyze_1D(int IterMin, int iterMax)
{ // Výpočetní část programu
    int scale;
    int i, j = 0;
    int o = 0;
    int n = Input.size();
    int npp = n;
    int iter = 1;
    int errors = 0;

    double range = 0;
    double s;
    double max_y = 0;
    double min_y = 1.7E+308;
    double prumer = 0;
    double RS = 0;

    double lenght;

    Output_x.clear();
    Output_RS.clear();
    Output_2.clear();
    Input.append(Input[n - 1]);

    for (iter =1; iter <= npp; iter = iter * 2)
        if (((iter > IterMin) and (iter < iterMax)) or iterMax == 0)
        {
        scale = n/iter;
            for (i = 0; i < iter; i++)
            {
                lenght = lenght + scale  + abs(Input[((i + 1) * scale)] - Input[i * scale]);
                for (j = 0; j < scale; j++)
                {
                     o = (i * scale) + j;
                     prumer = prumer + Input[o];
                   if (Input[o] > max_y)
                       max_y = Input[o];
                   if (Input[o] < min_y)
                        min_y = Input[o];
                 }
                 prumer = (prumer)/j;
                 for (j = 0; j < scale; j++)
                 {
                    o = (i * scale) + j;
                    s = s + ((prumer - Input[o])*(prumer - Input[o]));
                 }
                 range = max_y - min_y;
                 s = s/scale;
                 if (s > 0)
                     RS = RS + range/s;
                 else
                     errors++;

                 s = 0;
                 prumer = 0;
                 max_y = 0;
                 min_y = 1.7E+308;
             }
         RS = RS/(iter-errors);
         Output_x.append(log2(iter));
         Output_RS.append(log2(RS));
         Output_2.append(log2(lenght));
         RS = 0;
         lenght = 0;
         errors = 0;
         }
}
void RSAnalyzer::Analyze_2D(int IterMin, int iterMax)
{ // Výpočetní část programu

    int width;
    int height;
    int LineSize;
    QColor value;
    int sum_value; 

    int scale_x;
    int scale_y;
    int i, j, k, l = 0;
    int iter = 1;
    int errors = 0;

    long range = 0;
    long int s = 0;
    double s2 =0;
    long max = 0;
    long min = 2147483647;
    long prumer = 0;

    long p_prumer = 0;
    long count = 0;

    width = Image.width();
    height = Image.height();
    LineSize = Image.bytesPerLine();

    QVector<int> RGB;
    RGB.reserve(width*height*3);

    Output_x.clear();
    Output_RS.clear();
    Output_2.clear();

    int npp;
    if (width > height)
        npp = height;
    else
        npp = width;

    for (iter = 1; iter <= npp; iter = iter * 2)
        if (((iter > IterMin) and (iter < iterMax)) or iterMax == 0)
        {
             scale_x = width/iter;
             scale_y = height/iter;
             for (i = 0; i < iter; i++)
                for (j = 0; j < iter; j++)
                {
                    for (k = 0; k < scale_x; k++)
                        for (l = 0; l < scale_y; l++)
                        {
                            value = Image.pixelColor((i * scale_x) + k,(j * scale_y) + l);
                            sum_value = value.red() + value.green() + value.blue();
                            prumer = prumer + sum_value;
                            if (sum_value > max)
                                max = sum_value;
                            if (sum_value < min)
                                min = sum_value;
                         }
                    prumer = prumer/(scale_x * scale_y);

                    for (k = 0; k < scale_x; k++)
                        for (l = 0; l < scale_y; l++)
                        {
                            value = Image.pixelColor((i * scale_x) + k,(j * scale_y) + l);
                            sum_value = value.red() + value.green() + value.blue();
                            s = s + ((prumer - sum_value) * (prumer - sum_value));
                        }
                    if (p_prumer == 0)
                        p_prumer = prumer;
                    if (prumer >= p_prumer)
                        count++;

                range = max - min;
                s = s/(scale_x * scale_y);
                if (s > 0)
                    s2 = s2 + double(range)/double(s);
                else
                    errors++;
                prumer = 0;
                s = 0;
                max = 0;
                min = 0;
                }
        s2 = s2/((iter*iter)-errors);
        Output_x.append(double(log2(iter*iter)));
        Output_RS.append(log2(s2));
        Output_2.append(log2(count));
        s2 = 0;
        errors = 0;
        count = 0;
    }
}
