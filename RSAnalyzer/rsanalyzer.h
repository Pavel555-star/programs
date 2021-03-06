#ifndef RSANALYZER_H
#define RSANALYZER_H
#include "QString"
#include "QByteArray"
#include "QFile"
#include "QVector"
#include "QImage"

class RSAnalyzer
{
public:
    RSAnalyzer();

    int R = 1;  // Proměnná určující důležitost červené barvy pro fraktální analýzu
    int G = 1;  // Proměnná určující důležitost zelené barvy pro fraktální analýzu
    int B = 1;  // Proměnná určující důležitost modré barvy pro fraktální analýzu
    int CountBytes = 2; // Proměnná určující počet bytů  čísla v binárních souborech

    bool Load_CSV(QString inFileName, int Count_Deliminator); // Funkce pro načtení CSV souborů
    bool Load_Binary(QString InputFile); // Funkce R/S načtení binárních souborů (audio, ...)
    bool Load_Picture(QString InputFile); // Funkce R/S pro načtení obrázků
    bool Save_File(QString outFileName, int format_x, int format_y); // Funkce pro ukládání souborů ve formátu CSV

    void Analyze_Curve(int IterMin, int iterMax); // Funkce pro R/S analýzu CSV souborů a binárních souborů (audio, ...)
    void Analyze_Picture(int IterMin, int iterMax); // Funkce R/S pro analýzu obrázků

    QVector<double> Input; // Kontejnér vstupních hodnot
    QVector<double> Input_2; // Kontejnér vstupních hodnot 2
    QVector<double> Input_3; // Kontejnér vstupních hodnot 3

    QVector<double> Output_x; // Kontejnér výstupních hodnot logaritmů počtu půlení n
    QVector<double> Output_RS; // Kontejnér výstupních hodnot R/S
    QVector<double> Output_2; // Kontejnér výstupních hodnot R/S



private:

    int Dimension_Curve = 1; // Proměnná ve které je uložena dimenze křivky

    double *pInput; // Ukazatel na Input
    double *pInput_2; // Ukazatel na Input_2
    double *pInput_3; // Ukazatel na Input_3
    double *pOutput_x; // Ukazatel na Output_x
    double *pOutput_RS;// Ukazatel na Output_RS
    double *pOutput_2;// Ukazatel na Output_2

    QString inFileName;
    QImage Image;
};

#endif // RSANALYZER_H
