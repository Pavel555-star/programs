#include "mainwindow.h"
#include "rsanalyzer.h"

#include "string.h"
using namespace std;

#include <QApplication>
#include <QTranslator>


int main(int argc, char *argv[])
{


    if (argc < 4)
    {
        QApplication a(argc, argv); // Okno grafického rozhraní

        QTranslator translator;
        translator.load("localization_non_czech");
        a.installTranslator(&translator);
        MainWindow w;
        w.show();
        return a.exec();
    }
    else {
        {
            RSAnalyzer analyzer; // Ovládání aplikace parametry příkazové řádky
            int IterMin = 0;
            int IterMax = 0;
            int CountBytes = 2;
            int ColumnNumber = 1;


            string color;
            string bytes;
            string choice;
            string min = "";
            string max = "";
            string column = "";
            if (argc > 5)
                {
                min = argv[4];
                IterMin = stoi(min);
                max = argv[5];
                IterMax = stoi(max);
                }



            choice = argv[3];
            if (choice == "t") // Pro textové soubory
                {
                if (argc > 6)
                    {
                    column = argv[6];
                    ColumnNumber = stoi(column); // Číslo sloupce SCV souboru (1 - první sloupec)
                    }
                analyzer.Load_CSV(argv[1], ColumnNumber); // Načtení vstupního souboru
                analyzer.Analyze_Curve(IterMin, IterMax);
                }
            if (choice == "b") // Pro binární soubory
                {
                if (argc > 6)
                    {
                    bytes = argv[6];
                    CountBytes = stoi(bytes); // Počty bytů čísel binárního souboru
                    analyzer.CountBytes = CountBytes;
                    }
                analyzer.Load_Binary(argv[1]); // Načtení vstupního souboru
                analyzer.Analyze_Curve(IterMin, IterMax);
                }
            if (choice == "p") // Pro obrázky
                {
                if (argc > 8)
                    {
                    color = argv[6]; // Zahrnutí jednotlivých barev do analýzy
                    analyzer.R = stoi(color);
                    color = argv[7];
                    analyzer.G = stoi(color);
                    color = argv[8];
                    analyzer.B = stoi(color);
                    }
                analyzer.Load_Picture(argv[1]); // Načtení vstupního souboru
                analyzer.Analyze_Picture(IterMin, IterMax);
                }

            analyzer.Save_File(argv[2],  2, 0); // Uložení výstupního souboru
        }
    }
}
// Autor programu: Ing. Pavel Florián Ph.D. dovoluje použití zdrojového kódu této aplikace, pokud je splněna podmínka uvedení autorova jména.
