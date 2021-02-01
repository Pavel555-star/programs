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

            string color;
            string choice;
            string min = "";
            string max = "";
            if (argc > 5)
                {
                min = argv[4];
                IterMin = stoi(min);
                max = argv[5];
                IterMax = stoi(max);
                }
            if (argc > 8)
                {
                color = argv[6];
                analyzer.R = stoi(color);
                color = argv[7];
                analyzer.G = stoi(color);
                color = argv[8];
                analyzer.B = stoi(color);
                }

            choice = argv[3];
            if (choice == "t") // Pro textové soubory
                {
                analyzer.Load_CSV(argv[1]); // Načtení vstupního souboru
                analyzer.Analyze_1D(IterMin, IterMax);
                }
            if (choice == "b") // Pro binární soubory
                {
                analyzer.Load_Wav(argv[1]); // Načtení vstupního souboru
                analyzer.Analyze_1D(IterMin, IterMax);
                }
            if (choice == "p") // Pro obrázky
                {
                analyzer.Load_Image(argv[1]); // Načtení vstupního souboru
                analyzer.Analyze_2D(IterMin, IterMax);
                }

            analyzer.Save_File(argv[2],  2, 0); // Uložení výstupního souboru
        }
    }
}
// Autor programu: Ing. Pavel Florián Ph.D. dovoluje použití zdrojového kódu této aplikace, pokud je splněna podmínka uvedení autorova jména.
