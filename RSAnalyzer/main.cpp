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
            string min = "";
            string max = "";
            if (argc > 5)
                {
                min = argv[4];
                IterMin = stoi(min);
                max = argv[5];
                IterMax = stoi(max);
                }
            if (argv[3] = "text")
                {
                analyzer.Load_CSV(argv[1]);
                analyzer.Analyze_1D(IterMin, IterMax);
                }
            if (argv[3] = "binary")
                {
                analyzer.Load_Wav(argv[1]);
                analyzer.Analyze_1D(IterMin, IterMax);
                }
            if (argv[3] = "picture")
                {
                analyzer.Load_Image(argv[1]);
                analyzer.Analyze_2D(IterMin, IterMax);
                }
            analyzer.Save_File(argv[2],  2, 0);
        }
    }
}
// Autor programu: Ing. Pavel Florián Ph.D. dovoluje použití zdrojového kódu této aplikace, pokud je splněna podmínka uvedení autorova jména.
