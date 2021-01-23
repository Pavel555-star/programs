#include "mainwindow.h"
#include <QApplication>
#include <QTranslator>


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QTranslator translator;
    translator.load("localization_non_czech");
    a.installTranslator(&translator);

    MainWindow w;
    w.show();
    return a.exec();
}
// Autor programu: Ing. Pavel Florián Ph.D. dovoluje použití zdrojového kódu této aplikace, pokud je splněna podmínka uvedení autorova jména.
