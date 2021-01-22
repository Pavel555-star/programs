#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
// Autor programu: Ing. Pavel Florián Ph.D. dovoluje použití zdrojového kódu této aplikace, pokud je splněna podmínka uvedení autorova jména.
