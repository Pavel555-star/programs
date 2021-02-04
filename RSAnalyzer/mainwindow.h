#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QPushButton>
#include <QGridLayout>
#include <QLabel>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QStandardPaths>

#include <QString>
#include <QVector>
#include <QPixmap>

#include <QMessageBox>
#include <QWidget>
 #include <QDialog>
#include <QMainWindow>

#include "rsanalyzer.h"


class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    QWidget * central;
    QGridLayout *layout;
    QPushButton *csvButton;
    QPushButton *audioButton;
    QPushButton *pictureButton;
    QPushButton *settingsButton;
    QPushButton *helpButton;

    QLabel *TextLabel;
    QLabel *ImageLabel;

    int IterMin;
    int IterMax;
    int Column_Number;
    int R;
    int G;
    int B;
    int AudioFormat;
    QString File_Path;
    QString Output_Path;
    QString Result = "";
    RSAnalyzer analyzer;

private slots:
    void Open_CSV();
    void Open_Audio();
    void Open_Picture();
    void Settings();
    void Help();
    void Save_Dialog();
    void Result_String(QVector<double> x, QVector<double> y,  QVector<double> y_2);
};
#endif //MAINWINDOW_H
