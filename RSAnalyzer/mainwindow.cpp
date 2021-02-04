#include "mainwindow.h"
#include "rsanalyzer.h"

#include <QApplication>



MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    Output_Path = "output.csv";
    IterMin = 0;
    IterMax = 0;
    Column_Number = 1;
    R = 1;
    G = 1;
    B = 1;
    AudioFormat = 2;
    RSAnalyzer analyzer;


    resize(640, 480);
    setWindowTitle(tr("Program pro fraktální analýzu"));
    setToolTip(tr("Tento program analyzuje CSV, audio a obrazové soubory."));
    setWindowIcon(QIcon(":/img/MainIcon"));
    central = new QWidget(this);

    layout = new QGridLayout();
    central->setLayout(layout);

    csvButton = new QPushButton(tr("Otevřít textový/CSV soubor"));
    csvButton->setToolTip(tr("Vybere jednosloupcový CSV soubor pro fraktální analýzu R/S metodou a mřížkovou metodou."));

    audioButton = new QPushButton(tr("Otevřít audio soubor"));
    audioButton->setToolTip(tr("Vybere Audio soubor pro fraktální analýzu  R/S metodou a mřížkovou metodou."));

    pictureButton = new QPushButton(tr("Otevřít obrázek"));
    pictureButton->setToolTip(tr("Vybere obrázkový soubor pro fraktální analýzu  R/S metodou a mřížkovou metodou."));

    settingsButton = new QPushButton(tr("Nastavení analýzy"));
    pictureButton->setToolTip(tr("Zde můžete nastavit rozsah analýzy, zahrnutí barev do analýzy a výstupní soubor."));

    helpButton = new QPushButton(tr("Nápověda"));
    pictureButton->setToolTip(tr("Zobrazí informace o programu a jeho použití."));

    TextLabel = new QLabel(this);
    ImageLabel = new QLabel(this);

    layout->addWidget(csvButton, 0, 0, 1, 1);
    layout->addWidget(audioButton, 0, 1, 1, 1);
    layout->addWidget(pictureButton, 0, 2, 1, 1);
    layout->addWidget(settingsButton, 0, 3, 1, 1);
    layout->addWidget(helpButton, 0, 4, 1, 1);
    layout->addWidget(TextLabel, 1, 0, 1, 5);
    layout->addWidget(ImageLabel, 2, 0, 3, 5);


    setCentralWidget(central);
    connect(csvButton, SIGNAL(clicked()), this, SLOT(Open_CSV()));
    connect(audioButton, SIGNAL(clicked()), this, SLOT(Open_Audio()));
    connect(pictureButton, SIGNAL(clicked()), this, SLOT(Open_Picture()));
    connect(settingsButton, SIGNAL(clicked()), this, SLOT(Settings()));
    connect(helpButton, SIGNAL(clicked()), this, SLOT(Help()));



}

MainWindow::~MainWindow()
{
    if (csvButton != NULL) {delete csvButton;}
    if (audioButton != NULL) {delete audioButton;}
    if (pictureButton != NULL) {delete pictureButton;}


    if (layout != NULL) {delete layout;}
    if (central != NULL) {delete central;}
}

void MainWindow::Open_CSV()
{
        File_Path = QFileDialog::getOpenFileName(this,
            tr("Otevřít textový soubor"), QStandardPaths::standardLocations(QStandardPaths::DocumentsLocation).first(), tr("Textové soubory (*.csv *.txt *.inf, ...)"));
        analyzer.Load_CSV(File_Path, Column_Number);
        analyzer.Analyze_1D(IterMin, IterMax);
        analyzer.Save_File(Output_Path,  2, 0);
        Result_String(analyzer.Output_x, analyzer.Output_RS, analyzer.Output_2);
        TextLabel->setText(Result);
        resize(640, 480);


}

void MainWindow::Open_Audio()
{
    File_Path = QFileDialog::getOpenFileName(this,
        tr("Otevřít nekódovaný audio soubor"), QStandardPaths::standardLocations(QStandardPaths::MusicLocation).first(), tr("Audio soubory (*.wav *.raw *.pcm, ...)"));
    analyzer.CountBytes = AudioFormat;
    analyzer.Load_Wav(File_Path);
    analyzer.Analyze_1D(IterMin, IterMax);
    analyzer.Save_File(Output_Path,  2, 0);
    Result_String(analyzer.Output_x, analyzer.Output_RS, analyzer.Output_2);
    TextLabel->setText(Result);

}

void MainWindow::Open_Picture()
{
    File_Path = QFileDialog::getOpenFileName(this,
        tr("Otevřít obrázkový soubor"), QStandardPaths::standardLocations(QStandardPaths::MusicLocation).first(), tr("Obrázkové soubory (*.bmp *.gif *.jpg *.jpeg *.png *.pbm *.pgm *.ppm *.xbm *.xpm   ...)"));

    QPixmap pixmap(File_Path);
    int n = pixmap.width();
    int m = pixmap.height();
    int o = (n/640) + 1;
    pixmap = pixmap.scaled((n/o), (m/o), Qt::KeepAspectRatio, Qt::SmoothTransformation);
    ImageLabel->setPixmap(pixmap);

    analyzer.R = R;
    analyzer.G = G;
    analyzer.B = B;
    analyzer.Load_Image(File_Path);  
    analyzer.Analyze_2D(IterMin, IterMax);
    analyzer.Save_File(Output_Path,  2, 0);
    Result_String(analyzer.Output_x, analyzer.Output_RS, analyzer.Output_2);
    TextLabel->setText(Result);
}
void MainWindow::Settings()
{
    QDialog *settings = new QDialog;
    QGridLayout *set_grid = new QGridLayout();;
    settings->setLayout(set_grid);

    QLabel *Range = new QLabel(tr("Zvolte minimum a maximum intervalů/ délek stran čtverců: "));
    QLabel *Minimum = new QLabel(tr("Minimum: "));
    QLabel *Maximum = new QLabel(tr("Maximum: "));

    QLabel *Color = new QLabel(tr("Zvolte počty zahrnutí jednotlivých barev do analýzy: "));
    QLabel *Red = new QLabel(tr("Červená: "));
    QLabel *Green = new QLabel(tr("Zelená: "));
    QLabel *Blue = new QLabel(tr("Modrá: "));
    QLabel *Size = new QLabel(tr("Zvolte počty bytů čísel binárních souborů (1-8): "));
    QLabel *Column = new QLabel(tr("Zvolte číslo sloupce v CSV souboru: "));

    QDoubleSpinBox *min = new QDoubleSpinBox();
    QDoubleSpinBox *max = new QDoubleSpinBox();
    QDoubleSpinBox *size = new QDoubleSpinBox();
    QDoubleSpinBox *column = new QDoubleSpinBox();
    QComboBox *red = new QComboBox();
    QComboBox *green = new QComboBox();
    QComboBox *blue = new QComboBox();

    min->setRange(0, 1000000000);
    max->setRange(0, 1000000000);
    size->setRange(1, 8);
    column->setRange(1, 8);
    red->addItem("0");
    red->addItem("1");
    green->addItem("0");
    green->addItem("1");
    blue->addItem("0");
    blue->addItem("1");

    QPushButton *saveButton = new QPushButton(tr("Soubor s výsledky"));
    saveButton->setToolTip(tr("Zde můžete zvolit soubor pro výsupní data."));

    set_grid->addWidget(Range, 0, 0, 1, 1);
    set_grid->addWidget(Minimum, 1, 0, 1, 1);
    set_grid->addWidget(Maximum, 2, 0, 1, 1);
    set_grid->addWidget(Color, 3, 0, 1, 1);
    set_grid->addWidget(Red, 4, 0, 1, 1);
    set_grid->addWidget(Green, 5, 0, 1, 1);
    set_grid->addWidget(Blue, 6, 0, 1, 1);
    set_grid->addWidget(Size, 7, 0, 1, 1);
    set_grid->addWidget(Column, 8, 0, 1, 1);

    set_grid->addWidget(min, 1, 1, 1, 1);
    set_grid->addWidget(max, 2, 1, 1, 1);
    set_grid->addWidget(red, 4, 1, 1, 1);
    set_grid->addWidget(green, 5, 1, 1, 1);
    set_grid->addWidget(blue, 6, 1, 1, 1);
    set_grid->addWidget(size, 7, 1, 1, 1);
    set_grid->addWidget(column, 8, 1, 1, 1);

    set_grid->addWidget(saveButton, 9, 0, 1, 1);
    connect(saveButton, SIGNAL(clicked()), this, SLOT(Save_Dialog()));
    settings->exec();
    IterMin = min->value();
    IterMax = max->value();
    AudioFormat = size->value();
    Column_Number = column->value();
    if (red->currentText() == "0")
            R =0;
    if (red->currentText() == "1")
            R =1;
    if (green->currentText() == "0")
            G =0;
    if (green->currentText() == "1")
            G =1;
    if (blue->currentText() == "0")
            G =0;
    if (blue->currentText() == "1")
            G =1;
    if (Range != nullptr) delete Range;
    if (Minimum != nullptr) delete Minimum;
    if (Maximum != nullptr) delete Maximum;
    if (Color != nullptr) delete Color;
    if (Red != nullptr) delete Red;
    if (Green != nullptr) delete Green;
    if (Blue != nullptr) delete Blue;
    if (Size != nullptr) delete Size;
    if (Column != nullptr) delete Column;

    if (min != nullptr) delete min;
    if (max != nullptr) delete max;
    if (red != nullptr) delete red;
    if (green != nullptr) delete green;
    if (blue != nullptr) delete blue;
    if (size != nullptr) delete size;
    if (column != nullptr) delete column;

    if (saveButton != nullptr) delete saveButton;
}
void MainWindow::Help()
{
    QMessageBox message;
    message.about(this, tr ("Nápověda"), tr("Toto je program pro fraktální analýzu textových, binárních a obrázkových souborů.\n"
                                      "R/S analýza probíhá vydělením hodnoty rozsahu hodnot směrodatnou odchylkou (R/S) hodnot v daném intervalu. \n"
                                      "Analýza mřížkovou metodou probíhá výpočtem počtu čtverců (box count) na křivce hodnot v zadaném intervalu \n"
                                      "nebo počtu čtverců dané velikosti dosahujících minimálně průměrného jasu.\n"
                                      "Následně jsou vytvořeny logaritmické závislosti se základem 2 pro R/S a box count na\n"
                                      "celkovém počtu intervalů nebo čtverců (N) a uloženy ve formátu CSV do zvoleného souboru, nebo\n"
                                      "do souboru output.csv ve složce s programem. V nastavení je možné zvolit výstupní soubor,\n"
                                      " minimum a maximum počtu analyzovaných intervalů délek a barvy zahrnuté do analýzy obrázků.\n"
                                            "\n"
                                            "Program naprogramoval Ing. Pavel Florián Ph.D. v roce 2021."));
}
void MainWindow::Save_Dialog()
{
    Output_Path = QFileDialog::getSaveFileName(this,
        tr("Zvolit textový soubor"), QStandardPaths::standardLocations(QStandardPaths::MusicLocation).first(), tr("Textové soubory (*.csv *.txt *.inf, ...)"));
}

void MainWindow::Result_String(QVector<double> x, QVector<double> y, QVector<double> y_2)
{
    int i;
    int j = 0;
    int k;
    QString buffer1;
    QString Array = "";

    if (x.size() > 0)
        Result = Output_Path +":  [log2(N); log2(R/S); log2(box count)] \n";
    for (i = 0; i < x.size(); i++)
        {
            Array.append("[");
            buffer1.setNum(x[i]);
            Array.append(buffer1);
            Array.append("; ");
            buffer1.setNum(y[i]);
            Array.append(buffer1);
            Array.append("; ");
            buffer1.setNum(y_2[i]);
            Array.append(buffer1);
            Array.append("]");
            j = 24 -  Array.size();
            for (k = 0; k < j; k++)
                Array.append("  ");

            if ((i % 5) == 4)
            {
                Array.append("\n");
            }
            Result.append(Array);
            Array = "";
        }
}
