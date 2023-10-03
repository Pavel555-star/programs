#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <bits/stdc++.h>
#include "Huckel_calculations.h"
using namespace std; // compiler parameters: -pthread

int main(int argc, char *argv[])
{
int i, j;
static unsigned int order;
string path_file;
string line;

const double h = 6.62607015E-34;
const double c = 299796000;
const double e = 1.60217663E-19;
static double* matrix;

Huckel_calculations<double> Huckel_calculations;

if (argc == 1)
    {
    cout << "Zadejte řád Hückelovy matice: / Enter the order of Hückel matrix: / Geben Sie die Reihenfolge der Hückel Matrix ein: /"
     << " Введите порядок Хюккель матрицы:\n";
    cin >> order;
    
    matrix = new double[(order) * (order)];
    cout << "Zadejte elementy matice: / Enter the elements of matrix: / Eingeben Matrixelemente: / Ввести элементы матрицы:  \n";
    for (j = 0; j < order; j++)
        for (i = 0; i < order; i++)
            cin >> matrix[i + (order * j)];
    
    cout << "Zadali jste matici: / You have entered the matrix: / Du hast eine Matrix eingegeben: / Вы вошли в матрицу:  \n";     
        for (j = 0; j < order; j++)
            {  
            for (i = 0; i < order; i++)
                cout << matrix[i + (order * j)] <<" ";
                cout<<endl;
            }
    cout << endl;
    }
else // filename with matrix in first program command line parameter, first line - integer order of matrix, next lines - cell values
    {
    path_file = argv[1];
    ifstream myfile(path_file);
    if (myfile.is_open())
        {
        i = 0;
        getline(myfile,line);
        order = stod(line);
        matrix = new double[(order) * (order)];
        while (getline(myfile,line))
            {
            if (i < (order * order))
                matrix[i] = stod(line);
                i++;
            }
        myfile.close();
        }
    else 
        cout << "Unable to open input file";
    }
Huckel_calculations.Huckel_Determinant_check(order, matrix);
Huckel_calculations.Huckel_Determinant_solver(order, matrix);
Huckel_calculations.Huckel_Spectra_solver(); 
Huckel_calculations.Huckel_Wavefunction_Coefficient_solver(order, matrix);

cout << "Energy levels: / Energetivké hladiny: / Energieniveaus: / Уровни энергии:" << endl;
if (Huckel_calculations.determinants.size() > 0)
    for (i = Huckel_calculations.determinants.size() - 1; i >= 0; i--)
        cout << Huckel_calculations.determinants[i] << " β" <<endl;
cout << endl;
cout << "Spectral transitions: / Spektrální přechody: / Spektrale Übergänge: / Спектральные переходы:" << endl;
for (i = 0; i < Huckel_calculations.spectra.size(); i++)
    cout << Huckel_calculations.spectra[i] << " β, " << Huckel_calculations.spectra[i] * 0.52 << " eV, "
    << ((h * c * 1000000000)/(0.52 * Huckel_calculations.spectra[i] * e)) << " nm" <<endl;
cout << endl;

if (Huckel_calculations.cyclic == true)
    cout << "Molecule is cyclic. / Molekula je cyklická. / Das Molekül ist zyklisch. / Молекула циклическая" << endl;
else
    cout << "Molecule is linear. / Molekula je lineární. / Das Molekül ist linear. / Молекула линейна." << endl;
cout << endl;
cout << "Orders of atom: / Pořadí atomů: / Ordnung der Atome: / порядок атомов:" << endl;
for (i = 0; i < Huckel_calculations.sequence.size(); i++)
    cout << Huckel_calculations.sequence[i] <<endl;
cout << endl;
cout << "Values of their wavefunctions: / Hodnoty jejich vlnových funkcí: / Werte ihrer Wellenfunktionen: / значения их волновых функций:"
 << endl;
 
j = 0;
for (i = 0; i < Huckel_calculations.coefficients.size(); i++)
    {
    cout << Huckel_calculations.coefficients[i] <<endl;
    j++;
    if (j >= order)
        {
        j = 0;
        cout << endl;
        }
    }
return(0);
};
// Author of this program Ing. Pavel Florian Ph.D. licensed this source code under the the 3-Clause BSD License
