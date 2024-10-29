#include <iostream>
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include "basis_set_calculations_DFT.h"
#include "Huckel_calculations.h"
/*
#include "Visualization_3D.h"
*/
using namespace std;
// compiler parameters: -pthread -ffast-math -fno-finite-math-only

int main(int argc, char *argv[])
{
unsigned int i, j;
unsigned int count_electrons;
unsigned int count_pi_electrons;
unsigned int lenght_order;
unsigned int max_iterations;
unsigned int parameter_size;

double success;
double fidelity;
double Hamiltonian;

string input ="";
vector<double> values;
vector<double> spin_values;
vector<double> spin_density_matrix;
double* Huckel_matrix;
unsigned int* Huckel_matrix_order;
vector<unsigned int> Huckel_electrons_to_atom_numbers;

vector<unsigned int> number_electrons;
vector<unsigned int> n;
vector<unsigned int> l;
vector<int> m;
vector<double> spins;


basis_set_calculations_DFT<double> basis_set_calculations_DFT;
Huckel_calculations<double> Huckel_calculations;

cout<<"     ...                    ..     (@@@@@@@@@@@@@@(.                                                " <<endl;
cout<<"                                @@@@@@@@@@@@@@@@@@@@@@@@                                            " <<endl;
cout<<"                             @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#                                      " <<endl;
cout<<"                           &@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@               ..               " <<endl;
cout<<"                          @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@( @@@@       .. #@           " <<endl;
cout<<"                         @@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@ @@@@@@@@@@@@@       " <<endl;
cout<<"                        @@@@@@@@@@/@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@&@@@@@@@@@   " <<endl;
cout<<"                       @@@@@@@,@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@(.(@ " <<endl;
cout<<"                      @@@@@@@@@&@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@  @@@//@@@@@@" <<endl;
cout<<"  (@@@@              @@@@@@@@@ @@@@@@@@@@ @@@@@@@@@@@@@ @@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@ @// @@ @@@  " <<endl;
cout<<" @@@@@@@@           @@@@@@@@,@@@@@@@@@   %@@@@@@@@@ @@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@           @@    " <<endl;
cout<<"&@@@@@@@@        .@@@@@@@@,@@@&        @@@@@@@@@     @@@@@@@@@@#@@@@&@@@@@@@@@@@@                   " <<endl;
cout<<" @@@@@@@@@@@@@@@@@@@@@@@,@@@@@           @@@@@@@@@         @.@@@@@@@@@@ @@@@@@ @@ @@@/&             " <<endl;
cout<<"  (@@@@@@@@@@@@@@@@@@@@%@@@@@@@@@%          @@@@@@@@@@@ @                  @@@%@@@@@                " <<endl;
cout<<"     #@@@@@@@@@@@@@@&                                                                               " <<endl;
cout<<"                                                                                                    " <<endl;
cout<< endl;

cout<<"Program for basis set of wavefunctions calculation of atoms, ions and molecules." <<endl;
cout<<"Program pro výpočty základní sady vlnových funkcí atomů, iontů a molekul." <<endl;
cout<<"Programm zur Grundsatz von Wellenfunktionen Berechnung von Atomen, Ionen und Molekülen." <<endl;
cout<<"Программа для расчеты базовый набор волновых функций уровней атомов, ионов и молекул." <<endl;
cout<< endl;

if (argc > 9) // load list of excitations from parameters
    {
    for (i = 5; i < argc - 4; i+=5)
        {
        number_electrons.push_back(stoi(argv[i]));
        n.push_back(stoi(argv[i + 1]));
        l.push_back(stoi(argv[i + 2]));
        m.push_back(stoi(argv[i + 3]));
        spins.push_back(stod(argv[i + 4]));
        }
    }

if (argc < 4)
    {
    cout<<"Enter the input string in format: / Zadejte vstupní řetězec ve formátu: /" << endl;
    cout<<"Geben Sie eine Eingabezeichenfolge im Format ein: / Введите входную строку в формате:" <<endl;
    cout<<"atom1_atom2_[x2_y2_z2]_atom3_[x3_y3_z3]_{bondcount1}_{bondcount2}_(V1 eV)_(V2 eV)_(V3 eV)" <<endl;
    cout<<"<ligand1_configuration-number_highspin>" <<endl;
    cout<<"e.g. H_Li_[2.00_0.00_0.00]_{1_2_1}_(0.00)_(0.00) bondcount1 = atom1_atom2_bondcount  x2 = 2.00 Bohr_radius" << endl;
    cout<<"e.g. Fe_<1_1_1> configuration-number = 1-7 (1: orthombic, 2: pentagonal bipyramidal, 3: square antiprismatic, "
    << "4: square planar, 5: square pyramidal, 6: tetrahedral, 7: trigonal bipyramidal, highspin: 0 - low spin complexes" << endl;
    cout<<"atom1_atom2_[x2_y2_z2]_atom3_[x3_y3_z3]_{početvazeb1}_{početvazeb2}" << endl;
    cout<<"atom1_atom2_[x2_y2_z2]_atom3_[x3_y3_z3]_{zahlbindungen1}_{zahlbindungen2}" << endl;
    cout<<"атом1_атом2_[x2_y2_z2]_атом3_[x3_y3_z3]_{номерсвязь1}_{номерсвязь12}" << endl;
    cin>> input;
    
    cout<<"Enter the maximum count of iteration: / Zadejte maximální počet opakování: /" << endl;
    cout<<"Geben Sie die maximale Anzahl der Wiederholungen ein: / Введите максимальное количество повторений: /" << endl;
    cin>>max_iterations;
    
    cout<<"Enter desired precision (0-1): / Zadejte požadovanou přesnost (0-1): /" << endl;
    cout<<"Geben Sie die gewünschte Genauigkeit ein (0-1): / Введите желаемую точность (0-1): /" << endl;
    cin>>fidelity;
    
    cout<<"Enter the size factor of the wave functions (e.g. 50): /" <<endl;
    cout<<"Zadejte faktor velikosti vlnových funkcí (například 50): /" << endl;
    cout<<"Geben Sie den Betragsfaktor der Wellenfunktion ein (z. B. 50): /" <<endl;
    cout<<"Введите коэффициент величины волновой функции (например, 50): /" <<endl;
    cin>>lenght_order;
    }
else
    {
    input = argv[1];
    max_iterations = stoi(argv[2]);
    fidelity = stod(argv[3]);
    lenght_order = stoi(argv[4]);
    }
    
success = basis_set_calculations_DFT.String_to_advanced_parameters(input, lenght_order, false,
nullptr, nullptr, nullptr);
if (success == 0)
    {
    for (i = 0; i < number_electrons.size(); i++) // create excitations from command line parameters
        {
        if (n[i] > 0 and n[i] < 8)
            basis_set_calculations_DFT.results.n[number_electrons[i]] = n[i];
        if (l[i] < n[i])
            basis_set_calculations_DFT.results.l[number_electrons[i]] = l[i];
        if (abs(m[i]) <= l[i])
            basis_set_calculations_DFT.results.m[number_electrons[i]] = m[i];
        if (spins[i] == 0.5 or spins[i] == -0.5)
            basis_set_calculations_DFT.results.spins[number_electrons[i]] = spins[i];
        
        basis_set_calculations_DFT.results.spin_paired[number_electrons[i]] = spins[i];
        }
    Hamiltonian = basis_set_calculations_DFT.Execute_PBE_VQE(max_iterations, fidelity, lenght_order, false,
    &values, &spin_density_matrix, &spin_values);
    Huckel_matrix_order = new unsigned int[1];
    count_electrons = basis_set_calculations_DFT.results.n.size();
    count_pi_electrons = 0;
    for (i = 0; i < count_electrons; i++)
        if (basis_set_calculations_DFT.results.pi_bonding[i] >= 0)
            count_pi_electrons = count_pi_electrons + 1;
        
    Huckel_matrix = new double[count_pi_electrons * count_pi_electrons]; // Initializing Huckel matrix
    basis_set_calculations_DFT.Calculate_Huckel_Matrix(Huckel_matrix, Huckel_matrix_order,
    Huckel_electrons_to_atom_numbers);
    basis_set_calculations_DFT.Detect_symetry_information(&basis_set_calculations_DFT.symetry_axes_parameters,
    &basis_set_calculations_DFT.symetry_planes_parameters);
    }
else
    count_electrons = 0;

cout<<"Input string: / Vstupní řetězec: / eine Eingabezeichenfolge: / Входную строку: " << input <<endl;

cout<< "The number of iterations: / Počet pakování: / die Anzahl der Wiederholungen: / количество повторений: "
<< basis_set_calculations_DFT.iterations <<endl;

cout<< "Total Hamiltonian: / Celkový Hamiltonián: / Gesamt Hamiltonian: / общее гамильтониан: "
<< (Hamiltonian/basis_set_calculations_DFT.e) << " eV" <<endl;

cout<< "Electronic Hamiltonian: / Elektronový Hamiltonián: / Elektronisch Hamiltonian: / электронный гамильтониан: "
<< ((Hamiltonian - basis_set_calculations_DFT.nucleus_repulsive_energy)/basis_set_calculations_DFT.e) << " eV" <<endl;
cout<<endl;

cout<< "Center of mass: / Těžiště: / Schwerpunkt: / центр гравитации:"<< " x = "<< basis_set_calculations_DFT.x_center_of_mass
<< " y = "<< basis_set_calculations_DFT.y_center_of_mass << " z = " << basis_set_calculations_DFT.z_center_of_mass << " Bohr radius" <<endl;
cout<< "Molecule dipole moment: / Molekulový dipólový moment: / Molekulares Dipolmoment: / Молекулярный дипольный момент: "
<< basis_set_calculations_DFT.dipoles.sum_dipole_moment_Debye << " Debye "<<endl;

cout<< "Energy levels: / Energetické hladiny: / Energieniveaus: / уровни энергии:"<<endl;
for  (i = 0; i < count_electrons; i++)
    {
    cout<<"n = "<< basis_set_calculations_DFT.results.n[i] <<", l = "<< basis_set_calculations_DFT.results.l[i] <<
    ", m = " << basis_set_calculations_DFT.results.m[i] << ", s = " << basis_set_calculations_DFT.results.spins[i] <<
    ", Z = " << basis_set_calculations_DFT.results.Z[i]; 
    if (basis_set_calculations_DFT.results.bonding[i] >= 0)
        cout  << " J = " << basis_set_calculations_DFT.resonance_integral_matrix[(i * count_electrons) +
        basis_set_calculations_DFT.results.bonding[i]]/basis_set_calculations_DFT.e << " eV " << ", S = "
        << abs(basis_set_calculations_DFT.overlap_integral_matrix[(i * count_electrons) + basis_set_calculations_DFT.results.bonding[i]])
        << " -> " << basis_set_calculations_DFT.results.bonding[i];
    cout << endl;
    cout << " E = " << values[i]/basis_set_calculations_DFT.e <<" eV" << endl;
    cout << " wavefunction coefficient: " << basis_set_calculations_DFT.results.wavefunction_coefficients[i]
    << ", wavefunction lenght: " << (basis_set_calculations_DFT.results.efective_radius_base[i]/
    (basis_set_calculations_DFT.results.wavefunction_lenght_multipliers[i] * basis_set_calculations_DFT.Hartree_lenght)) << " Bohr radius"
    <<endl;
    cout << "shielding coefficient: / koeficient stínění: / Abschirmkoeffizient: / коэффициент экранирования: " <<
    basis_set_calculations_DFT.results.wavefunction_lenght_multipliers[i] <<endl;
    }
cout <<endl;

cout << "Z, " << "N, " << "C, " << "n, " << "l, " << "m, " << "s, " << "Z_eff/Z" <<endl;
for  (i = 0; i < count_electrons; i++)
    {
    cout<< basis_set_calculations_DFT.results.Z[i] << ", " << basis_set_calculations_DFT.results.count_electrons[i] << ", "
    << basis_set_calculations_DFT.results.wavefunction_coefficients[i] << ", " << basis_set_calculations_DFT.results.n[i] << ", "
    << basis_set_calculations_DFT.results.l[i] << ", " << basis_set_calculations_DFT.results.m[i] << ", "
    << basis_set_calculations_DFT.results.spins[i] << ", " << basis_set_calculations_DFT.results.wavefunction_lenght_multipliers[i] <<endl;
    }
/*
cout << "Spectra: / Spektra: / Spektren: / спектры:" <<endl;
for (i = 0; i < basis_set_calculations_DFT.electron_spectra.size(); i++)
    {
    cout << " E = " << (basis_set_calculations_DFT.electron_spectra[i]/basis_set_calculations_DFT.e) << " eV; " 
    << ", λ = " << (basis_set_calculations_DFT.h * basis_set_calculations_DFT.c / (basis_set_calculations_DFT.electron_spectra[i])
    * 1000000000)
    << " nm" <<endl;
    }
cout <<endl;
*/
cout << "Spin energy levels: / Spinové energetické hladiny: / Spin-Energieniveaus: / уровни спиновой энергии:" <<endl;
if (spin_values.size() > 0)
    for (i = spin_values.size() -1; i >= 0; i--)
        {
        cout << spin_values[i] <<endl;
        if (i == 0)
            break;
        }
cout <<endl;

if  (success == 0 and Huckel_matrix_order[0] >= 2
    and Huckel_calculations.Huckel_Determinant_check(Huckel_matrix_order[0], Huckel_matrix) == 0)
    {
    Huckel_calculations.Huckel_Determinant_solver(Huckel_matrix_order[0], Huckel_matrix);
    Huckel_calculations.Huckel_Spectra_solver();
    
    cout << "Energy levels: / Energetivké hladiny: / Energieniveaus: / Уровни энергии:" << endl;
    if (Huckel_calculations.determinants.size() > 0)
        for (i = 0; i < Huckel_calculations.determinants.size(); i++)
            {
            cout << Huckel_calculations.determinants[i] << " β" <<endl;
            if (i == 0)
            break;
            }
    cout << endl;

    cout << "Spectral transitions: / Spektrální přechody: / Spektrale Übergänge: / Спектральные переходы:" << endl;
    for (i = 0; i < Huckel_calculations.spectra.size(); i++)
        {
        cout << Huckel_calculations.spectra[i] << " β, " << Huckel_calculations.spectra[i] * 0.52 << " eV, "
        << ((basis_set_calculations_DFT.h * basis_set_calculations_DFT.c * 1000000000)/
        (0.52 * Huckel_calculations.spectra[i] * basis_set_calculations_DFT.e)) << " nm" <<endl;
        }
    cout << endl;
    
    if (Huckel_calculations.cyclic == true)
        cout << "Pi-stucture is cyclic. / Pi struktura je cyklická. / Die Pi-Struktur ist zyklisch. / Структура Пи циклическая." << endl;
    else
        cout << "Pi-stucture is non-cyclic. / Pi struktura je necyklická. / Die Pi-Struktur ist nichtzyklisch." 
        << "/ Структура Пи нециклический." << endl;
    cout << endl;
    }
if (success == 0)
    {
    if (basis_set_calculations_DFT.symetry_axes_parameters.symetry_center == true)
        cout << "Center of mass is the center of symmetry. / Těžiště je centrem symetrie." << endl
        << " /Der Schwerpunkt ist das Symmetriezentrum / Центр масс это центр симметрии." << endl;
    else
        cout << "Center of mass is not the center of symmetry. / Těžiště není centrem symetrie." << endl
        << " /Der Schwerpunkt ist nich das Symmetriezentrum / Центр масс не является центром симметрии." << endl;
    
    cout << "Number of axes of symetry: / Počet os symetrie: / Anzahl der Symmetrieachsen: / Количество осей симметрии: "
    << basis_set_calculations_DFT.symetry_axes_parameters.u_x.size() << endl;
    cout << "Number of planes of symetry: / Počet rovin symetrie: / Anzahl der Symmetrieebenen: / Количество плоскостей симметрии: "
    << basis_set_calculations_DFT.symetry_planes_parameters.a.size() << endl;
    
    if (basis_set_calculations_DFT.symetry_axes_parameters.linear == true)
        cout << "Stucture is linear. / Struktura je lineární. / Die Struktur ist linear. / Структура линейна." << endl;
    else
        cout << "Stucture is non-linear. / Struktura je nelineární. / Die Struktur ist nichtlinear. / Структура нелинейная." << endl;
        
    if (basis_set_calculations_DFT.symetry_planes_parameters.planar == true)
        cout << "Stucture is planar. / Struktura je planární. / Die Struktur ist planar. / Структура плоская." << endl;
    else
        cout << "Stucture is non-planar. / Struktura není planární. / Die Struktur ist nicht planar. / структура не плоская." << endl;
    
    if (Huckel_matrix_order[0] >= 2)
        delete Huckel_matrix;
        
    delete Huckel_matrix_order;
    }
cout << "Korelační energie elektronů: / Correlation energies of electrons: / Elektronenkorrelationsenergie: / электронная корреляционная энергия:" << endl;
for  (i = 0; i < basis_set_calculations_DFT.correlation_energies.size(); i++)
    {
    cout << basis_set_calculations_DFT.correlation_energies[i]/basis_set_calculations_DFT.e << " eV" << endl;
    }
cout << "Výměnné energie elektronů: / Exchange energies of electrons: / обменная: / электронная обменная энергия:" << endl;
for  (i = 0; i < basis_set_calculations_DFT.exchange_energies.size(); i++)
    {
    cout << basis_set_calculations_DFT.exchange_energies[i]/basis_set_calculations_DFT.e << " eV" << endl;
    }
/*
Visualization_3D.Compute_densities(&basis_set_calculations_DFT.results.probabilities, &basis_set_calculations_DFT.index_atoms, lenght_order,
&basis_set_calculations_DFT.results.spins, &basis_set_calculations_DFT.results.spin_paired, false);
Visualization_3D.Generate_coordinates(basis_set_calculations_DFT.vector_lenght, &basis_set_calculations_DFT.results.x,
&basis_set_calculations_DFT.results.y, &basis_set_calculations_DFT.results.z);
Visualization_3D.Compute_densities_weights_95(&Visualization_3D.atoms_densities_list, &Visualization_3D.weights_95, lenght_order, true);
Visualization_3D.Generate_2D_cross_section(0 , 0);
average = Visualization_3D.averages[0];
x = Visualization_3D.x_pixels_list[0];
y = Visualization_3D.y_pixels_list[0];

overall_density = Visualization_3D.Cross_sections_2D[0];

for (i = 0; i < y; i++)
    {
    for (j = 0; j < x; j++)
        {
        if (overall_density[j + (i * x)] > 0)
            cout << "#";
        else
            cout << " ";
        }
    cout << endl;
    }
*/
return(0);
};
/*
Author of this source code Ing. Pavel Florian Ph.D. licensed this source code under the the Apache License:

Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/
*/
