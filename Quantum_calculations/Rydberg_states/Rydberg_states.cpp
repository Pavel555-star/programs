#include <iostream>
#include <string>
#include <thread> // compiler options: -pthread  -lpthread
using namespace std;

template <class T>
class Quantum_calculations
{
public:
    Quantum_calculations(T Relative_permitivity, T Kronecker_delta)
    {
    T Er = Relative_permitivity;
    T dp = Kronecker_delta; // quantum defect of p-excitons
    if (dp < 0)
        dp =0;
    if (dp >= 1)
        dp = 0.999;
    }
    const double h =  6.62607015E-34; // m2·kg·s-1; J·s
    const double me = 9.10938356E-31; // kg
    const double e = -1.602176634E-19; // C
    const double E0 = 8.8541878128E-12; // m-3·kg-1·s4·A2
    const double c =  299792458; // m·s-1
    const double pi = 3.14460551103;
    T Rydberg_energies(T energy_band_gap, T mass_hole, int count_n, T energies[])
    // hodnota energie v objemovém polovodiči, relativní permitivita, efektivní hmotnost děr, maximum hlavního kvantového čísla, výsledky
    {
    
    unsigned int i, j;

    if (reduced_mass == 1)
        reduced_mass = (me *  mass_hole)/(1 + mass_hole);
    T Hartree_energy = (reduced_mass * e * e * e * e)/(4 * E0 * E0 * Er * Er * h * h); //H = mce^4/(ε0·εr·h)^2; H - 27.21161 eV
    for (i = 0; i < (count_n - 7); i += 8)
    {   // E = Eg + m·e^4/(ε0·εr·h)^2 · 1/(n-dp)
        energies[i] = energy_band_gap - (Hartree_energy/((i + 2 - dp)*( i + 2 - dp)));
        energies[i + 1] = energy_band_gap - (Hartree_energy/((i + 3 - dp)*(i + 3 - dp)));
        energies[i + 2] = energy_band_gap - (Hartree_energy/((i + 4 - dp)*(i + 4 - dp)));
        energies[i + 3] = energy_band_gap - (Hartree_energy/((i + 5 - dp)*(i + 5 - dp)));
        energies[i + 4] = energy_band_gap - (Hartree_energy/((i + 6 - dp)*(i + 6 - dp)));
        energies[i + 5] = energy_band_gap - (Hartree_energy/((i + 7 - dp)*(i + 7 - dp)));
        energies[i + 6] = energy_band_gap - (Hartree_energy/((i + 8 - dp)*(i + 8 - dp)));
        energies[i + 7] = energy_band_gap - (Hartree_energy/((i + 9 - dp)*(i + 9 - dp)));
    }
    j = i;
    for (i = j; i < count_n; i++)
        energies[i] = energy_band_gap - (Hartree_energy/((i + 2 - dp)*( i + 2 - dp)));
    return 0;
    }
    T Rydberg_shift_energies(T lenght, T dimension, T mass_hole,  int j, int max_n, T energy_shifts[])
    // rozměr polovodičové struktury, euklidovská dimenze objektu, počet kroků hlavního kvantového čísla, číslo momentu hybnosti, výsledky
    { // |l+s| ≥ j ≥ l+s
    unsigned int i, k;
    T a_material;
    T const_part;
    T exciton_radius[8];
    T reduced_radius[8];

    if (reduced_mass == 1.00)
        reduced_mass = (me *  mass_hole)/(1 + mass_hole); // mr
    a_material = a_Bohr * 1.054 * Er * me/reduced_mass; // a_m = a_Bohr·1.054·Er·me/reduced_mass
    const_part = (h * h * pi * pi)/(2 * reduced_mass); // (h^2·pi^2)/(2·mr)

    for (i = 0; i < (max_n - 7); i += 8)
    {
        exciton_radius[0] = a_material * 0.5 * (3 * (i + 2) * (i + 2) - 2); // r = a_Bohr·1.054·Er·me/reduced_mass·0.5·(3n^2 -l(l+1))
        exciton_radius[1] = a_material * 0.5 * (3 * (i + 3) * (i + 3) - 2); // p: l = 1
        exciton_radius[2] = a_material * 0.5 * (3 * (i + 4) * (i + 4) - 2);
        exciton_radius[3] = a_material * 0.5 * (3 * (i + 5) * (i + 5) - 2);
        exciton_radius[4] = a_material * 0.5 * (3 * (i + 6) * (i + 6) - 2);
        exciton_radius[5] = a_material * 0.5 * (3 * (i + 7) * (i + 7) - 2);
        exciton_radius[6] = a_material * 0.5 * (3 * (i + 8) * (i + 8) - 2);
        exciton_radius[7] = a_material * 0.5 * (3 * (i + 9) * (i + 9) - 2);
        reduced_radius[0] = lenght - (2 * exciton_radius[0]); // dL = L - 2r
        reduced_radius[1] = lenght - (2 * exciton_radius[1]);
        reduced_radius[2] = lenght - (2 * exciton_radius[2]);
        reduced_radius[3] = lenght - (2 * exciton_radius[3]);
        reduced_radius[4] = lenght - (2 * exciton_radius[4]);
        reduced_radius[5] = lenght - (2 * exciton_radius[5]);
        reduced_radius[6] = lenght - (2 * exciton_radius[6]);
        reduced_radius[7] = lenght - (2 * exciton_radius[7]);
        energy_shifts[i] = const_part * (double(j)/reduced_radius[0]) * (double(j)/reduced_radius[0]) * (3 - dimension); // (h^2·pi^2)/(2·mr)·(j/dL)^2
        energy_shifts[i + 1] = const_part * (double(j)/reduced_radius[1]) * (double(j)/reduced_radius[1]) * (3 - dimension);
        energy_shifts[i + 2] = const_part * (double(j)/reduced_radius[2]) * (double(j)/reduced_radius[2]) * (3 - dimension);
        energy_shifts[i + 3] = const_part * (double(j)/reduced_radius[3]) * (double(j)/reduced_radius[3]) * (3 - dimension);
        energy_shifts[i + 4] = const_part * (double(j)/reduced_radius[4]) * (double(j)/reduced_radius[4]) * (3 - dimension);
        energy_shifts[i + 5] = const_part * (double(j)/reduced_radius[5]) * (double(j)/reduced_radius[5]) * (3 - dimension);
        energy_shifts[i + 6] = const_part * (double(j)/reduced_radius[6]) * (double(j)/reduced_radius[6]) * (3 - dimension);
        energy_shifts[i + 7] = const_part * (double(j)/reduced_radius[7]) * (double(j)/reduced_radius[7]) * (3 - dimension);
    }
    k = i;
    for (i = k; i < max_n; i++)
    {// Cm = 3 R (hv/kt)^2 exp(hv/kt)/((exp(hv/kt)-1)^2)
        exciton_radius[0] = a_material * 0.5 * (3 * (i + 2) * (i + 2) - 2);
        reduced_radius[0] = lenght - 2 * exciton_radius[0];
        energy_shifts[i] = const_part * (j/reduced_radius[0]) * (j/reduced_radius[0]);
    }
    return 0;
    }
private:
    T a_Bohr = (E0 * h * h)/(pi * me * e * e); // Bohr radius
    T Er = 1.00;
    T reduced_mass = 1.00;
    T dp;
};
int main(int argc, char *argv[]) {
double  band_gap;
double  relative_permitivity; // Cu2O - 9.8
double  mass_hole; // Cu2O - 0.68
double lenght;
double dimension;
double Kronecker_delta; // Cu2O - 0.23, 2D quantum wells - 0.5
double Hartree_energies[24];
double differences_1[24];
double differences_3[24];
double energies[24];
const double h =  6.62607015E-34; // m2·kg·s-1; J·s
const double c =  299792458; // m·s-1
const double e = -1.602176634E-19; // J
int i = 24;
int k;
if (argc > 6)
{
    band_gap = stod(argv[1]);
    band_gap *= -e;
    relative_permitivity = stod(argv[2]);    
    mass_hole = stod(argv[3]);
    lenght = stod(argv[4]);
    lenght /= 100000000.00;
    dimension = stod(argv[5]);
    Kronecker_delta = stod(argv[6]);
}
else
{
cout << "Zadejte energii zakázaného pásu (eV): / Enter the energy of the band gap (eV): / Geben Sie die Energie der Bandlücke ein (eV):  / Введите энергию запрещенной зоны (эв): \n";
cin >> band_gap;
band_gap *= -e;
cout << "Zadejte relativní permitivitu: / Enter the relative permitivity: / Geben Sie die relative Permittivität ein:  / Введите относительную диэлектрическую проницаемость: \n";
cin >> relative_permitivity;
cout << "Zadejte efektivní hmotnost díry: / Enter the effective weight of the hole: / Geben Sie das effektive Gewicht des Lochs ein:  / Введите эффективный вес отверстия: \n";
cin >> mass_hole;
cout << "Zadejte velikost částic (nm): / Enter the particle size (nm): / Geben Sie die Partikelgröße (nm) ein: / Введите размер частиц (нм):  \n";
cin >> lenght;
lenght /= 1000000000.00;
cout << "Zadejte dimenzi částic: / Enter the particle dimension: / Geben Sie die Partikeldimension ein: / Введите размерность частицы: \n";
cin >> dimension;

cout << "Zadejte Kroneckerovu deltu (0-1): / Enter the Kronecker delta: / Geben Sie die Kronecker delta ein (0-1): / Введите дельта Кронекера (0-1): \n";
cin >> Kronecker_delta;
}
Quantum_calculations<double> * Quantum_calculationsPtr = new Quantum_calculations<double>(relative_permitivity, Kronecker_delta);
thread t1(&Quantum_calculations<double>::Rydberg_energies, Quantum_calculationsPtr, band_gap, mass_hole, i, Hartree_energies);
thread t2(&Quantum_calculations<double>::Rydberg_shift_energies, Quantum_calculationsPtr, lenght, dimension, mass_hole, 1, 24,
differences_1);
thread t3(&Quantum_calculations<double>::Rydberg_shift_energies, Quantum_calculationsPtr, lenght, dimension, mass_hole, 3, 24,
differences_3);
t1.join();
t2.join();
t3.join();
delete Quantum_calculationsPtr;

for (k = 0; k < i; k++)
{
    cout << "n = " << (k + 2) << " E = " << (Hartree_energies[k]/(-e)) << " eV " << "f= " << (Hartree_energies[k]/h) << " Hz "
    << "l= " << (h * c / Hartree_energies[k] * 1000000000.00) << " nm \n";
}
cout << "j = 1: \n";
for (k = 0; k < i; k++)
{
    cout << "n = " << (k + 2) << " dE = " << (differences_1[k]/(-e)) << " eV " << "f= " << ((Hartree_energies[k] + differences_1[k])/h) 
    << " Hz " << "l= " << (h * c / ((Hartree_energies[k] + differences_1[k])) * 1000000000.00) << " nm \n";
}
cout << "j = 3: \n";
for (k = 0; k < i; k++)
{
    cout << "n = " << (k + 2) << " dE = " << (differences_3[k]/(-e)) << " eV " << "f= " << ((Hartree_energies[k] + differences_3[k])/h) 
    << " Hz " << "l= " << (h * c / ((Hartree_energies[k] + differences_3[k])) * 1000000000.00) << " nm \n";
}
return 0;
} // Autor Ing. Pavel Florián Ph.D. dovoluje šíření tohoto programu jen s uvedením svého jména.
