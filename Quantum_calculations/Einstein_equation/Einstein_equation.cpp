#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <thread> // parameters: -pthread  -lpthread
using namespace std;
template <class T>
class Einstein_calculations
{
public:
    T Einstein_molar_capacity(T frequency, T empirical_coeficient, T T_min, T T_max, long long iter, T values[], T average[])
    {    
    long long i, j;
    T constant = h * frequency/(R/Na); // hv/k
    T step = ((T_max-T_min)/iter);
    T multistep = 8 * step;
    T sum = 0;
    T temperature[8] = {T_min, T_min + step, T_min + 2 * step, T_min + 3 * step,  T_min + 4 * step,  T_min + 5 * step,  T_min + 6 * step,  T_min + 7 * step};
    T ratio[8];
    T exponent[8];

    for (i = 0; i < (iter-7); i+=8)
    {
    ratio[0] = constant/temperature[0]; // hv/kt
    ratio[1] = constant/temperature[1];
    ratio[2] = constant/temperature[2];
    ratio[3] = constant/temperature[3];
    ratio[4] = constant/temperature[4];
    ratio[5] = constant/temperature[5];
    ratio[6] = constant/temperature[6];
    ratio[7] = constant/temperature[7];
    exponent[0] = exp(ratio[0]); // exp(hv/kt)
    exponent[1] = exp(ratio[1]);
    exponent[2] = exp(ratio[2]);
    exponent[3] = exp(ratio[3]);
    exponent[4] = exp(ratio[4]);
    exponent[5] = exp(ratio[5]);
    exponent[6] = exp(ratio[6]);
    exponent[7] = exp(ratio[7]);
    // Cm = 3 R (hv/kt)^2 exp(hv/kt)/((exp(hv/kt)-1)^2)
    values[i] = 3 * R * ratio[0] * ratio[0] * exponent[0]/((exponent[0] - 1) * (exponent[0] - 1)) * empirical_coeficient;
    values[i + 1] = 3 * R * ratio[1] * ratio[1] * exponent[1]/((exponent[1] - 1) * (exponent[1] - 1)) * empirical_coeficient;
    values[i + 2] = 3 * R * ratio[2] * ratio[2] * exponent[2]/((exponent[2] - 1) * (exponent[2] - 1)) * empirical_coeficient;
    values[i + 3] = 3 * R * ratio[3] * ratio[3] * exponent[3]/((exponent[3] - 1) * (exponent[3] - 1)) * empirical_coeficient;
    values[i + 4] = 3 * R * ratio[4] * ratio[4] * exponent[4]/((exponent[3] - 1) * (exponent[4] - 1)) * empirical_coeficient;
    values[i + 5] = 3 * R * ratio[5] * ratio[5] * exponent[5]/((exponent[5] - 1) * (exponent[5] - 1)) * empirical_coeficient;
    values[i + 6] = 3 * R * ratio[6] * ratio[6] * exponent[6]/((exponent[6] - 1) * (exponent[6] - 1)) * empirical_coeficient;
    values[i + 7] = 3 * R * ratio[7] * ratio[7] * exponent[7]/((exponent[7] - 1) * (exponent[7] - 1)) * empirical_coeficient;
    sum += values[i];
    sum += values[i + 1];
    sum += values[i + 2];
    sum += values[i + 3];
    sum += values[i + 4];
    sum += values[i + 5];
    sum += values[i + 6];
    sum += values[i + 7];
    temperature[0] += multistep;
    temperature[1] += multistep;
    temperature[2] += multistep;
    temperature[3] += multistep;
    temperature[4] += multistep;
    temperature[5] += multistep;
    temperature[6] += multistep;
    temperature[7] += multistep;
    }
    j = i;
    for (i = j; i < iter; i++)
    {
    ratio[i] = constant/temperature[i-j]; // hv/kt
    exponent[i] = exp(ratio[i]); // exp(hv/kt)
    // Cm = 3 R (hv/kt)^2 exp(hv/kt)/((exp(hv/kt)-1)^2)
    values[i] = 3 * R * ratio[i] * ratio[i] * exponent[i]/((exponent[i] - 1) * (exponent[i] - 1)) * empirical_coeficient;
    sum += values[i];
    }
    average[0] = sum/double(iter);
    return 0;
    }
    T Einstein_spectral_molar_capacity(T frequency, T T_min, T T_max, long long iter,T partial_values[], T values[], T average[])
    {

    long long i, j;
    long long k;
    T constant = h * frequency/(R/Na); // hv/k
    T step = ((T_max-T_min)/iter);
    T multistep = 8 * step;
    T sum; 
    
    T ratio[8];
    T exponent[8];
    for (k = 1; k <= spectrum_size; k++)
    {
    constant = h * frequency * (double(k)/double(spectrum_size))/(R/Na); // hv/k
    T temperature[8] = {T_min, T_min + step, T_min + 2 * step, T_min + 3 * step,  T_min + 4 * step,  T_min + 5 * step,  T_min +
     6 * step,  T_min + 7 * step};
    for (i = 0; i < (iter-7); i+=8)
    {
    ratio[0] = constant/temperature[0]; // hv/kt
    ratio[1] = constant/temperature[1];
    ratio[2] = constant/temperature[2];
    ratio[3] = constant/temperature[3];
    ratio[4] = constant/temperature[4];
    ratio[5] = constant/temperature[5];
    ratio[6] = constant/temperature[6];
    ratio[7] = constant/temperature[7];
    exponent[0] = exp(ratio[0]); // exp(hv/kt)
    exponent[1] = exp(ratio[1]);
    exponent[2] = exp(ratio[2]);
    exponent[3] = exp(ratio[3]);
    exponent[4] = exp(ratio[4]);
    exponent[5] = exp(ratio[5]);
    exponent[6] = exp(ratio[6]);
    exponent[7] = exp(ratio[7]);
    // Cm = 3 R (hv/kt)^2 exp(hv/kt)/((exp(hv/kt)-1)^2)
    partial_values[(k-1) + ((i) * spectrum_size)] = spectrum[k - 1]
     * 3 * R * ratio[0] * ratio[0] * exponent[0]/((exponent[0] - 1) * (exponent[0] - 1));
    partial_values[(k-1) + ((i + 1) * spectrum_size)] = spectrum[k - 1]
     * 3 * R * ratio[1] * ratio[1] * exponent[1]/((exponent[1] - 1) * (exponent[1] - 1));
    partial_values[(k-1) + ((i + 2) * spectrum_size)] = spectrum[k - 1]
     * 3 * R * ratio[2] * ratio[2] * exponent[2]/((exponent[2] - 1) * (exponent[2] - 1));
    partial_values[(k-1) + ((i + 3) * spectrum_size)] = spectrum[k - 1]
     * 3 * R * ratio[3] * ratio[3] * exponent[3]/((exponent[3] - 1) * (exponent[3] - 1));
    partial_values[(k-1) + ((i + 4) * spectrum_size)] = spectrum[k - 1]
     * 3 * R * ratio[4] * ratio[4] * exponent[4]/((exponent[3] - 1) * (exponent[4] - 1));
    partial_values[(k-1) + ((i + 5) * spectrum_size)] = spectrum[k - 1]
     * 3 * R * ratio[5] * ratio[5] * exponent[5]/((exponent[5] - 1) * (exponent[5] - 1));
    partial_values[(k-1) + ((i + 6) * spectrum_size)] = spectrum[k - 1]
     * 3 * R * ratio[6] * ratio[6] * exponent[6]/((exponent[6] - 1) * (exponent[6] - 1));
    partial_values[(k-1) + ((i + 7) * spectrum_size)] = spectrum[k - 1]
     * 3 * R * ratio[7] * ratio[7] * exponent[7]/((exponent[7] - 1) * (exponent[7] - 1));
    sum += partial_values[(k-1) + ((i) * spectrum_size)];
    sum += partial_values[(k-1) + ((i + 1) * spectrum_size)];
    sum += partial_values[(k-1) + ((i + 2) * spectrum_size)];
    sum += partial_values[(k-1) + ((i + 3) * spectrum_size)];
    sum += partial_values[(k-1) + ((i + 4) * spectrum_size)];
    sum += partial_values[(k-1) + ((i + 5) * spectrum_size)];
    sum += partial_values[(k-1) + ((i + 6) * spectrum_size)];
    sum += partial_values[(k-1) + ((i + 7) * spectrum_size)];
    temperature[0] += multistep;
    temperature[1] += multistep;
    temperature[2] += multistep;
    temperature[3] += multistep;
    temperature[4] += multistep;
    temperature[5] += multistep;
    temperature[6] += multistep;
    temperature[7] += multistep;
    }
    j = i;
    for (i = j; i < iter; i++)
    {
    ratio[i] = constant/temperature[i]; // hv/kt
    exponent[i] = exp(ratio[i]); // exp(hv/kt)
    // Cm = 3 R (hv/kt)^2 exp(hv/kt)/((exp(hv/kt)-1)^2)
    partial_values[(k - 1) + ((i) * spectrum_size)] = 3 * R * ratio[i] * ratio[i] * exponent[i]/((exponent[i] - 1) * (exponent[i] - 1));
    sum += partial_values[(k - 1) + ((i) * spectrum_size)];
    }
    }
    sum = sum/double(iter);
    average[0] = sum;
    {
    for (i = 0; i < (iter-7); i+=8)
    {
        values[i] = 0;
        values[i + 1] = 0;
        values[i + 2] = 0;
        values[i + 3] = 0;
        values[i + 4] = 0;
        values[i + 5] = 0;
        values[i + 6] = 0;
        values[i + 7] = 0;
        for (k = 0; k < (spectrum_size); k++)
        {
            values[i] += partial_values[k + (i * spectrum_size)];
            values[i + 1] += partial_values[k + ((i + 1) * spectrum_size)];
            values[i + 2] += partial_values[k + ((i + 2) * spectrum_size)];
            values[i + 3] += partial_values[k + ((i + 3) * spectrum_size)];
            values[i + 4] += partial_values[k + ((i + 4) * spectrum_size)];
            values[i + 5] += partial_values[k + ((i + 5) * spectrum_size)];
            values[i + 6] += partial_values[k + ((i + 6) * spectrum_size)];
            values[i + 7] += partial_values[k + ((i + 7) * spectrum_size)];
        }
    }
    j = i;
    for (i = j; i < iter; i++)
    {
        values[i] = 0;
        for (k = 0; k < (spectrum_size); k++)
        {  
        values[i] += partial_values[k + (i * spectrum_size)];
        }
    }
    }    
    return 0;
    }

    T Load_CSV()
    {
    int i = 0;
    string line;
    ifstream myfile("spectrum.csv");
    if (myfile.is_open())
    {
    while (getline(myfile,line))
    {
        if (i < spectrum_size)
            this->spectrum[i] = stod(line);
        i++;
    }
    myfile.close();
    }
    else cout << "Unable to open file";
    return 0;
    }
    const double R =  8.31446261815324; // J·K−1·mol−1
    const double Na = 6.02214076E23; // mol-1
    const double h =  6.62607015E-34; // kg·m2·s-1; J·s
    const long long spectrum_size = 100;
    double spectrum[100];   
};
int main(int argc, char *argv[]) {
double  number;
double coeficient;
double matrix_X1[200000];
double matrix_X2[200000];
double matrix_X3[200000];
double matrix_X4[200000];
double matrix_1[2000];
double matrix_2[2000];
double matrix_3[2000];
double matrix_4[2000];
double average_1[1];
double average_2[1];
double average_3[1];
double average_4[1];
int i;

if (argc > 2)
{
    number = stod(argv[1]);
    coeficient = stod(argv[2]);}
else
{
cout << "Zadejte frekvenci oscilátoru (Hz)a empirický koeficient (nebo nulu k spektrální analýze): \n";
cout << "Enter the frequency of oscilator (Hz) and empirical coeficient (or zero to spectral analysis): \n";
cout << "Geben Sie die Oszillatorfrequenz (Hz) und empirischer Koeffizient ein (oder Null für Spektralanalyse): \n";
cout << "Введите  частоту генератора (Гц) и эмпирический коэффициент (или ноль для спектрального анализа ): \n";
cin >> number;
cin >> coeficient;
}
Einstein_calculations<double> * Einstein_calculationsPtr = new Einstein_calculations<double>;
if (coeficient != 0)
{
thread t1(&Einstein_calculations<double>::Einstein_molar_capacity, Einstein_calculationsPtr, number, coeficient, 0.001, 99.9, 2000, matrix_1, average_1);
thread t2(&Einstein_calculations<double>::Einstein_molar_capacity, Einstein_calculationsPtr, number, coeficient, 100.0, 199.9, 2000, matrix_2, average_2);
thread t3(&Einstein_calculations<double>::Einstein_molar_capacity, Einstein_calculationsPtr, number, coeficient, 200.0, 299.9, 2000, matrix_3, average_3);
thread t4(&Einstein_calculations<double>::Einstein_molar_capacity, Einstein_calculationsPtr, number, coeficient, 300.0, 399.9, 2000, matrix_4, average_4);
t1.join();
t2.join();
t3.join();
t4.join();
}
else
{
Einstein_calculationsPtr->Load_CSV();
thread t1(&Einstein_calculations<double>::Einstein_spectral_molar_capacity, Einstein_calculationsPtr, number, 0.001, 99.9, 2000, matrix_X1, matrix_1, average_1);
thread t2(&Einstein_calculations<double>::Einstein_spectral_molar_capacity, Einstein_calculationsPtr, number, 100.0, 199.9, 2000, matrix_X2, matrix_2, average_2);
thread t3(&Einstein_calculations<double>::Einstein_spectral_molar_capacity, Einstein_calculationsPtr, number, 200.0, 299.9, 2000, matrix_X3, matrix_3, average_3);
thread t4(&Einstein_calculations<double>::Einstein_spectral_molar_capacity, Einstein_calculationsPtr, number, 300.0, 399.9, 2000, matrix_X4, matrix_4, average_4);
t1.join();
t2.join();
t3.join();
t4.join();
}
delete Einstein_calculationsPtr;

for (i = 0; i < 2000; i++)
    cout << "Cm = " << matrix_1[i] << " J·mol−1·K−1" << " T = " << (double(i)/20) << " K \n";
for (i = 0; i < 2000; i++)
    cout << "Cm = " << matrix_2[i] << " J·mol−1·K−1" << " T = " << (100 + double(i)/20) << " K \n";
for (i = 0; i < 2000; i++)
    cout << "Cm = " << matrix_3[i] << " J·mol−1·K−1" << " T = " << (200 + double(i)/20) << " K \n";
for (i = 0; i < 2000; i++)
    cout << "Cm = " << matrix_4[i] << " J·mol−1·K−1" << " T = " << (300 + double(i)/20) << " K \n";
cout << endl;
cout << "Cm = " << average_1[0] << " J·mol−1·K−1" << " T = " << "  0-100" << " K \n";
cout << "Cm = " << average_2[0] << " J·mol−1·K−1" << " T = " << "100-200" << " K \n";
cout << "Cm = " << average_3[0] << " J·mol−1·K−1" << " T = " << "200-300" << " K \n";
cout << "Cm = " << average_4[0] << " J·mol−1·K−1" << " T = " << "300-400" << " K \n";
return 0; }// Autor Ing. Pavel Florián Ph.D. dovoluje šíření tohoto programu jen s uvedením svého jména.
