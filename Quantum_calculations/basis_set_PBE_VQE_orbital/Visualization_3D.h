#include <thread> 
#include <vector>
#include <bits/stdc++.h>
using namespace std;

# ifndef VISUALIZATION_3D_H
# define VISUALIZATION_3D_H

template <typename T>
class Visualization_3D
{
private:
    const unsigned int max_atoms = 128;
    const vector<T> ratios = {0.95};
    unsigned int side_pixels = 0;
    T cube_side_half_lenght = 0;
    T pixels_to_unit = 0;
    T x_min = 0;
    T y_min = 0;
    T z_min = 0;
    T x_max = 0;
    T y_max = 0;
    T z_max = 0;
    
    T Compute_density_thread(const vector<T*>* densities, T* atom_density, unsigned int begin, unsigned int end, unsigned int atom_size,
    const vector<T>* spins, const vector<int>* spin_paired, bool spin_density);
    T Compute_densities_weights_thread(T* density, T ratio, T* weight_pointer, unsigned int learning_cycles, unsigned int size_order);
public:
    vector<T*> atoms_densities_list;
    vector<array<T, 3>> coordinates;
    vector<T> weights_95;
    vector<T*> Cross_sections_2D;
    vector<unsigned int> x_pixels_list;
    vector<unsigned int> y_pixels_list;
    vector<T> averages;

    T Generate_coordinates(T side_half_lenght, const vector<T>* x, const vector <T>* y, const vector<T>* z);
    T Compute_densities(const vector<T*>* densities, vector<unsigned int>* index, unsigned int size_order,
    const vector<T>* spins, const vector<int>* spin_paired, bool spin_density);
    T Compute_densities_weights_95(const vector<T*>* densities, vector<T>* weights, unsigned int size_order);
    T Generate_2D_cross_section(unsigned int direction, T coordinate);
    T Clear();
    ~Visualization_3D();
};
template <typename T>
T Visualization_3D<T>::Compute_density_thread(const vector<T*>* densities, T* atom_density, unsigned int begin, unsigned int end, unsigned int atom_size, const vector<T>* spins, const vector<int>* spin_paired, bool spin_density)
{
unsigned int i, j;
T* pointer_to_density;

// Vectorisation code
for (i = 0; i + 7 < atom_size; i++)
    {
    atom_density[i] = 0;
    atom_density[i + 1] = 0;
    atom_density[i + 2] = 0;
    atom_density[i + 3] = 0;
    atom_density[i + 4] = 0;
    atom_density[i + 5] = 0;
    atom_density[i + 6] = 0;
    atom_density[i + 7] = 0;
    }
for (i = atom_size - (atom_size % 8); i < atom_size; i++)
    {
    atom_density[i] = 0;
    }
for (i = begin; i <= end; i++)
    {
    if (spin_density == false or spin_paired->operator[](i) == -1 or spin_paired->operator[](i) < begin
    or spin_paired->operator[](i) > end)
        {
        pointer_to_density = densities->operator[](i); // Unpaired and bonding electrons are included for spin_paired[i] < begin
        if (spin_density == true and spins->operator[](i) == -0.5) // For spin density map and negative spins
            {
            for (j = 0; j + 7 < atom_size; j+=8)
                {
                atom_density[j] = atom_density[j] - pointer_to_density[j];
                atom_density[j + 1] = atom_density[j + 1] - pointer_to_density[j + 1];
                atom_density[j + 2] = atom_density[j + 2] - pointer_to_density[j + 2];
                atom_density[j + 3] = atom_density[j + 3] - pointer_to_density[j + 3];
                atom_density[j + 4] = atom_density[j + 4] - pointer_to_density[j + 4];
                atom_density[j + 5] = atom_density[j + 5] - pointer_to_density[j + 5];
                atom_density[j + 6] = atom_density[j + 6] - pointer_to_density[j + 6];
                atom_density[j + 7] = atom_density[j + 7] - pointer_to_density[j + 7];
                }
            for (j = atom_size - (atom_size % 8); j < atom_size; j++)
                {
                atom_density[j] = atom_density[j] - pointer_to_density[j];
                }
            }
        else
            {
            for (j = 0; j + 7 < atom_size; j+=8)
                {
                atom_density[j] = atom_density[j] + pointer_to_density[j];
                atom_density[j + 1] = atom_density[j + 1] + pointer_to_density[j + 1];
                atom_density[j + 2] = atom_density[j + 2] + pointer_to_density[j + 2];
                atom_density[j + 3] = atom_density[j + 3] + pointer_to_density[j + 3];
                atom_density[j + 4] = atom_density[j + 4] + pointer_to_density[j + 4];
                atom_density[j + 5] = atom_density[j + 5] + pointer_to_density[j + 5];
                atom_density[j + 6] = atom_density[j + 6] + pointer_to_density[j + 6];
                atom_density[j + 7] = atom_density[j + 7] + pointer_to_density[j + 7];
                }
            for (j = atom_size - (atom_size % 8); j < atom_size; j++)
                {
                atom_density[j] = atom_density[j] + pointer_to_density[j];
                }
            }
        }
    }// End of vectorisation code
return(0);
}
template <typename T>
T Visualization_3D<T>::Compute_densities_weights_thread(T* density, T ratio, T* weight_pointer, unsigned int learning_cycles, unsigned int size_order)
{ // perceptron function for values in 2D cross section for probability in weights 
unsigned int i, j;
unsigned int size;
T weight;
T sum;
T prediction;

size = (2 * size_order + 1) * (2 * size_order + 1) * (2 * size_order + 1);
if (ratio > 1)
    ratio = 1;
if (ratio < 0)
    ratio = 0;
    
if (ratio == 0) // weight = max value for ratio == 0
    {
    weight = 0;
    for (i = 0; i < size - (size % 8); i+=8)
        {
        if (abs(density[i]) > weight)
            weight = abs(density[i]);
        if (abs(density[i + 1]) > weight)
            weight = abs(density[i + 1]);
        if (abs(density[i + 2]) > weight)
            weight = abs(density[i + 2]);
        if (abs(density[i + 3]) > weight)
            weight = abs(density[i + 3]);
        if (abs(density[i + 4]) > weight)
            weight = abs(density[i + 4]);
        if (abs(density[i + 5]) > weight)
            weight = abs(density[i + 5]);
        if (abs(density[i + 6]) > weight)
            weight = abs(density[i + 6]);
        if (abs(density[i + 7]) > weight)
            weight = abs(density[i + 7]);
        }
    for (i = size - (size % 8); i < size; i++)
        if (abs(density[i]) > weight)
            weight = abs(density[i]);
    }
if (ratio == 1) // weight = min value for ratio == 1
    {
    weight = density[0];
    for (i = 0; i < size - (size % 8); i+=8)
        {
        if (abs(density[i]) < weight)
            weight = abs(density[i]);
        if (abs(density[i + 1]) < weight)
            weight = abs(density[i + 1]);
        if (abs(density[i + 2]) < weight)
            weight = abs(density[i + 2]);
        if (abs(density[i + 3]) < weight)
            weight = abs(density[i + 3]);
        if (abs(density[i + 4]) < weight)
            weight = abs(density[i + 4]);
        if (abs(density[i + 5]) < weight)
            weight = abs(density[i + 5]);
        if (abs(density[i + 6]) < weight)
            weight = abs(density[i + 6]);
        if (abs(density[i + 7]) < weight)
            weight = abs(density[i + 7]);
        }
    for (i = size - (size % 8); i < size; i++)
        {
        if (abs(density[i]) < weight)
            weight = abs(density[i]);
        }
    }
if (ratio > 0 and ratio < 1)
    { // compute sum initialized value of weight
    sum = 0;
    for (i = 0; i < size - (size % 8); i+=8)
        {
        sum = sum + abs(density[i]);
        sum = sum + abs(density[i + 1]);
        sum = sum + abs(density[i + 2]);
        sum = sum + abs(density[i + 3]);
        sum = sum + abs(density[i + 4]);
        sum = sum + abs(density[i + 5]);
        sum = sum + abs(density[i + 6]);
        sum = sum + abs(density[i + 7]);
        }
    for (i = size - (size % 8); i < size; i++)
        if (abs(density[i]) > weight)
            sum = sum + abs(density[i]);
            
    weight = sum/size; // setting average as initial weight
    for (i = 0; i < learning_cycles; i++) // learning loop
        {
        prediction = 0;
        for (j = 0; j < size - (size % 8); j+=8)
            {
            if (abs(density[j]) < weight)
                prediction = prediction + abs(density[j]);
            if (abs(density[j + 1]) < weight)
                prediction = prediction + abs(density[j + 1]);
            if (abs(density[j + 1]) < weight)
                prediction = prediction + abs(density[j + 2]);
            if (abs(density[j + 2]) < weight)
                prediction = prediction + abs(density[j + 1]);
            if (abs(density[j + 3]) < weight)
                prediction = prediction + abs(density[j + 3]);
            if (abs(density[j + 4]) < weight)
                prediction = prediction + abs(density[j + 4]);
            if (abs(density[j + 5]) < weight)
                prediction = prediction + abs(density[j + 5]);
            if (abs(density[j + 6]) < weight)
                prediction = prediction + abs(density[j + 6]);
            if (abs(density[j + 7]) < weight)
                prediction = prediction + abs(density[j + 7]);
            }
        for (j = size - (size % 8); j < size; j++)
            {
            if (abs(density[j]) < weight)
                prediction = prediction + abs(density[j]);
            }
        weight = weight + (sum * ratio - prediction)/size;
        weight = weight + (weight - prediction) * sum/size;
        }
    }
weight_pointer[0] = weight;
return(0);
}
template <typename T>
T Visualization_3D<T>::Generate_coordinates(T side_half_lenght, const vector<T>* x, const vector <T>* y, const vector<T>* z)
{ // Generate coordinates from x, y and z vectors
unsigned int i;
unsigned int size_x, size_y, size_z;
array<T, 3> coordinates_trinity;

size_x = x->size();
size_y = y->size();
size_z = z->size();

if (side_half_lenght <= 0) // Check the input parameters
    return(-1);
    
if (size_x != size_y or size_x != size_z)
    return(-1);
    
cube_side_half_lenght = side_half_lenght; // Set the cube_side_half_lenght

for (i = 0; i < size_x; i++) // Fill vector of atom coordinates and set min and max values
    {
    coordinates_trinity[0] = x->operator[](i);
    coordinates_trinity[1] = y->operator[](i);
    coordinates_trinity[2] = z->operator[](i);
    coordinates.push_back(coordinates_trinity);
    if (x->operator[](i) < x_min)
        x_min = x->operator[](i);
    if (y->operator[](i) < y_min)
        y_min = y->operator[](i);
    if (z->operator[](i) < z_min)
        z_min = z->operator[](i);
    }
x_min = x_min - side_half_lenght;
y_min = y_min - side_half_lenght;
z_min = z_min - side_half_lenght;
x_max = x_max + side_half_lenght;
y_max = y_max + side_half_lenght;
z_max = z_max + side_half_lenght;
return(0);
}
template <typename T>
T Visualization_3D<T>::Compute_densities(const vector<T*>* densities, vector<unsigned int>* index, unsigned int size_order,
const vector<T>* spins, const vector<int>* spin_paired, bool spin_density)
{ // Compute densities of wavefunctions/probabilities for atoms in 1 thread for each atoms.
unsigned int i, j;
unsigned int begin, end;
unsigned int count_atoms;
unsigned int size_atoms;
unsigned int begins_atoms[max_atoms];
unsigned int ends_atoms[max_atoms];
T* pointers_to_atoms[max_atoms]; 

thread t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
bool t9_flag, t10_flag, t11_flag, t12_flag, t13_flag, t14_flag, t15_flag;

t9_flag = false;
t10_flag = false;
t11_flag = false;
t12_flag = false;
t13_flag = false;
t14_flag = false;
t15_flag = false;

count_atoms = index->size();
if (count_atoms > max_atoms) // Control, that avoid throwing an exception.
    return(-1);
size_atoms = (size_order * 2 + 1) * (size_order * 2 + 1) * (size_order * 2 + 1);
side_pixels = (size_order * 2 + 1);

i = 0;
while (i < count_atoms) // Filling begins_atoms and ends_atoms array.
    {
    begins_atoms[i] = index->operator[](i);
    if (i + 1 < count_atoms)
        ends_atoms[i] = index->operator[](i + 1) -1;
    
    else
        ends_atoms[i] = densities->size() -1;
    i++;
    }
try
    {
    for (i = 0; i < count_atoms; i++) // Allocating memory for cubes
        {
        pointers_to_atoms[i] = new T[size_atoms];
        atoms_densities_list.push_back(pointers_to_atoms[i]);
        }
    }
catch(error_t)
    {
    Clear();
    cout << "Error with alocating memory" << endl;
    return(-1);
    }
// Multithreading code 
for (i = 0; (i + 7) < count_atoms; i = i + 8)
        { // multithreading code
        t1 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[i],
        begins_atoms[i], ends_atoms[i], size_atoms, spins, spin_paired, spin_density);
        t2 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[i + 1],
        begins_atoms[i + 1], ends_atoms[i + 1], size_atoms, spins, spin_paired, spin_density);
        t3 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[i + 2],
        begins_atoms[i + 2], ends_atoms[i + 2], size_atoms, spins, spin_paired, spin_density);
        t4 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[i + 3],
        begins_atoms[i + 3], ends_atoms[i + 3], size_atoms, spins, spin_paired, spin_density);
        t5 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[i + 4],
        begins_atoms[i + 4], ends_atoms[i + 4], size_atoms, spins, spin_paired, spin_density);
        t6 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[i + 5],
        begins_atoms[i + 5], ends_atoms[i + 5], size_atoms, spins, spin_paired, spin_density);
        t7 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[i + 6],
        begins_atoms[i + 6], ends_atoms[i + 6], size_atoms, spins, spin_paired, spin_density);
        t8 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[i + 7],
        begins_atoms[i + 7], ends_atoms[i + 7], size_atoms, spins, spin_paired, spin_density);
        t1.join();
        t2.join();
        t3.join();
        t4.join();
        t5.join();
        t6.join();
        t7.join();
        t8.join();
        }
    if (count_atoms % 8 >= 7)
        {
        t9 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[count_atoms - 7],
        begins_atoms[count_atoms - 7], ends_atoms[count_atoms - 7], size_atoms, spins, spin_paired, spin_density);
        t9_flag = true;
        }
    if (count_atoms % 8 >= 6)
        {
        t10 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[count_atoms - 6],
        begins_atoms[count_atoms - 6], ends_atoms[count_atoms - 6], size_atoms, spins, spin_paired, spin_density);
        t10_flag = true;
        }
    if (count_atoms % 8 >= 5)
        {
        t11 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[count_atoms - 5],
        begins_atoms[count_atoms - 5], ends_atoms[count_atoms - 5], size_atoms, spins, spin_paired, spin_density);
        t11_flag = true;
        }
    if (count_atoms % 8 >= 4)
        {
        t12 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[count_atoms - 4],
        begins_atoms[count_atoms - 4], ends_atoms[count_atoms - 4], size_atoms, spins, spin_paired, spin_density);
        t12_flag = true;
        }
    if (count_atoms % 8 >= 3)
        {
        t13 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[count_atoms - 3],
        begins_atoms[count_atoms - 3], ends_atoms[count_atoms - 3], size_atoms, spins, spin_paired, spin_density);
        t13_flag = true;
        }
    if (count_atoms % 8 >= 2)
        {
        t14 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[count_atoms - 2],
        begins_atoms[count_atoms - 2], ends_atoms[count_atoms - 2], size_atoms, spins, spin_paired, spin_density);
        t14_flag = true;
        }
    if (count_atoms % 8 >= 1)
        {
        t15 = thread(&Visualization_3D::Compute_density_thread, this, densities, pointers_to_atoms[count_atoms - 1],
        begins_atoms[count_atoms - 1], ends_atoms[count_atoms - 1], size_atoms, spins, spin_paired, spin_density);
        t15_flag = true;
        }
    if (t9_flag == true)
        {
        t9.join();
        t9_flag = false;
        }
    if (t10_flag == true)
        {
        t10.join();
        t10_flag = false;
        }
    if (t11_flag == true)
        {
        t11.join();
        t11_flag = false;
        }
    if (t12_flag == true)
        {
        t12.join();
        t12_flag = false;
        }
    if (t13_flag == true)
        {
        t13.join();
        t13_flag = false;
        }
    if (t14_flag == true)
        {
        t14.join();
        t14_flag = false;
        }
    if (t15_flag == true)
        {
        t15.join();
        t15_flag = false;
        }
// End of multithreading code  
return(0);  
}
template <typename T>
T Visualization_3D<T>::Compute_densities_weights_95(const vector<T*>* densities, vector<T>* weights, unsigned int size_order)
{
unsigned int i;
unsigned int densities_size;
T* densities_array[max_atoms];
T weights_array[max_atoms];
T* weights_array_p = &weights_array[0];

thread t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
bool t24_flag, t25_flag, t26_flag, t27_flag, t28_flag, t29_flag, t30_flag;

t24_flag = false;
t25_flag = false;
t26_flag = false;
t27_flag = false;
t28_flag = false;
t29_flag = false;
t30_flag = false;
densities_size = densities->size();
for (i = 0; i < densities_size; i++)
    {
    densities_array[i] = densities->operator[](i);
    }
if (OpenCL == false or not is_same<T, double>::value or weights_95.size() != densities_size)
    {
    if (weights->size() > 0)
        weights->clear();
    
    weights->reserve(densities_size);
    for (i = 0; (i + 7) < densities_size; i = i + 8)
        {    
        // multithreading code
        t16 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[i], 0.95, weights_array_p + i,
        20, size_order);
        t17 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[i + 1], 0.95,
        weights_array_p + i + 1, 20, size_order);
        t18 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[i + 2], 0.95,
        weights_array_p + i + 2, 20, size_order);
        t19 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[i + 3], 0.95,
        weights_array_p + i + 3, 20, size_order);
        t20 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[i + 4], 0.95,
        weights_array_p + i + 4, 20, size_order);
        t21 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[i + 5], 0.95,
        weights_array_p + i + 5, 20, size_order);
        t22 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[i + 6], 0.95,
        weights_array_p + i + 6, 20, size_order);
        t23 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[i + 7], 0.95,
        weights_array_p + i + 7, 20, size_order);
        t16.join();
        t17.join();
        t18.join();
        t19.join();
        t20.join();
        t21.join();
        t22.join();
        t23.join();
        }
    if (densities_size % 8 >= 7)
        {
        t24 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[densities_size - 7], 0.95,
        weights_array_p + densities_size - 7, 20, size_order);
        t24_flag = true;
        }
    if (densities_size % 8 >= 6)
        {
        t25 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[densities_size - 6], 0.95,
        weights_array_p + densities_size - 6, 20, size_order);
        t25_flag = true;
        }
    if (densities_size % 8 >= 5)
        {
        t26 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[densities_size - 5], 0.95,
        weights_array_p + densities_size - 5, 20, size_order);
        t26_flag = true;
        }
    if (densities_size % 8 >= 4)
        {
        t27 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[densities_size - 4], 0.95,
        weights_array_p + densities_size - 4, 20, size_order);
        t27_flag = true;
        }
    if (densities_size % 8 >= 3)
        {
        t28 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[densities_size - 3], 0.95,
        weights_array_p + densities_size - 3, 20, size_order);
        t28_flag = true;
        }
    if (densities_size % 8 >= 2)
        {
        t29 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[densities_size - 2], 0.95,
        weights_array_p + densities_size - 2, 20, size_order);
        t29_flag = true;
        }
    if (densities_size % 8 >= 1)
        {
        t30 = thread(&Visualization_3D::Compute_densities_weights_thread, this, densities_array[densities_size - 1], 0.95,
        weights_array_p + densities_size - 1, 20, size_order);
        t30_flag = true;
        }
    if (t24_flag == true)
        {
        t24.join();
        t24_flag = false;
        }
    if (t25_flag == true)
        {
        t25.join();
        t25_flag = false;
        }
    if (t26_flag == true)
        {
        t26.join();
        t26_flag = false;
        }
    if (t27_flag == true)
        {
        t27.join();
        t27_flag = false;
        }
    if (t28_flag == true)
        {
        t28.join();
        t28_flag = false;
        }
    if (t29_flag == true)
        {
        t29.join();
        t29_flag = false;
        }
    if (t30_flag == true)
        {
        t30.join();
        t30_flag = false;
        }
    // End of multithreading code
    for (i = 0; i < densities_size; i++)
        {
        weights->push_back(weights_array[i]);
        }
    }
return(0);
}
template <typename T>
T Visualization_3D<T>::Generate_2D_cross_section(unsigned int direction, T coordinate)
{ // Generate cross-section for picture representation
unsigned int i, j, k;
unsigned int x_pixels, y_pixels;
unsigned int count_atoms;

unsigned int side_cube;
unsigned int cube_coordinate;
unsigned int x_cross_section_min, x_cross_section_max;
unsigned int y_cross_section_min, y_cross_section_max;
unsigned int x_cross_section_diff, y_cross_section_diff;
unsigned int picture_size;
T pixel[8];
T average;
T* pointer_to_cube;
T* picture;
T weight;

count_atoms = atoms_densities_list.size(); // Initialize
side_cube = side_pixels;
average = 0;

// Check the input parameters
if (direction > 2)
    return(-1);
    
if (side_cube == 0 or cube_side_half_lenght == 0)
    return(-1);
    
switch (direction)
    {
    case 0: // y-z cross section
        if (coordinate < x_min or coordinate > x_max)
            return(-1);
        break;
    case 1: // x-z cross section
        if (coordinate < y_min or coordinate > y_max)
            return(-1);
        break;
    case 2: // x-y cross section
        if (coordinate < z_min or coordinate > z_max)
            return(-1);
        break;
    }
pixels_to_unit = T(side_pixels)/(2 * cube_side_half_lenght); // Set number of pixels to lenght unit
switch (direction) // Set dimensions of 2D cross-section in pixels
    {
    case 0: // y-z cross section {
        x_pixels = (unsigned int) ((y_max - y_min) * pixels_to_unit);
        y_pixels = (unsigned int) ((z_max - z_min) * pixels_to_unit);
        break;
        }
    case 1: // x-z cross section {
        x_pixels = (unsigned int) ((x_max - x_min) * pixels_to_unit);
        y_pixels = (unsigned int) ((z_max - z_min) * pixels_to_unit);
        break;
        }
    case 2: // x-y cross section {
        x_pixels = (unsigned int) ((x_max - x_min) * pixels_to_unit);
        y_pixels = (unsigned int) ((y_max - y_min) * pixels_to_unit);
        break;
        }
    }
try // Allocate memory for cross-section
    {
    picture_size = x_pixels * y_pixels;
    picture = new T[picture_size];
    x_pixels_list.push_back(x_pixels);
    y_pixels_list.push_back(y_pixels);
    }
catch(error_t)
    {
    cout << "Error in memory allocation" << endl;
    return(-1);
    }
// Vectorisation code
// Fill alocated memory by zeros
for (i = 0; i + 7 < picture_size; i+=8)
    {
    picture[i] = 0;
    picture[i + 1] = 0;
    picture[i + 2] = 0;
    picture[i + 3] = 0;
    picture[i + 4] = 0;
    picture[i + 5] = 0;
    picture[i + 6] = 0;
    picture[i + 7] = 0;
    }
for (i = picture_size - (picture_size % 8); i < picture_size; i++)
    picture[i] = 0;
// End of vectorisation code

for (i = 0; i < count_atoms; i++) 
    {
    if (weights_95.size() > i)
        weight = weights_95[i];
    else
        weight = 0;
    // Check for possibility to adding atoms
    if (direction == 0 and (coordinates[i][0] < x_min or coordinates[i][0] > x_max)) // y-z cross section
        continue;
    if (direction == 1 and (coordinates[i][1] < y_min or coordinates[i][1] > y_max)) // x-z cross section
        continue;
    if (direction == 2 and (coordinates[i][2] < z_min or coordinates[i][2] > z_max)) // x-y cross section
        continue;
        
    // Preparing variables for adding atoms to cross-section
    cube_coordinate = side_cube/2 + (unsigned int) (coordinate - coordinates[i][direction] * pixels_to_unit);
    if (cube_coordinate > side_cube) // Check for overflow
        cube_coordinate = side_cube;
        
    pointer_to_cube = atoms_densities_list[i];
    
    if (direction == 0) // y-z cross section
        {
        x_cross_section_min = (unsigned int) ((coordinates[i][1] - cube_side_half_lenght - y_min) * pixels_to_unit);
        if (x_cross_section_min < y_min * pixels_to_unit) // Check for overflow
            x_cross_section_min = y_min * pixels_to_unit;
        
        y_cross_section_min = (unsigned int) ((coordinates[i][2] - cube_side_half_lenght - z_min) * pixels_to_unit);
        if (y_cross_section_min < z_min * pixels_to_unit) // Check for overflow
            y_cross_section_min = z_min * pixels_to_unit;
        }
    if (direction == 1) // x-z cross section
        {
        x_cross_section_min = (unsigned int) ((coordinates[i][0] - cube_side_half_lenght - x_min) * pixels_to_unit);
        if (x_cross_section_min < 0) // Check for overflow
            x_cross_section_min = 0;
        
        y_cross_section_min = (unsigned int) ((coordinates[i][2] - cube_side_half_lenght - z_min) * pixels_to_unit);
        if (y_cross_section_min < 0) // Check for overflow
            y_cross_section_min = 0;
        }
    if (direction == 2) // x-y cross section
        {
        x_cross_section_min = (unsigned int) ((coordinates[i][0] - cube_side_half_lenght - x_min) * pixels_to_unit);
        if (x_cross_section_min < 0) // Check for overflow
            x_cross_section_min = 0;
        
        y_cross_section_min = (unsigned int) ((coordinates[i][1] - cube_side_half_lenght - y_min) * pixels_to_unit);
        if (y_cross_section_min < 0) // Check for overflow
            y_cross_section_min = 0;
        }
    x_cross_section_max = x_cross_section_min + side_pixels;
    if (x_cross_section_max > x_pixels) // Check for overflow
        x_cross_section_max = x_pixels;
    
    y_cross_section_max = y_cross_section_min + side_pixels;
    if (y_cross_section_max > y_max * pixels_to_unit) // Check for overflow
        y_cross_section_max = y_pixels;
    
    x_cross_section_diff = x_cross_section_max - x_cross_section_min;
    y_cross_section_diff = y_cross_section_max - y_cross_section_min;
    
    // Adding atoms to cross-section
    if (direction == 0) // y-z cross section
        {
        for (j = 0; j < x_cross_section_diff; j++)
            for (k = 0; k < y_cross_section_diff; k++)
                {
                pixel[0] = pointer_to_cube[cube_coordinate + (j * side_cube) + (k * side_cube * side_cube)];
                if (abs(pixel[0]) >= weight)
                    {
                    picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[0];
                    average += pixel[0];
                    }
                }
        }
    if (direction == 1) // x-z cross section
        {
        for (j = 0; j < x_cross_section_diff; j++)
            for (k = 0; k < y_cross_section_diff; k++)
                {
                pixel[0] = pointer_to_cube[j + (cube_coordinate * side_cube) + (k * side_cube * side_cube)];
                if (abs(pixel[0]) >= weight)
                    {
                    picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[0];
                    average += pixel[0];
                    }
                }
        }
    // Vectorisation code
    if (direction == 2) // x-y cross section
        {
        for (j = 0; j < x_cross_section_diff; j++)
            {
            for (k = 0; k + 7 < y_cross_section_diff; k+=8)
                {
                pixel[0] = pointer_to_cube[j + (k * side_cube) + (cube_coordinate * side_cube * side_cube)];
                pixel[1] = pointer_to_cube[j + (k * side_cube) + (cube_coordinate * side_cube * side_cube) + 1];
                pixel[2] = pointer_to_cube[j + (k * side_cube) + (cube_coordinate * side_cube * side_cube) + 2];
                pixel[3] = pointer_to_cube[j + (k * side_cube) + (cube_coordinate * side_cube * side_cube) + 3];
                pixel[4] = pointer_to_cube[j + (k * side_cube) + (cube_coordinate * side_cube * side_cube) + 4];
                pixel[5] = pointer_to_cube[j + (k * side_cube) + (cube_coordinate * side_cube * side_cube) + 5];
                pixel[6] = pointer_to_cube[j + (k * side_cube) + (cube_coordinate * side_cube * side_cube) + 6];
                pixel[7] = pointer_to_cube[j + (k * side_cube) + (cube_coordinate * side_cube * side_cube) + 7];
                if (weight > 0)
                    {
                    if (abs(pixel[0]) >= weight)
                        {
                        picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[0];
                        average += pixel[0];
                        }
                    if (abs(pixel[1]) >= weight)
                        {
                        picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[1];
                        average += pixel[1];
                        }
                    if (abs(pixel[2]) >= weight)
                        {
                        picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[2];
                        average += pixel[2];
                        }
                    if (abs(pixel[3]) >= weight)
                        {
                        picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[3];
                        average += pixel[3];
                        }
                    if (abs(pixel[4]) >= weight)
                        {
                        picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[4];
                        average += pixel[4];
                        }
                    if (abs(pixel[5]) >= weight)
                        {
                        picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[5];
                        average += pixel[5];
                        }
                    if (abs(pixel[6]) >= weight)
                        {
                        picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[6];
                        average += pixel[6];
                        }
                    if (abs(pixel[7]) >= weight)
                        {
                        picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[7];
                        average += pixel[7];
                        }
                    }
                else
                    {
                    picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[0];
                    picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels) + 1] += pixel[1];
                    picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels) + 2] += pixel[2];
                    picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels) + 3] += pixel[3];
                    picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels) + 4] += pixel[4];
                    picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels) + 5] += pixel[5];
                    picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels) + 6] += pixel[6];
                    picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels) + 7] += pixel[7];
                    average += pixel[0];
                    average += pixel[1];
                    average += pixel[2];
                    average += pixel[3];
                    average += pixel[4];
                    average += pixel[5];
                    average += pixel[6];
                    average += pixel[7];
                    }
                }
            for (k = y_cross_section_diff - (y_cross_section_diff % 8); k < y_cross_section_diff; k++)
                {
                pixel[0] = pointer_to_cube[j + (k * side_cube) + (cube_coordinate * side_cube * side_cube)];
                if (abs(pixel[0]) >= weight)
                    {
                    picture[(x_cross_section_min + j) + ((y_cross_section_min + k) * x_pixels)] += pixel[0];
                    average += pixel[0];
                    }
                }
            }
        } // End of vectorisation code
    }
average = average/T(picture_size);
averages.push_back(average);
Cross_sections_2D.push_back(picture);
return(0);
}
template <typename T>
T Visualization_3D<T>::Clear()
{  // Delete pointers in atoms_list, clear class vectors and set default coordinates.
unsigned int i;
unsigned int size;

size = atoms_densities_list.size();
for (i = 0; i < size; i++)
    delete[] atoms_densities_list[i];
    
size = Cross_sections_2D.size();
for (i = 0; i < size; i++)
    delete[] Cross_sections_2D[i];

atoms_densities_list.clear();
coordinates.clear();
x_pixels_list.clear();
y_pixels_list.clear();
averages.clear();
weights_95.clear();
cube_side_half_lenght = 0;
side_pixels = 0;
pixels_to_unit = 0;
x_min = 0;
y_min = 0;
z_min = 0;
x_max = 0;
y_max = 0;
z_max = 0;
return(0);
}
template <typename T> // Destructor of class instance
Visualization_3D<T>::~Visualization_3D() {
Clear();}
#endif // VISUALIZATION_3D_H
/*
Author of this source code Ing. Pavel Florian Ph.D. licensed this source code under the the Apache License:

Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/
*/
