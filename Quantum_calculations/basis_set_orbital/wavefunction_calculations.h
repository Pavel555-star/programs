#include <iostream>
#include <string>
#include <cmath>
#include <thread> // compiler parameters: -pthread -ffast-math -fno-finite-math-only
#include <vector>
#include <bits/stdc++.h>
#include "basis_set_calculations.h"
using namespace std;

template <typename T>
class Wavefunction_calculations: public basis_set_calculations<T>
{
public:
    unsigned int count_atoms = 0;
    vector <unsigned int> electron_to_atom_numbers;
    vector <unsigned int> index_atoms;
    vector<T>* energies = nullptr;
    T x_center_of_mass = 0;
    T y_center_of_mass = 0;
    T z_center_of_mass = 0;
    
    struct dipole_moment {
    vector <T> partial_charges;
    vector <T> x;
    vector <T> y;
    vector <T> z;
    vector <T> x_dipole_moments;
    vector <T> y_dipole_moments;
    vector <T> z_dipole_moments;
    T sum_dipole_moment = 0; // in Bohr radius * elementary charge
    T sum_dipole_moment_Debye = 0; // in Debye
    } dipoles;
    
    struct central_cations {
    vector <unsigned int> atom_numbers;
    vector <unsigned int> configuration_numbers;
    vector <unsigned int> high_spin;
    } central_cation_atoms;
    
    struct symetry_axes {
    vector <T> u_x; // x = u_x * t + x_center
    vector <T> u_y; // y = u_y * t + y_center
    vector <T> u_z; // z = u_z * t + z_center
    T x_center = 0; // coordinates of mass center
    T y_center = 0;
    T z_center = 0;
    bool symetry_center = false;
    bool linear = false;
    } symetry_axes_parameters;
    
    struct symetry_planes {
    vector <T> a; // a * x + b * y + c * z + d = 0
    vector <T> b;
    vector <T> c;
    vector <T> d;
    bool planar = false;
    } symetry_planes_parameters;
    
private:
    T Create_index_atoms();
    T Detect_center_of_mass();
    T Detect_dipole_moments(dipole_moment *dipole_moment);
    T Create_crystal_field(central_cations *central_cations);
public:
    T String_to_advanced_parameters(string UI_input, unsigned int size_order,
    bool extern_coordinates, vector<T>* x_2, vector<T>* y_2, vector<T>* z_2);
    T Execute_calculation(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool dealocate,
    vector<T>* values, vector<T>* spin_density_vector,  vector<T>* spin_values);
    T Calculate_Huckel_Matrix(T* Huckel_matrix, unsigned int* Huckel_matrix_order, vector<unsigned int> atom_numbers);
    T Detect_symetry_information(symetry_axes *symetry_axes, symetry_planes *symetry_planes);
    T D3_vector_multiply(T* a, T* b, T* c);
    T Clear();
    ~Wavefunction_calculations();
};
template <typename T>
T Wavefunction_calculations<T>::Create_index_atoms()
    {
    unsigned int i;
    unsigned int sum_electrons;
    unsigned int count_bonds;
    T* x;
    T* y;
    T* z;
    x = this->results.x.data();
    y = this->results.y.data();
    z = this->results.z.data();
    sum_electrons = this->results.n.size();
    i = 0;
    count_atoms = 1;
    
    index_atoms.push_back(0);
    electron_to_atom_numbers.push_back(count_atoms);
    
    for (i = 1; i < sum_electrons; i++)
        {
        if (x[i] != x[i-1] or y[i] != y[i-1] or z[i] != z[i-1])
            {
            index_atoms.push_back(i);
            count_atoms++;
            }
        electron_to_atom_numbers.push_back(count_atoms);
        }
    return(0);
    }
template <typename T>
T Wavefunction_calculations<T>::Detect_center_of_mass()
    {
    unsigned int i;
    unsigned int Z;
    unsigned int sum_Z;
    T sum_x;
    T sum_y;
    T sum_z;
    
    Z = 0;
    sum_Z = 0;
    sum_x = 0;
    sum_y = 0;
    sum_z = 0;
    
    for (i = 0; i < count_atoms; i++)
        {
        Z = this->results.Z[index_atoms[i]];
        sum_Z = sum_Z + Z;
        sum_x = sum_x + (this->results.x[index_atoms[i]] * Z);
        sum_y = sum_y + (this->results.y[index_atoms[i]] * Z);
        sum_z = sum_z + (this->results.z[index_atoms[i]] * Z);
        }
    x_center_of_mass = sum_x/sum_Z;
    y_center_of_mass = sum_y/sum_Z;
    z_center_of_mass = sum_z/sum_Z;
    return(0);
    }
template <typename T>
T Wavefunction_calculations<T>::Detect_dipole_moments(dipole_moment *dipole_moment)
    {
    unsigned int i, j;
    unsigned int count_electrons;
    unsigned int atom_begin, next_atom_begin;
    
    T charge;
    T ion_charge;
    T wavefunction_coefficient;
    T x, y, z;
    T x_dipole_moment, y_dipole_moment, z_dipole_moment;
    T sum_x_dipole_moment, sum_y_dipole_moment, sum_z_dipole_moment;
    
    charge = 0;
    sum_x_dipole_moment = 0;
    sum_y_dipole_moment = 0;
    sum_z_dipole_moment = 0;
    count_electrons = electron_to_atom_numbers.size();
    
    for (i = 0; i < count_atoms; i++)
        {
        atom_begin = index_atoms[i];
        if ((i + 1) < count_atoms)
            next_atom_begin = index_atoms[i + 1];
        else
            next_atom_begin = count_electrons;
        
        charge = 0;
        ion_charge = 0;
        for (j = atom_begin; j < next_atom_begin; j++)
            {
            wavefunction_coefficient = this->results.wavefunction_coefficients[j]; // including polarization of bonds
            if (wavefunction_coefficient != 1)
                charge = charge - (wavefunction_coefficient * wavefunction_coefficient - 1.00);
                
            if (ion_charge == 0 and this->results.reduced_Z[j] != this->results.count_electrons[j])
                ion_charge = this->results.reduced_Z[j] - this->results.count_electrons[j];
            }
        charge = charge + ion_charge;
        if (charge != 0)
            {
            x = this->results.x[index_atoms[i]] - x_center_of_mass;
            y = this->results.y[index_atoms[i]] - y_center_of_mass;
            z = this->results.z[index_atoms[i]] - z_center_of_mass;
            x_dipole_moment = x * charge;
            y_dipole_moment = y * charge;
            z_dipole_moment = z * charge;
            sum_x_dipole_moment = sum_x_dipole_moment + x_dipole_moment;
            sum_y_dipole_moment = sum_y_dipole_moment + y_dipole_moment;
            sum_z_dipole_moment = sum_z_dipole_moment + z_dipole_moment;
            
            dipole_moment->partial_charges.push_back(charge);
            dipole_moment->x.push_back(x);
            dipole_moment->y.push_back(y);
            dipole_moment->z.push_back(z);
            dipole_moment->x_dipole_moments.push_back(x_dipole_moment);
            dipole_moment->y_dipole_moments.push_back(y_dipole_moment);
            dipole_moment->z_dipole_moments.push_back(z_dipole_moment);
            }
        }
    dipole_moment->sum_dipole_moment = sqrt((sum_x_dipole_moment * sum_x_dipole_moment)
    + (sum_y_dipole_moment * sum_y_dipole_moment) + (sum_z_dipole_moment * sum_z_dipole_moment));
    dipole_moment->sum_dipole_moment_Debye = dipole_moment->sum_dipole_moment * 2.54;
    return(0);
    }
template <typename T>
T Wavefunction_calculations<T>::Create_crystal_field(central_cations *central_cations)
    {
    unsigned int i, j, k;
    unsigned int atom_begin, next_atom_begin;
    unsigned int size;
    unsigned int* n;
    unsigned int* l;
    int* m;
    
    T atom_numbers[this->max_electrons];
    T configuration_numbers[this->max_electrons];
    T* s;
    T d_energy_min, d_energy_max, d_energy_difference;
    bool high_spin;
    
    int d_order_1_h[10] = {-2, 0, 2, -1, 1, -2, 0, 2, -1, 1}; // Octahedral shape
    int d_order_1_l[10] = {-2, -2, 0, 0, 2, 2, -1, -1, 1, 1};
    int d_order_2_h[10] = {-2, 0, 2, -1, 1, -2, 0, 2, -1, 1}; // Pentagonal bipyramidal shape
    int d_order_2_l[10] = {-2, -2, 0, 0, 2, 2, -1, -1, 1, 1};
    int d_order_3_h[10] = {1, 2, -1, -2, 0, 1, 2, -1, -2, 0}; // Square antiprismatic shape
    int d_order_3_l[10] = {1, 1, 2, 2, -1, -1, -2, -2, 0, 0};
    int d_order_4_h[10] = {-2, 0, 1, 2, -1, -2, 0, 1, 2, -1}; // Square planar shape
    int d_order_4_l[10] = {-2, -2, 0, 0, 1, 1, 2, 2, -1, -1};
    int d_order_5_h[10] = {-2, 0, 2, 1, -1, -2, 0, 2, 1, -1}; // Square pyramidal shape
    int d_order_5_l[10] = {-2, -2, 0, 0, 2, 2, 1, 1, -1, -1};
    int d_order_6_h[10] = {-1, 1, -2, 0, 2, -1, 1, -2, 0, 2}; // Tetrahedral shape
    int d_order_6_l[10] = {-1, -1, 1, 1, -2, -2, 0, 0, 2, 2};
    int d_order_7_h[10] = {-2, 0, 2, -1, 1, -2, 0, 2, -1, 1}; // Trigonal bipyramidal shape
    int d_order_7_l[10] = {-2, -2, 0, 0, 2, 2, -1, -1, 1, 1};
    T d_s_order_h[10] = {0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5, -0.5};
    T d_s_order_l[10] = {0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5};
    
    int* pointer_to_m_order;
    T* pointer_to_s_order;
    size = central_cations->atom_numbers.size();
    for (i = 0; i < size; i++) // Copying vectors to stack
        {
        atom_numbers[i] = central_cations->atom_numbers[i];
        configuration_numbers[i] = central_cations->configuration_numbers[i];
        }
    n = this->results.n.data();
    l = this->results.l.data();
    m = this->results.m.data();
    s = this->results.spins.data();
        
    for (i = 0; i < size; i++) // For every central cation
        {
        if (central_cations->atom_numbers[i] > 0)
            {
            atom_begin = index_atoms[central_cations->atom_numbers[i] -1]; // searching atom range in this->results
            if ((i + 1) < count_atoms)
                next_atom_begin = index_atoms[central_cations->atom_numbers[i + 1] -1];
            else
                next_atom_begin = this->results.n.size();
            }
        if (central_cations->high_spin[i] == 0)
            high_spin = false;
        else
            high_spin = true;
            
        if (configuration_numbers[i] == 0)
            continue; // 0 - no crystal field theory corrections
                
        if (configuration_numbers[i] == 1 and high_spin == true) // Octahedral shape
            {
            pointer_to_m_order = &d_order_1_h[0];
            pointer_to_s_order = &d_s_order_h[0];
            }
        if (configuration_numbers[i] == 1 and high_spin == false) // Octahedral shape
            {
            pointer_to_m_order = &d_order_1_l[0];
            pointer_to_s_order = &d_s_order_l[0];
            }
        if (configuration_numbers[i] == 2 and high_spin == true) // Pentagonal bipyramidal shape
            {
            pointer_to_m_order = &d_order_2_h[0];
            pointer_to_s_order = &d_s_order_h[0];
            }
        if (configuration_numbers[i] == 2 and high_spin == false) // Pentagonal bipyramidal shape
            {
            pointer_to_m_order = &d_order_2_l[0];
            pointer_to_s_order = &d_s_order_l[0];
            }
        if (configuration_numbers[i] == 3 and high_spin == true) // Square antiprismatic shape
            {
            pointer_to_m_order = &d_order_3_h[0];
            pointer_to_s_order = &d_s_order_h[0];
            }
        if (configuration_numbers[i] == 3 and high_spin == false) // Square antiprismatic shape
            {
            pointer_to_m_order = &d_order_3_l[0];
            pointer_to_s_order = &d_s_order_l[0];
            }
        if (configuration_numbers[i] == 4 and high_spin == true) // Square planar shape
            {
            pointer_to_m_order = &d_order_4_h[0];
            pointer_to_s_order = &d_s_order_h[0];
            }
        if (configuration_numbers[i] == 4 and high_spin == false) // Square planar shape
            {
            pointer_to_m_order = &d_order_4_l[0];
            pointer_to_s_order = &d_s_order_l[0];
            }
        if (configuration_numbers[i] == 5 and high_spin == true) // Square pyramidal shape
            {
            pointer_to_m_order = &d_order_5_h[0];
            pointer_to_s_order = &d_s_order_h[0];
            }
        if (configuration_numbers[i] == 5 and high_spin == false) // Square pyramidal shape
            {
            pointer_to_m_order = &d_order_5_l[0];
            pointer_to_s_order = &d_s_order_l[0];
            }
        if (configuration_numbers[i] == 6 and high_spin == true) // Tetrahedral shape
            {
            pointer_to_m_order = &d_order_6_h[0];
            pointer_to_s_order = &d_s_order_h[0];
            }
        if (configuration_numbers[i] == 6 and high_spin == false) // Tetrahedral shape
            {
            pointer_to_m_order = &d_order_6_l[0];
            pointer_to_s_order = &d_s_order_l[0];
            }
        if (configuration_numbers[i] == 7 and high_spin == true) // Trigonal bipyramidal shape
            {
            pointer_to_m_order = &d_order_7_h[0];
            pointer_to_s_order = &d_s_order_h[0];
            }
        if (configuration_numbers[i] == 7 and high_spin == false) // Trigonal bipyramidal shape
            {
            pointer_to_m_order = &d_order_7_l[0];
            pointer_to_s_order = &d_s_order_l[0];
            }
        k = 0; // Applying m and s acording to model
        for (j = atom_begin; j < next_atom_begin; j++) // Copying vectors to stack
            {
            if (l[j] == 2)
                {
                m[j] = pointer_to_m_order[k];
                s[j] = pointer_to_s_order[k];
                k++;
                }
            }
        }
    return(0);
    }
template <typename T>
T Wavefunction_calculations<T>::String_to_advanced_parameters(string UI_input, unsigned int size_order,
    bool extern_coordinates, vector<T>* x_2, vector<T>* y_2, vector<T>* z_2)
    { // kompletně přepsat
    unsigned int i, j, k;
    unsigned int last_readed_index;
    unsigned int input_size;
    unsigned int count_excitations;
    unsigned int count_hybridization;
    unsigned int count_pi_bonding;
    unsigned int count_central_cations;
    
    bool read_switch;
    string input;
    string character;
    string central_cations_string[this->max_atoms * 3];
    vector<unsigned int> x;
    vector<T> y;
    vector<T> z;
    vector<T> potentials;
    const string numbers = "+-0123456789.";
    const string begin_brackets = "<";
    const string end_brackets = ">";
    
    read_switch = true;
    character = " ";
    input_size = UI_input.size();
    input = UI_input;
    count_excitations = 0;
    count_pi_bonding = 0;
    count_hybridization = 0;
    count_central_cations = 0;
    
    if (this->String_to_list_electrons(UI_input, size_order, extern_coordinates, x_2, y_2, z_2) == -1)
        return(-1); // generate list of electons and check of first input string by function from parrent class
     
    read_switch = false;
    last_readed_index = 10000;
    for (i = 0; i < input_size; i++) // reading list of central cations
        {
        for (j = 0; j < numbers.size(); j++) // check for input[i] in list of numbers
            if (input[i] == numbers[j] and read_switch == true)
                {
                if ((i - last_readed_index) != 1)
                    count_central_cations++;
                        
                character[0] = input[i];
                central_cations_string[count_central_cations - 1].append(character);
                last_readed_index = i;
                }
        if (input[i] == begin_brackets[0])
            read_switch = true;
        if (input[i] == end_brackets[0])
            read_switch = false;
        }
    for (i = 0; i < count_central_cations; i++)
        {
        if (i % 3 == 0)
            central_cation_atoms.atom_numbers.push_back(stoi(central_cations_string[i]));
        if (i % 3 == 1)
            central_cation_atoms.configuration_numbers.push_back(stoi(central_cations_string[i]));
        if (i % 3 == 2)
            central_cation_atoms.high_spin.push_back(stoi(central_cations_string[i]));
        }
    return(0);
    }
template <typename T>
T Wavefunction_calculations<T>::Execute_calculation(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order,
bool dealocate, vector<T>* values, vector<T>* spin_density_vector,  vector<T>* spin_values)
    {
    T Hamiltonian;
    bool repetition_flag;
    
    energies = values;
    if (Create_index_atoms() == -1)
       return(-1);
    if (Detect_center_of_mass() == -1)
       return(-1);
    if (Detect_dipole_moments(&dipoles) == -1)
       return(-1);
    if (central_cation_atoms.atom_numbers.size() > 0)
        Create_crystal_field(&central_cation_atoms);
       
    Hamiltonian = this->Calculate(max_iterations, minimal_fidelity, size_order, dealocate, values, spin_density_vector,
    spin_values);
    return(Hamiltonian);
    }
template <typename T>
T Wavefunction_calculations<T>::Calculate_Huckel_Matrix(T* Huckel_matrix, unsigned int* Huckel_matrix_order,
vector<unsigned int> atom_numbers)
    {
    unsigned int i, j;
    unsigned int count_pi_electrons;
    unsigned int index_pi_electrons[this->max_electrons];
    int* antibonding = this->results.antibonding.data();
    int* pi_bonding = this->results.pi_bonding.data();
    unsigned int n_down_electrons[7] = {0, 2, 10, 18, 36, 54, 86};
    unsigned int l_up_n[4] = {0, 0, 1, 2};
    unsigned int matrix_order = this->results.n.size();
    
    T layer_coefficient;
    T ionization_energy;
    T alpha, beta;
    T average_diagonal;
    
    alpha = -11.4 * this->e;
    beta = 0.52 * this->e;
    count_pi_electrons = 0;
    average_diagonal = 0;
    atom_numbers.clear();
    
    for (i = 0; i < matrix_order; i++) // Creating list of pi electrons for creating Huckel matrix
        {
        if (pi_bonding[i] >= 0 and antibonding[i] < 0)
            {
            index_pi_electrons[count_pi_electrons] = i;
            count_pi_electrons++;
            }
        }
    Huckel_matrix_order[0] = count_pi_electrons;
   
    for (i = 0; i < count_pi_electrons * count_pi_electrons; i++)
        Huckel_matrix[i] = 0;
        
    for (i = 0; i < count_pi_electrons; i++) // Copying resonance integrals
        {
        for (j = 0; j < count_pi_electrons; j++)
            {
            if (j != i and electron_to_atom_numbers[index_pi_electrons[i]] != electron_to_atom_numbers[index_pi_electrons[j]])
                Huckel_matrix[(i * count_pi_electrons) + j] = Huckel_matrix[(i * count_pi_electrons) + j]
                + abs(this->resonance_integral_matrix[(index_pi_electrons[i] * matrix_order) +  index_pi_electrons[j]]);
            }
        }
    for (i = 0; i < count_pi_electrons; i++) // Copying nuclear_atraction and coulombic integrals
        {
        ionization_energy = this->nuclear_atraction_integral_matrix[index_pi_electrons[i] * (matrix_order + 1)]; 
        for (j = 0; j < count_pi_electrons; j++)
            ionization_energy = ionization_energy
            + this->coulombic_integral_matrix[(index_pi_electrons[i] * matrix_order) + index_pi_electrons[j]];
            
        ionization_energy = ionization_energy;
        Huckel_matrix[i * (count_pi_electrons + 1)] = Huckel_matrix[i * (count_pi_electrons + 1)] + ionization_energy;
        }
    for (i = 0; i < count_pi_electrons; i++) // diagonal transition to shape (alpha - alpha max) / beta
        {
        Huckel_matrix[i * (count_pi_electrons + 1)] = Huckel_matrix[i * (count_pi_electrons + 1)] - alpha;
        if (beta != 0) 
            Huckel_matrix[i * (count_pi_electrons + 1)] = Huckel_matrix[i * (count_pi_electrons + 1)]/beta;
        }
    
    for (i = 0; i < count_pi_electrons; i++) // finding average of diagonal
        average_diagonal = average_diagonal + Huckel_matrix[i * (count_pi_electrons + 1)];
    
    average_diagonal = average_diagonal/count_pi_electrons; 
    for (i = 0; i < count_pi_electrons; i++) // diagonal shift
        Huckel_matrix[i * (count_pi_electrons + 1)] = Huckel_matrix[i * (count_pi_electrons + 1)] - average_diagonal;
            
    if (beta != 0) // non-diagonal transition to shape beta/beta max
        for (i = 0; i < count_pi_electrons; i++)
            for (j = 0; j < count_pi_electrons; j++)
                if (i != j)
                    Huckel_matrix[i + (j * count_pi_electrons)] = Huckel_matrix[i + (j * count_pi_electrons)]/beta;
    
    for (i = 0; i < count_pi_electrons; i++) // Copying atom_numbers to vector
        atom_numbers.push_back(electron_to_atom_numbers[index_pi_electrons[i]]);
        
    return(0);
    }
template <typename T>
T Wavefunction_calculations<T>::Detect_symetry_information(symetry_axes *symetry_axes, symetry_planes *symetry_planes)
    {
    unsigned int i, j, k, l, m;
    T x[this->max_atoms];
    T y[this->max_atoms];
    T z[this->max_atoms];
    T Z[this->max_atoms];
    T x_center; // coordinates of mass center
    T y_center;
    T z_center;
    // axes
    T u_x_0[this->max_atoms];
    T u_y_0[this->max_atoms];
    T u_z_0[this->max_atoms];
    T u_x[this->max_atoms]; // (x + y + z) * t + T = 0
    T u_y[this->max_atoms];
    T u_z[this->max_atoms];
    vector<T> u_x_2(this->max_atoms * this->max_atoms);
    vector<T> u_y_2(this->max_atoms * this->max_atoms);
    vector<T> u_z_2(this->max_atoms * this->max_atoms);
    unsigned int count_axes = 0;
    // planes
    T a[this->max_atoms]; // a * x + b * y + c * z + d = 0
    T b[this->max_atoms];
    T c[this->max_atoms];
    T d[this->max_atoms];
    // mirroring points
    T x_mirroring[this->max_atoms];
    T y_mirroring[this->max_atoms];
    T z_mirroring[this->max_atoms];
    T t[this->max_atoms];
    T deviation = 0.01;
    bool center_symetry;
    bool axe_symetry, axe_symetry_flag;
    // planes arrays
    vector<T> a_2(this->max_atoms * this->max_atoms);
    vector<T> b_2(this->max_atoms * this->max_atoms);
    vector<T> c_2(this->max_atoms * this->max_atoms);
    vector<T> d_2(this->max_atoms * this->max_atoms);
    
    unsigned int count_planes = 0;
    bool plane_symetry, plane_symetry_flag;
    bool previous_finded;
    T ratio;
    // for detecting symetry information from location of atoms
    T u_A[3];
    T u_B[3];
    T u_C[3];
    T u_D[3];
    T u_A_B[3];
    T u_A_C[3];
    T planar_equation[4];
    
    for (i = 0; i < count_atoms; i++)
        { // Copy coordinates and proton numbers to stack
        x[i] = this->results.x[index_atoms[i]];
        y[i] = this->results.y[index_atoms[i]];
        z[i] = this->results.z[index_atoms[i]];
        Z[i] = this->results.Z[index_atoms[i]];
        }
    x_center = x_center_of_mass; // For center of mass and lines from center of mass are center parameters center of mass
    y_center = y_center_of_mass;
    z_center = z_center_of_mass;
    
    // Here begins detecting symetry information from location of atoms and mass center of system
    for (i = 0; i < count_atoms; i++) // Fill axes paramaters between center and atoms -differenced between atoms and center of mass
        {
        u_x_0[i] = x[i] - x_center;
        u_y_0[i] = y[i] - y_center;
        u_z_0[i] = z[i] - z_center;
        }
    // search for symetry center
    center_symetry = true;
    for (i = 0; i + 7 < count_atoms; i+=8) // Compute planes between lines and points
        {
        // computing mirroring points
        x_mirroring[i] = - x[i] + 2 * (x_center);
        x_mirroring[i + 1] = - x[i + 1] + 2 * (x_center);
        x_mirroring[i + 2] = - x[i + 2] + 2 * (x_center);
        x_mirroring[i + 3] = - x[i + 3] + 2 * (x_center);
        x_mirroring[i + 4] = - x[i + 4] + 2 * (x_center);
        x_mirroring[i + 5] = - x[i + 5] + 2 * (x_center);
        x_mirroring[i + 6] = - x[i + 6] + 2 * (x_center);
        x_mirroring[i + 7] = - x[i + 7] + 2 * (x_center);
        y_mirroring[i] = - y[i] + 2 * (y_center);
        y_mirroring[i + 1] = - y[i + 1] + 2 * (y_center);
        y_mirroring[i + 2] = - y[i + 2] + 2 * (y_center);
        y_mirroring[i + 3] = - y[i + 3] + 2 * (y_center);
        y_mirroring[i + 4] = - y[i + 4] + 2 * (y_center);
        y_mirroring[i + 5] = - y[i + 5] + 2 * (y_center);
        y_mirroring[i + 6] = - y[i + 6] + 2 * (y_center);
        y_mirroring[i + 7] = - y[i + 7] + 2 * (y_center);
        z_mirroring[i] = - z[i] + 2 * (z_center);
        z_mirroring[i + 1] = - z[i + 1] + 2 * (z_center);
        z_mirroring[i + 2] = - z[i + 2] + 2 * (z_center);
        z_mirroring[i + 3] = - z[i + 3] + 2 * (z_center);
        z_mirroring[i + 4] = - z[i + 4] + 2 * (z_center);
        z_mirroring[i + 5] = - z[i + 5] + 2 * (z_center);
        z_mirroring[i + 6] = - z[i + 6] + 2 * (z_center);
        z_mirroring[i + 7] = - z[i + 7] + 2 * (z_center);
            
        for (j = i; j <  i + 7; j++) // search mirroring points in other points
            {
            axe_symetry_flag = false;
            for (k = 0; k < count_atoms; k++)
                {
                if (x[k] < x_mirroring[j] + deviation  and x[k] > x_mirroring[j] - deviation
                and y[k] < y_mirroring[j] + deviation  and y[k] > y_mirroring[j] - deviation
                and z[k] < z_mirroring[j] + deviation  and z[k] > z_mirroring[j] - deviation
                and Z[k] == Z[j])
                    {
                    axe_symetry_flag = true; // found
                    break;
                    }
                }
            if (axe_symetry_flag == false)
                {
                center_symetry = false;
                break;
                }
            }
        if (axe_symetry == false)
           break;
        }
    for (i = (count_atoms - count_atoms % 8); i < count_atoms; i++) // Compute planes between lines and points
        {
        x_mirroring[i] = - x[i] + 2 * (x_center);
        y_mirroring[i] = - y[i] + 2 * (y_center);
        z_mirroring[i] = - z[i] + 2 * (z_center);
            
        for (j = i; j <  i + 1; j++)
            {
            axe_symetry_flag = false;
            for (k = 0; k < count_atoms; k++)
                {
                if (x[k] < x_mirroring[j] + deviation  and x[k] > x_mirroring[j] - deviation
                and y[k] < y_mirroring[j] + deviation  and y[k] > y_mirroring[j] - deviation
                and z[k] < z_mirroring[j] + deviation  and z[k] > z_mirroring[j] - deviation
                and Z[k] == Z[j])
                    {
                    axe_symetry_flag = true; // found
                    break;
                    }
                }
            if (axe_symetry_flag == false)
                {
                center_symetry = false;
                break;
                }
            }
        if (center_symetry == false)
           break;
        }
    for (i = 0; i < count_atoms; i++) // search for symetry axes
        {
        for (j = i; j + 7 < count_atoms; j++)
            {
            u_x[j] = (u_x_0[j] + u_x_0[i])/2;
            u_x[j + 1] = (u_x_0[j + 1] + u_x_0[i])/2;
            u_x[j + 2] = (u_x_0[j + 2] + u_x_0[i])/2;
            u_x[j + 3] = (u_x_0[j + 3] + u_x_0[i])/2;
            u_x[j + 4] = (u_x_0[j + 4] + u_x_0[i])/2;
            u_x[j + 5] = (u_x_0[j + 5] + u_x_0[i])/2;
            u_x[j + 6] = (u_x_0[j + 6] + u_x_0[i])/2;
            u_x[j + 7] = (u_x_0[j + 7] + u_x_0[i])/2;
            u_y[j] = (u_y_0[j] + u_y_0[i])/2;
            u_y[j + 1] = (u_y_0[j + 1] + u_y_0[i])/2;
            u_y[j + 2] = (u_y_0[j + 2] + u_y_0[i])/2;
            u_y[j + 3] = (u_y_0[j + 3] + u_y_0[i])/2;
            u_y[j + 4] = (u_y_0[j + 4] + u_y_0[i])/2;
            u_y[j + 5] = (u_y_0[j + 5] + u_y_0[i])/2;
            u_y[j + 6] = (u_y_0[j + 6] + u_y_0[i])/2;
            u_y[j + 7] = (u_y_0[j + 7] + u_y_0[i])/2;
            u_z[j] = (u_z_0[j] + u_z_0[i])/2;
            u_z[j + 1] = (u_z_0[j + 1] + u_z_0[i])/2;
            u_z[j + 2] = (u_z_0[j + 2] + u_z_0[i])/2;
            u_z[j + 3] = (u_z_0[j + 3] + u_z_0[i])/2;
            u_z[j + 4] = (u_z_0[j + 4] + u_z_0[i])/2;
            u_z[j + 5] = (u_z_0[j + 5] + u_z_0[i])/2;
            u_z[j + 6] = (u_z_0[j + 6] + u_z_0[i])/2;
            u_z[j + 7] = (u_z_0[j + 7] + u_z_0[i])/2;
            }
        for (j = count_atoms - (count_atoms % 8); j < count_atoms; j++)
            {
            u_x[j] = (u_x_0[j] + u_x_0[i])/2;
            u_y[j] = (u_y_0[j] + u_y_0[i])/2;
            u_z[j] = (u_z_0[j] + u_z_0[i])/2;
            }
        for (j = i; j < count_atoms; j++)
            {
            axe_symetry = true;
            for (k = 0; k + 7 < count_atoms; k+=8) // Compute planes between lines and points
                { // compute d parameter of plane equation a = u_x[j], b = u_y[j], and c = u_z[j] for equation ax + by + cz + d = 0
                d[k] = -(u_x[j] * x[k] + u_y[j] * y[k] + u_z[j] * z[k]); 
                d[k + 1] = -(u_x[j] * x[k + 1] + u_y[j] * y[k + 1] + u_z[j] * z[k + 1]);
                d[k + 2] = -(u_x[j] * x[k + 2] + u_y[j] * y[k + 2] + u_z[j] * z[k + 2]);
                d[k + 3] = -(u_x[j] * x[k + 3] + u_y[j] * y[k + 3] + u_z[j] * z[k + 3]);
                d[k + 4] = -(u_x[j] * x[k + 4] + u_y[j] * y[k + 4] + u_z[j] * z[k + 4]);
                d[k + 5] = -(u_x[j] * x[k + 5] + u_y[j] * y[k + 5] + u_z[j] * z[k + 5]);
                d[k + 6] = -(u_x[j] * x[k + 6] + u_y[j] * y[k + 6] + u_z[j] * z[k + 6]);
                d[k + 7] = -(u_x[j] * x[k + 7] + u_y[j] * y[k + 7] + u_z[j] * z[k + 7]);
                // substitution of line equation to plane equation and computing t parameter
                t[k] = -((u_x[j] * x_center) + (u_y[j] * y_center) + (u_z[j] * z_center) + d[k])/
                (u_x[j] * u_x[j] + u_y[j] * u_y[j] + u_z[j] * u_z[j]);
                t[k + 1] = -((u_x[j] * x_center) + (u_y[j] * y_center) + (u_z[j] * z_center) + d[k + 1])/
                (u_x[j] * u_x[j] + u_y[j] * u_y[j] + u_z[j] * u_z[j]);
                t[k + 2] = -((u_x[j] * x_center) + (u_y[j] * y_center) + (u_z[j] * z_center) + d[k + 2])/
                (u_x[j] * u_x[j] + u_y[j] * u_y[j] + u_z[j] * u_z[j]);
                t[k + 3] = -((u_x[j] * x_center) + (u_y[j] * y_center) + (u_z[j] * z_center) + d[k + 3])/
                (u_x[j] * u_x[j] + u_y[j] * u_y[j] + u_z[j] * u_z[j]);
                t[k + 4] = -((u_x[j] * x_center) + (u_y[j] * y_center) + (u_z[j] * z_center) + d[k + 4])/
                (u_x[j] * u_x[j] + u_y[j] * u_y[j] + u_z[j] * u_z[j]);
                t[k + 5] = -((u_x[j] * x_center) + (u_y[j] * y_center) + (u_z[j] * z_center) + d[k + 5])/
                (u_x[j] * u_x[j] + u_y[j] * u_y[j] + u_z[j] * u_z[j]);
                t[k + 6] = -((u_x[j] * x_center) + (u_y[j] * y_center) + (u_z[j] * z_center) + d[k + 6])/
                (u_x[j] * u_x[j] + u_y[j] * u_y[j] + u_z[j] * u_z[j]);
                t[k + 7] = -((u_x[j] * x_center) + (u_y[j] * y_center) + (u_z[j] * z_center) + d[k + 7])/
                (u_x[j] * u_x[j] + u_y[j] * u_y[j] + u_z[j] * u_z[j]);
                // computing mirroring points
                x_mirroring[k] = - x[k] + 2 * (u_x[j] * t[k] + x_center);
                x_mirroring[k + 1] = - x[k + 1] + 2 * (u_x[j] * t[k + 1] + x_center);
                x_mirroring[k + 2] = - x[k + 2] + 2 * (u_x[j] * t[k + 2] + x_center);
                x_mirroring[k + 3] = - x[k + 3] + 2 * (u_x[j] * t[k + 3] + x_center);
                x_mirroring[k + 4] = - x[k + 4] + 2 * (u_x[j] * t[k + 4] + x_center);
                x_mirroring[k + 5] = - x[k + 5] + 2 * (u_x[j] * t[k + 5] + x_center);
                x_mirroring[k + 6] = - x[k + 6] + 2 * (u_x[j] * t[k + 6] + x_center);
                x_mirroring[k + 7] = - x[k + 7] + 2 * (u_x[j] * t[k + 7] + x_center);
                y_mirroring[k] = - y[k] + 2 * (u_y[j] * t[k] + y_center);
                y_mirroring[k + 1] = - y[k + 1] + 2 * (u_y[j] * t[k + 1] + y_center);
                y_mirroring[k + 2] = - y[k + 2] + 2 * (u_y[j] * t[k + 2] + y_center);
                y_mirroring[k + 3] = - y[k + 3] + 2 * (u_y[j] * t[k + 3] + y_center);
                y_mirroring[k + 4] = - y[k + 4] + 2 * (u_y[j] * t[k + 4] + y_center);
                y_mirroring[k + 5] = - y[k + 5] + 2 * (u_y[j] * t[k + 5] + y_center);
                y_mirroring[k + 6] = - y[k + 6] + 2 * (u_y[j] * t[k + 6] + y_center);
                y_mirroring[k + 7] = - y[k + 7] + 2 * (u_y[j] * t[k + 7] + y_center);
                z_mirroring[k] = - z[k] + 2 * (u_z[j] * t[k] + z_center);
                z_mirroring[k + 1] = - z[k + 1] + 2 * (u_z[j] * t[k + 1] + z_center);
                z_mirroring[k + 2] = - z[k + 2] + 2 * (u_z[j] * t[k + 2] + z_center);
                z_mirroring[k + 3] = - z[k + 3] + 2 * (u_z[j] * t[k + 3] + z_center);
                z_mirroring[k + 4] = - z[k + 4] + 2 * (u_z[j] * t[k + 4] + z_center);
                z_mirroring[k + 5] = - z[k + 5] + 2 * (u_z[j] * t[k + 5] + z_center);
                z_mirroring[k + 6] = - z[k + 6] + 2 * (u_z[j] * t[k + 6] + z_center);
                z_mirroring[k + 7] = - z[k + 7] + 2 * (u_z[j] * t[k + 7] + z_center);
            
                for (l = k; l <  k + 8; l++) // search mirroring points in other points
                    {
                    axe_symetry_flag = false;
                    if (j == i)
                        continue;
                
                    for (m = 0; m < count_atoms; m++)
                        {
                        if (x[m] < x_mirroring[l] + deviation  and x[m] > x_mirroring[l] - deviation
                        and y[m] < y_mirroring[l] + deviation  and y[m] > y_mirroring[l] - deviation
                        and z[m] < z_mirroring[l] + deviation  and z[m] > z_mirroring[l] - deviation
                        and Z[l] == Z[k])
                            {
                            axe_symetry_flag = true; // found
                            break;
                            }
                        }
                    if (axe_symetry_flag == false)
                        {
                        axe_symetry = false;
                        break;
                        }
                    }
                if (axe_symetry == false)
                   break;
                }
            for (k = (count_atoms - count_atoms % 8); k < count_atoms; k++) // Compute planes between lines and points
                {
                d[k] = -(u_x[j] * x[k] + u_y[j] * y[k] + u_z[j] * z[k]);
                t[k] = -((u_x[j] * x_center) + (u_x[j] * y_center) + (u_x[j] * z_center) + d[k])/
                (u_x[j] * u_x[j] + u_y[j] * u_y[j] + u_z[j] * u_z[j]);
                x_mirroring[k] = - x[k] + 2 * (u_x[j] * t[k] + x_center);
                y_mirroring[k] = - y[k] + 2 * (u_y[j] * t[k] + y_center);
                z_mirroring[k] = - z[k] + 2 * (u_z[j] * t[k] + z_center);
                
                for (l = k; l <  k + 1; l++)
                    {
                    axe_symetry_flag = false;
                    for (m = 0; m < count_atoms; m++)
                        {
                        if (x[m] < x_mirroring[l] + deviation  and x[m] > x_mirroring[l] - deviation
                        and y[m] < y_mirroring[l] + deviation  and y[m] > y_mirroring[l] - deviation
                        and z[m] < z_mirroring[l] + deviation  and z[m] > z_mirroring[l] - deviation
                        and Z[l] == Z[k])
                            {
                            axe_symetry_flag = true; // found
                            break;
                            }
                        }
                    if (axe_symetry_flag == false)
                        {
                        axe_symetry = false;
                        break;
                        }
                    }
                if (axe_symetry == false)
                   break;
                }
            if (axe_symetry == true)
                {
                u_x_2[count_axes] = u_x[j];
                u_y_2[count_axes] = u_y[j];
                u_z_2[count_axes] = u_z[j];
                count_axes++;
                }
            }
        }
    for (i = 0; i < count_atoms; i++) // search for symetry planes
        {
        x_center = x[i]; // For lines between points are center parameters first point
        y_center = y[i];
        z_center = z[i];
        
        for (j = 0; j + 7 < count_atoms; j+=8) // Fill lines paramaters between two atoms - differences between atoms
            {
            a[j] = x[j] - x_center;
            a[j + 1] = x[j + 1] - x_center;
            a[j + 2] = x[j + 2] - x_center;
            a[j + 3] = x[j + 3] - x_center;
            a[j + 4] = x[j + 4] - x_center;
            a[j + 5] = x[j + 5] - x_center;
            a[j + 6] = x[j + 6] - x_center;
            a[j + 7] = x[j + 7] - x_center;
            b[j] = y[j] - y_center;
            b[j + 1] = y[j + 1] - y_center;
            b[j + 2] = y[j + 2] - y_center;
            b[j + 3] = y[j + 3] - y_center;
            b[j + 4] = y[j + 4] - y_center;
            b[j + 5] = y[j + 5] - y_center;
            b[j + 6] = y[j + 6] - y_center;
            b[j + 7] = y[j + 7] - y_center;
            c[j] = z[j] - z_center;
            c[j + 1] = z[j + 1] - z_center;
            c[j + 2] = z[j + 2] - z_center;
            c[j + 3] = z[j + 3] - z_center;
            c[j + 4] = z[j + 4] - z_center;
            c[j + 5] = z[j + 5] - z_center;
            c[j + 6] = z[j + 6] - z_center;
            c[j + 7] = z[j + 7] - z_center;
            }
        for (j = count_atoms - (count_atoms % 8); j < count_atoms; j++) 
            { 
            a[j] = x[j] - x_center;
            b[j] = y[j] - y_center;
            c[j] = z[j] - z_center;
            }
        for (j = i + 1; j + 7 < count_atoms; j += 8) // Compute planes between doubles of points - substitution to lines equations
            { //  a, b and c from line equations parameters are substitued  by a, b and c 
            d[j] = -(a[j] * (x_center + 0.5 * a[j]) + b[j] * (y_center + 0.5 * b[j]) + c[j] * (z_center + 0.5 * c[j]));
            d[j + 1] = -(a[j + 1] * (x_center + 0.5 * a[j + 1]) + b[j + 1] * (y_center + 0.5 * b[j + 1])
            + c[j + 1] * (z_center + 0.5 * c[j + 1]));
            d[j + 2] = -(a[j + 2] * (x_center + 0.5 * a[j + 2]) + b[j + 2] * (y_center + 0.5 * b[j + 2])
            + c[j + 2] * (z_center + 0.5 * c[j + 2]));
            d[j + 3] = -(a[j + 3] + (x_center * 0.5 * a[j + 3]) + b[j + 3] * (y_center + 0.5 * b[j + 3])
            + c[j + 3] * (z_center + 0.5 * c[j + 3]));
            d[j + 4] = -(a[j + 4] * (x_center + 0.5 * a[j + 4]) + b[j + 4] * (y_center + 0.5 * b[j + 4])
            + c[j + 4] * (z_center + 0.5 * c[j + 4]));
            d[j + 5] = -(a[j + 5] * (x_center + 0.5 * a[j + 5]) + b[j + 5] * (y_center + 0.5 * b[j + 5])
            + c[j + 5] * (z_center + 0.5 * c[j + 5]));
            d[j + 6] = -(a[j + 6] * (x_center + 0.5 * a[j + 6]) + b[j + 6] * (y_center + 0.5 * b[j + 6])
            + c[j + 6] * (z_center + 0.5 * c[j + 6]));
            d[j + 7] = -(a[j + 7] * (x_center + 0.5 * a[j + 7]) + b[j + 7] * (y_center + 0.5 * b[j + 7])
            + c[j + 7] * (z_center + 0.5 * c[j + 7]));
            }
        for (j = count_atoms - (count_atoms % 8); j < count_atoms; j++)
            d[j] = -(a[j] * (x_center + 0.5 * a[j]) + b[j] * (y_center + 0.5 * b[j]) + c[j] * (z_center + 0.5 * c[j]));
        
        for (j = i + 1; j < count_atoms; j++)
            {
            plane_symetry = true;
            for (k = 0; k + 7 < count_atoms; k += 8) // Verify symetry of planes
                { // From points and planes construct curves and substitute to plane equations to compute mirroring points
                // a, b and c parameters are substitued from line equations by a, uy and c
                // computing t parameter
                t[k] = -(x[k] * a[j] + y[k] * b[j] + z[k] * c[j] + d[j])/(a[j] * a[j] + b[j] * b[j] + c[j] * c[j]);
                t[k + 1] = -(x[k + 1] * a[j] + y[k + 1] * b[j] + z[k + 1] * c[j] + d[j])/
                (a[j] * a[j] + b[j] * b[j] + c[j] * c[j]);
                t[k + 2] = -(x[k + 2] * a[j] + y[k + 2] * b[j] + z[k + 2] * c[j] + d[j])/
                (a[j] * a[j] + b[j] * b[j] + c[j] * c[j]);
                t[k + 3] = -(x[k + 3] * a[j] + y[k + 3] * b[j] + z[k + 3] * c[j] + d[j])/
                (a[j] * a[j] + b[j] * b[j] + c[j] * c[j]);
                t[k + 4] = -(x[k + 4] * a[j] + y[k + 4] * b[j] + z[k + 4] * c[j] + d[j])/
                (a[j] * a[j] + b[j] * b[j] + c[j] * c[j]);
                t[k + 5] = -(x[k + 5] * a[j] + y[k + 5] * b[j] + z[k + 5] * c[j] + d[j])/
                (a[j] * a[j] + b[j] * b[j] + c[j] * c[j]);
                t[k + 6] = -(x[k + 6] * a[j] + y[k + 6] * b[j] + z[k + 6] * c[j] + d[j])/
                (a[j] * a[j] + b[j] * b[j] + c[j] * c[j]);
                t[k + 7] = -(x[k + 7] * a[j] + y[k + 7] * b[j] + z[k + 7] * c[j] + d[j])/
                (a[j] * a[j] + b[j] * b[j] + c[j] * c[j]);
                // computing mirroring points
                x_mirroring[k] = x[k] + a[j] * 2 * t[k];
                x_mirroring[k + 1] = x[k + 1] + a[j] * 2 * t[k + 1];
                x_mirroring[k + 2] = x[k + 2] + a[j] * 2 * t[k + 2];
                x_mirroring[k + 3] = x[k + 3] + a[j] * 2 * t[k + 3];
                x_mirroring[k + 4] = x[k + 4] + a[j] * 2 * t[k + 4];
                x_mirroring[k + 5] = x[k + 5] + a[j] * 2 * t[k + 5];
                x_mirroring[k + 6] = x[k + 6] + a[j] * 2 * t[k + 6];
                x_mirroring[k + 7] = x[k + 7] + a[j] * 2 * t[k + 7];
                y_mirroring[k] = y[k] + b[j] * 2 * t[k];
                y_mirroring[k + 1] = y[k + 1] + b[j] * 2 * t[k + 1];
                y_mirroring[k + 2] = y[k + 2] + b[j] * 2 * t[k + 2];
                y_mirroring[k + 3] = y[k + 3] + b[j] * 2 * t[k + 3];
                y_mirroring[k + 4] = y[k + 4] + b[j] * 2 * t[k + 4];
                y_mirroring[k + 5] = y[k + 5] + b[j] * 2 * t[k + 5];
                y_mirroring[k + 6] = y[k + 6] + b[j] * 2 * t[k + 6];
                y_mirroring[k + 7] = y[k + 7] + b[j] * 2 * t[k + 7];
                z_mirroring[k] = z[k] + c[j] * 2 * t[k];
                z_mirroring[k + 1] = z[k + 1] + c[j] * 2 * t[k + 1];
                z_mirroring[k + 2] = z[k + 2] + c[j] * 2 * t[k + 2];
                z_mirroring[k + 3] = z[k + 3] + c[j] * 2 * t[k + 3];
                z_mirroring[k + 4] = z[k + 4] + c[j] * 2 * t[k + 4];
                z_mirroring[k + 5] = z[k + 5] + c[j] * 2 * t[k + 5];
                z_mirroring[k + 6] = z[k + 6] + c[j] * 2 * t[k + 6];
                z_mirroring[k + 7] = z[k + 7] + c[j] * 2 * t[k + 7];
                // control routine
                for (l = k; l <  k + 8; l++) // search mirroring points in other points
                    {
                    if (l == i or l == j)
                        continue;
                    
                    plane_symetry_flag = false;
                    for (m = 0; m < count_atoms; m++)
                        {
                        if (x[m] < x_mirroring[l] + deviation  and x[m] > x_mirroring[l] - deviation
                        and y[m] < y_mirroring[l] + deviation  and y[m] > y_mirroring[l] - deviation
                        and z[m] < z_mirroring[l] + deviation  and z[m] > z_mirroring[l] - deviation
                        and Z[m] == Z[l])
                            {
                            plane_symetry_flag = true; // found
                            break;
                            }
                        }
                    if (plane_symetry_flag == false)
                        {
                        plane_symetry = false;
                        break; // break actual control
                        }
                    }
                if (plane_symetry == false)
                    break; // break cyclus of control
                }
            if (plane_symetry == false)
                break; // not continue to control
                
            for (k = count_atoms - (count_atoms % 8); k < count_atoms; k ++) // Verify symetry of planes
                { // From points and planes construct curves and substitute to plane equations to compute mirroring points
                // a, b and c parameters are substitued from line equations by a, uy and c
                // compute t parameter
                t[k] = -(x[k] * a[j] + y[k] * b[j] + z[k] * c[j] + d[j])/(a[j] * a[j] + b[j] * b[j] + c[j] * c[j]);
                // compute mirroring points
                x_mirroring[k] = x[k] + a[j] * 2 * t[k];
                y_mirroring[k] = y[k] + b[j] * 2 * t[k];
                z_mirroring[k] = z[k] + c[j] * 2 * t[k];
                // control routine
                plane_symetry = true;
                for (l = k; l <  k + 1; l++) // searh mirroring points in other points
                    {
                    if (l == i or l == j)
                        continue;
                    
                    plane_symetry_flag = false;
                    
                    for (m = 0; m < count_atoms; m++)
                        {
                        if (x[m] < x_mirroring[l] + deviation  and x[m] > x_mirroring[l] - deviation
                        and y[m] < y_mirroring[l] + deviation  and y[m] > y_mirroring[l] - deviation
                        and z[m] < z_mirroring[l] + deviation  and z[m] > z_mirroring[l] - deviation
                        and Z[m] == Z[l])
                            {
                            plane_symetry_flag = true; // found
                            break;
                            }
                        }
                    if (plane_symetry_flag == false)
                        {
                        plane_symetry = false;
                        break;
                        }
                    }
                if (plane_symetry == false)
                   break;
                }
            if (plane_symetry == true)
                {
                a_2[count_planes] = a[j];
                b_2[count_planes] = b[j];
                c_2[count_planes] = c[j];
                d_2[count_planes] = d[j];
                count_planes++;
                }
            } 
        }
    symetry_axes->x_center = x_center;
    symetry_axes->y_center = y_center;
    symetry_axes->z_center = z_center;
    symetry_axes->symetry_center = center_symetry;
    for (i = 0; i < count_axes; i++) // copying symetric axes to vectors
        {
        previous_finded = false;
        for (j = 0; j < i; j++)
            {
            if (u_x_2[j] != 0)
                ratio = u_x_2[i]/u_x_2[j];
            else
                {
                if (u_y_2[j] != 0)
                    ratio = u_y_2[i]/u_y_2[j];
                else
                    ratio = u_z_2[i]/u_z_2[j];
                }
            if ((u_x_2[i] > (u_x_2[j] - deviation) * ratio and u_x_2[i] < (u_x_2[j] + deviation) * ratio
            and u_y_2[i] > (u_y_2[j] - deviation) * ratio and u_y_2[i] < (u_y_2[j] + deviation) * ratio
            and u_z_2[i] > (u_z_2[j] - deviation) * ratio and u_z_2[i] < (u_z_2[j] + deviation) * ratio)
            )
                {
                previous_finded = true;
                break;
                }
            }
        if (previous_finded == true)
            continue;
            
        symetry_axes->u_x.push_back(u_x_2[i]);
        symetry_axes->u_y.push_back(u_y_2[i]);
        symetry_axes->u_z.push_back(u_z_2[i]);
        }
    for (i = 0; i < count_planes; i++) // copying symetric axes to vectors
        {
        previous_finded = false;
        for (j = 0; j < i; j++)
            {
            if (a_2[j] != 0)
                    ratio = a_2[i]/a_2[j];
                else
                    {
                    if (b_2[j] != 0)
                        ratio = b_2[i]/b_2[j];
                    else
                        ratio = c_2[i]/c_2[j];
                    }
            if (a_2[i] > a_2[j] * ratio - deviation and a_2[i] < a_2[j] * ratio + deviation
            and b_2[i] > b_2[j] * ratio - deviation and b_2[i] < b_2[j] * ratio + deviation
            and c_2[i] > c_2[j] * ratio - deviation and c_2[i] < c_2[j] * ratio + deviation)
                {
                previous_finded = true;
                break;
                }
            }
        if (previous_finded == true)
            continue;
        
        symetry_planes->a.push_back(a_2[i]);
        symetry_planes->b.push_back(b_2[i]);
        symetry_planes->c.push_back(c_2[i]);
        symetry_planes->d.push_back(d_2[i]);
        }
    // Here begins detecting symetry information from location of atoms
    // First two atoms
    u_A[0] = x[0];
    u_A[1] = y[0];
    u_A[2] = z[0];
    u_B[0] = x[1];
    u_B[1] = y[1];
    u_B[2] = z[1];
    u_A_B[0] = u_A[0] - u_B[0];
    u_A_B[1] = u_A[1] - u_B[1];
    u_A_B[2] = u_A[2] - u_B[2];
    
    if (count_atoms < 2)
        symetry_axes_parameters.linear = false;
    else
        symetry_axes_parameters.linear = true;
        
    for (i = 2; i < count_atoms; i++)
        {
        u_C[0] = x[i];
        u_C[1] = y[i];
        u_C[2] = z[i];
        u_A_C[0] = u_A[0] - u_C[0];
        u_A_C[1] = u_A[1] - u_C[1];
        u_A_C[2] = u_A[2] - u_C[2];
        if (u_A_C[0]/u_A_B[0] < (u_A_C[1]/u_A_B[1] - deviation) or u_A_C[0]/u_A_B[0] > (u_A_C[1]/u_A_B[1] + deviation)
        or  u_A_C[0]/u_A_B[0] < (u_A_C[2]/u_A_B[2] - deviation) or u_A_C[0]/u_A_B[0] > (u_A_C[2]/u_A_B[2] + deviation))
            symetry_axes_parameters.linear = false;
        }
    // Construction curve equation from normal vector of three atoms and verifying planarity
    if (symetry_axes_parameters.linear == false and count_atoms >= 3)
        {
        u_C[0] = x[2];
        u_C[1] = y[2];
        u_C[2] = z[2];
        u_A_C[0] = u_A[0] - u_C[0];
        u_A_C[1] = u_A[1] - u_C[1];
        u_A_C[2] = u_A[2] - u_C[2];
        D3_vector_multiply(&u_A_B[0], &u_A_C[0], &planar_equation[0]);
        planar_equation[3] = planar_equation[0] * u_A[0] + planar_equation[1] * u_A[1] + planar_equation[2] * u_A[2];
        symetry_planes_parameters.planar = false;
        i = 2;
        
        while (i < count_atoms)
            {
            if (abs(x[i] * planar_equation[0] + y[i] * planar_equation[1] + z[i] * planar_equation[2] + planar_equation[3]) > deviation)
                {
                symetry_planes_parameters.planar = false;
                break;
                }
            i++;
            }
        if (i == count_atoms)
            symetry_planes_parameters.planar = true;
        }
    if (symetry_planes_parameters.planar == true)
        {
        symetry_planes_parameters.a.push_back(planar_equation[0]);
        symetry_planes_parameters.b.push_back(planar_equation[0]);
        symetry_planes_parameters.c.push_back(planar_equation[0]);
        symetry_planes_parameters.d.push_back(planar_equation[0]);
        }
    return(0);
    }
template <typename T>
T Wavefunction_calculations<T>::D3_vector_multiply(T* a, T* b, T* c)
    {c[0] = (a[1] * b[2]) - (a[2] * b[1]); c[1] = (a[2] * b[0]) - (a[0] * b[2]); c[2] = (a[0] * b[1]) - (a[1] * b[0]);
    return(0);}
template <typename T>
T Wavefunction_calculations<T>::Clear()
    {
    count_atoms = 0;
    x_center_of_mass = 0;
    y_center_of_mass = 0;
    z_center_of_mass = 0;
    electron_to_atom_numbers.clear();
    index_atoms.clear();
    energies = nullptr;
    dipoles.partial_charges.clear();
    dipoles.x.clear();
    dipoles.y.clear();
    dipoles.z.clear();
    dipoles.x_dipole_moments.clear();
    dipoles.y_dipole_moments.clear();
    dipoles.z_dipole_moments.clear();
    dipoles.sum_dipole_moment = 0;
    dipoles.sum_dipole_moment_Debye = 0;
    central_cation_atoms.atom_numbers.clear();
    central_cation_atoms.configuration_numbers.clear();
    central_cation_atoms.high_spin.clear();
    symetry_axes_parameters.u_x.clear();
    symetry_axes_parameters.u_y.clear();
    symetry_axes_parameters.u_z.clear();
    symetry_axes_parameters.x_center = 0;
    symetry_axes_parameters.y_center = 0;
    symetry_axes_parameters.z_center = 0;
    symetry_axes_parameters.symetry_center = 0;
    symetry_planes_parameters.a.clear();
    symetry_planes_parameters.b.clear();
    symetry_planes_parameters.c.clear();
    symetry_planes_parameters.d.clear();
    return(0);
    }
template <typename T>
Wavefunction_calculations<T>::~Wavefunction_calculations(){
    Clear();}
/*
  Author of this source code Ing. Pavel Florian Ph.D. licensed this source code under the the 3-Clause BSD License:
  
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimerin the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/ 
