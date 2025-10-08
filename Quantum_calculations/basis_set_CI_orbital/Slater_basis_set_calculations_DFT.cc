#include "Slater_basis_set_calculations_DFT.h"
using namespace std;

template <typename T>
int Slater_basis_set_calculations_DFT<T>::Create_index_atoms()
    {
    unsigned int i;
    unsigned int sum_electrons = this->results.n.size();
    unsigned int count_bonds;
    T* x = this->results.x.data();
    T* y = this->results.y.data();
    T* z = this->results.z.data();
    
    index_atoms.clear();
    electron_to_atom_numbers.clear();

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
int Slater_basis_set_calculations_DFT<T>::Detect_center_of_mass()
    {
    unsigned int i;
    unsigned int Z = 0;
    unsigned int sum_Z = 0;
    T sum_x = 0;
    T sum_y = 0;
    T sum_z = 0;
    
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
int Slater_basis_set_calculations_DFT<T>::Detect_dipole_moments(dipole_moment *dipole_moment)
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
    
    dipole_moment->partial_charges.clear();
    dipole_moment->x.clear();
    dipole_moment->y.clear();
    dipole_moment->z.clear();
    dipole_moment->x_dipole_moments.clear();
    dipole_moment->y_dipole_moments.clear();
    dipole_moment->z_dipole_moments.clear();
    
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
int Slater_basis_set_calculations_DFT<T>::Create_crystal_field(central_cations *central_cations)
    {
    unsigned int i, j, k;
    unsigned int atom_begin, next_atom_begin;
    unsigned int size;
    unsigned int* n;
    unsigned int* l;
    int* m;
    array<T, this->max_electrons> atom_numbers;
    array<T, this->max_electrons> configuration_numbers;
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
int Slater_basis_set_calculations_DFT<T>::Allocate_densities(unsigned int size_order)
    {
    unsigned int i;
    unsigned int density_size;
    T* density = nullptr;
    
    density_size = (2 * size_order + 1) * (2 * size_order + 1) * (2 * size_order + 1);
    try {
        for (i = 0; i < index_atoms.size(); i++)
            {
            density = new T[density_size];
            atoms_spin_densities.push_back(density);
            }
        }
    catch (int){
        for (i = 0; i < atoms_spin_densities.size(); i++)
            delete [] atoms_spin_densities[i];
        
        return(-1);
        }
    try {
        for (i = 0; i < index_atoms.size(); i++)
            {
            density = new T[density_size * 2];
            atoms_electron_densities.push_back(density);
            }
        }
    catch (int){
        for (i = 0; i < atoms_electron_densities.size(); i++)
            delete [] atoms_electron_densities[i];
        
        return(-1);
        }
    try {
        for (i = 0; i < index_atoms.size(); i++)
            {
            density = new T[density_size * 2];
            atoms_Fi_densities.push_back(density);
            }
        }
    catch (int) {
        for (i = 0; i < atoms_Fi_densities.size(); i++)
            delete [] atoms_Fi_densities[i];
        
        return(-1);
        }
    try {
        for (i = 0; i < index_atoms.size(); i++)
            {
            density = new T[density_size * 2];
            atoms_gradients_densities.push_back(density);
            }
        }
    catch (int) {
        for (i = 0; i < atoms_gradients_densities.size(); i++)
            delete [] atoms_gradients_densities[i];
        return(-1);
        }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Compute_density_thread(T** densities, T* atom_density, unsigned int begin, unsigned int end, unsigned int size_order, vector<T> spins, vector<int> spin_paired,
vector<unsigned int> x_range, vector<unsigned int> y_range, vector<unsigned int> z_range, bool spin_density)
    {
    unsigned int i, j, k, l;
    unsigned int j_min, k_min, l_min;
    unsigned int j_max, k_max, l_max;
    unsigned int side;
    unsigned int size;
    unsigned int pre_index;
    unsigned int multiplicity = 1;
    T* pointer_to_density;
    
    side = 2 * size_order + 1;
    size = side * side * side;
    // Vectorisation code
    // memset to zeros;
    for (i = 0; i + 7 < size; i++)
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
    for (i = size - (size % 8); i < size; i++)
        atom_density[i] = 0;
    
    for (i = begin; i <= end; i++)
        {
        j_min = size_order - z_range[i];
        k_min = size_order - y_range[i];
        l_min = size_order - x_range[i];
        j_max = size_order + z_range[i] + 1;
        k_max = size_order + y_range[i] + 1;
        l_max = size_order + x_range[i] + 1;
        
        if (spin_density == false or spin_paired[i] == -1 or spin_paired[i] < begin
        or spin_paired[i] > end)
            {
            pointer_to_density = densities[i]; // Unpaired and bonding electrons are included for spin_paired[i] < begin
            // spin-paired orbitals in closed-shell method are added 1x
            if (multiplicity == 2)
                {
                multiplicity = 1;
                continue;
                }
            if (i < end and spin_density == false)
                if (densities[i] == densities[i + 1])
                    {
                    multiplicity = 2;
                    if (spin_density == true)
                        continue;
                    }
            if (spin_density == true and spins[i] == -0.5) // For spin density map and negative spins
                {
                for (j = j_min; j < j_max; j++)
                    for (k = k_min; k < k_max; k++)
                        {
                        pre_index = j * side * side + k * side;
                        for (l = l_min; l + 7 < l_max; l+=8)
                            {
                            atom_density[pre_index + l] -= pointer_to_density[pre_index + l];
                            atom_density[pre_index + l + 1] -= pointer_to_density[pre_index + l + 1];
                            atom_density[pre_index + l + 2] -= pointer_to_density[pre_index + l + 2];
                            atom_density[pre_index + l + 3] -= pointer_to_density[pre_index + l + 3];
                            atom_density[pre_index + l + 4] -= pointer_to_density[pre_index + l + 4];
                            atom_density[pre_index + l + 5] -= pointer_to_density[pre_index + l + 5];
                            atom_density[pre_index + l + 6] -= pointer_to_density[pre_index + l + 6];
                            atom_density[pre_index + l + 7] -= pointer_to_density[pre_index + l + 7];
                            }
                        for (l = l_max - (l_max % 8); l < l_max; l++)
                            {
                            atom_density[pre_index + l] -= pointer_to_density[pre_index + l];
                            }
                        }
                }
            else
                {
                if (multiplicity == 1)
                    {
                    for (j = j_min; j < j_max; j++)
                        for (k = k_min; k < k_max; k++)
                            {
                            pre_index = j * side * side + k * side;
                            for (l = l_min; l + 7 < l_max; l+=8)
                                {
                                atom_density[pre_index + l] += pointer_to_density[pre_index + l];
                                atom_density[pre_index + l + 1] += pointer_to_density[pre_index + l + 1];
                                atom_density[pre_index + l + 2] += pointer_to_density[pre_index + l + 2];
                                atom_density[pre_index + l + 3] += pointer_to_density[pre_index + l + 3];
                                atom_density[pre_index + l + 4] += pointer_to_density[pre_index + l + 4];
                                atom_density[pre_index + l + 5] += pointer_to_density[pre_index + l + 5];
                                atom_density[pre_index + l + 6] += pointer_to_density[pre_index + l + 6];
                                atom_density[pre_index + l + 7] += pointer_to_density[pre_index + l + 7];
                                }
                            for (l = l_max - (l_max % 8); l < l_max; l++)
                                {
                                atom_density[pre_index + l] += pointer_to_density[pre_index + l];
                                }
                        }
                    }
                else
                    {
                    for (j = j_min; j < j_max; j++)
                        for (k = k_min; k < k_max; k++)
                            {
                            pre_index = j * side * side + k * side;
                            for (l = l_min; l + 7 < l_max; l+=8)
                                {
                                atom_density[pre_index + l] += pointer_to_density[pre_index + l] * multiplicity;
                                atom_density[pre_index + l + 1] += pointer_to_density[pre_index + l + 1] * multiplicity;
                                atom_density[pre_index + l + 2] += pointer_to_density[pre_index + l + 2] * multiplicity;
                                atom_density[pre_index + l + 3] += pointer_to_density[pre_index + l + 3] * multiplicity;
                                atom_density[pre_index + l + 4] += pointer_to_density[pre_index + l + 4] * multiplicity;
                                atom_density[pre_index + l + 5] += pointer_to_density[pre_index + l + 5] * multiplicity;
                                atom_density[pre_index + l + 6] += pointer_to_density[pre_index + l + 6] * multiplicity;
                                atom_density[pre_index + l + 7] += pointer_to_density[pre_index + l + 7] * multiplicity;
                                }
                            for (l = l_max - (l_max % 8); l < l_max; l++)
                                {
                                atom_density[pre_index + l] += pointer_to_density[pre_index + l] * multiplicity;
                                }
                            }
                    }
                }
            }
        }
    // End of vectorisation code
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Compute_densities(vector<T*> densities, vector<T*> atoms_densities,
vector<unsigned int> index, unsigned int size_order, vector<T> spins, vector<int> spin_paired,
vector<unsigned int> x_range, vector<unsigned int> y_range, vector<unsigned int> z_range, bool spin_density)
    { // Compute densities of wavefunctions/probabilities for atoms in 1 thread for each atoms.
    unsigned int i, j;
    unsigned int begin, end;
    unsigned int count_atoms;
    unsigned int size_atoms;
    array<unsigned int, this->max_atoms> begins_atoms;
    array<unsigned int, this->max_atoms> ends_atoms;
    array<T*, this->max_atoms> pointers_to_atoms;
    array<T*, this->max_atoms> pointers_to_probabilities;

    count_atoms = index.size();
    if (count_atoms > this->max_atoms) // Control, that avoid throwing an exception.
        return(-1);
    
    size_atoms = (size_order * 2 + 1) * (size_order * 2 + 1) * (size_order * 2 + 1);
    i = 0;
    while (i < count_atoms) // Filling begins_atoms and ends_atoms array.
        {
        begins_atoms[i] = index[i];
        if (i + 1 < count_atoms)
            ends_atoms[i] = index[i + 1];
        else
            ends_atoms[i] = densities.size() -1;
        i++;
       }
   for (i = 0; i < count_atoms; i++)
        pointers_to_atoms[i] = atoms_densities[i];
        
   for (i = 0; i < densities.size(); i++)
       pointers_to_probabilities[i] = densities[i];
   // Multithreading code
   #pragma omp parallel
        {
        #pragma omp for
        for (i = 0; i < count_atoms; i++)
            {
            Compute_density_thread(&pointers_to_probabilities[0], pointers_to_atoms[i],
            begins_atoms[i], ends_atoms[i], size_order,  spins, spin_paired, x_range,
            y_range, z_range, spin_density);
            }
        }
    // End of multithreading code  
    return(0);  
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Correct_densities(vector<T*> atoms_densities_list, vector<T>  x, vector <T> y, vector <T> z,
unsigned int size_order)
    {
    unsigned int i, j, k, l, m;
    unsigned int side;
    unsigned int x_side, y_side, z_side;
    unsigned int x_contraction, y_contraction, z_contraction;
    unsigned int x_1_min, y_1_min, z_1_min;
    unsigned int x_2_min, y_2_min, z_2_min;
    unsigned int x_condition, x_condition_2;
    int X, Y, Z;
    T* densities_1;
    T* densities_2;
    
    side = (2 * size_order + 1);
    for (i = 0; i < atoms_densities_list.size(); i++)
        for (j = 0; j < atoms_densities_list.size(); j++)
            {
            densities_1 = atoms_densities_list[i];
            densities_2 = atoms_densities_list[j];
            
            X = (x[index_atoms[j]] - x[index_atoms[i]]) * size_order/this->vector_lenght;
            Y = (y[index_atoms[j]] - y[index_atoms[i]]) * size_order/this->vector_lenght;
            Z = (z[index_atoms[j]] - z[index_atoms[i]]) * size_order/this->vector_lenght;
            x_contraction = abs(X);
            y_contraction = abs(Y);
            z_contraction = abs(Z);
            if (x_contraction > side)
                x_contraction = side;
            if (y_contraction > side)
                y_contraction = side;
            if (z_contraction > side)
                z_contraction = side;
    
            x_side = side - x_contraction;
            y_side = side - y_contraction;
            z_side = side - z_contraction;
            if (x_contraction >= side or y_contraction >= side or z_contraction >= side)
                continue;
        
            if (X < 0)
                {
                x_1_min = 0;
                x_2_min = x_contraction;
                }
            else
                {
                x_1_min = x_contraction;
                x_2_min = 0;
                }
            if (Y < 0)
                {
                y_1_min = 0;
                y_2_min = y_contraction;
                }
            else
                {
                y_1_min = y_contraction;
                y_2_min = 0;
                }
            if (Z < 0)
                {
                z_1_min = 0;
                z_2_min = z_contraction;
               }
            else
                {
                z_1_min = z_contraction;
                z_2_min = 0;
                }
            x_condition = side - x_contraction;   
            x_condition_2 = x_condition - (x_condition % 8);
            // vectorisation code
            for (k = 0; k < (side - z_contraction); k++)
                {
                for (l = 0; l < (side - y_contraction); l++)
                    {
                    for (m = 0; m + 7 < x_condition; m+=8)
                        {
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min)] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min)] +
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min)];
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min)] = 
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min)];
                        
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 1] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 1] +
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 1];
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 1] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 1];
                        
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 2] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 2] +
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 2];
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 2] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 2];
                        
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 3] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 3] +
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 3];
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 3] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 3];
                        
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 4] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 4] +
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 4];
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 4] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 4];
                        
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 5] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 5] +
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 5];
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 5] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 5];
                        
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 6] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 6] +
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 6];
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 6] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 6];
                        
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 7] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 7] +
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 7];
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min) + 7] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min) + 7];
                        }
                    for (m = x_condition_2; m < x_condition; m++)
                        {
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min)] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min)] +
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min)];
                        densities_2[((k + z_2_min) * side * side) + ((l + y_2_min) * side) + (m + x_2_min)] =
                        densities_1[((k + z_1_min) * side * side) + ((l + y_1_min) * side) + (m + x_1_min)];
                        }
                    }
                }
            }
    // end of vectorisation code
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Compute_Fi_thread(T* spin_density, T* electron_density, T* Fi_density, unsigned int size_order)
    {
    unsigned int i;
    unsigned int size;
    T pow_1, pow_2, pow_3, pow_4, pow_5, pow_6, pow_7, pow_8;
    
    size = (2 * size_order + 1) * (2 * size_order + 1) * (2 * size_order + 1);
    // Vectorisation code
    for (i = 0; i + 7 < size; i+= 8)
        {
        pow_1 = spin_density[i]/electron_density[i];
        pow_2 = spin_density[i + 1]/electron_density[i + 1];
        pow_3 = spin_density[i + 2]/electron_density[i + 2];
        pow_4 = spin_density[i + 3]/electron_density[i + 3];
        pow_5 = spin_density[i + 4]/electron_density[i + 4];
        pow_6 = spin_density[i + 5]/electron_density[i + 5];
        pow_7 = spin_density[i + 6]/electron_density[i + 6];
        pow_8 = spin_density[i + 7]/electron_density[i + 7];
        Fi_density[i] = (pow(1 + pow_1, 2.0/3.0) + pow(1 - pow_1, 2.0/3.0))/2;
        Fi_density[i + 1] = (pow(1 + pow_2, 2.0/3.0) + pow(1 - pow_2, 2.0/3.0))/2;
        Fi_density[i + 2] = (pow(1 + pow_3, 2.0/3.0) + pow(1 - pow_3, 2.0/3.0))/2;
        Fi_density[i + 3] = (pow(1 + pow_4, 2.0/3.0) + pow(1 - pow_4, 2.0/3.0))/2;
        Fi_density[i + 4] = (pow(1 + pow_5, 2.0/3.0) + pow(1 - pow_5, 2.0/3.0))/2;
        Fi_density[i + 5] = (pow(1 + pow_6, 2.0/3.0) + pow(1 - pow_6, 2.0/3.0))/2;
        Fi_density[i + 6] = (pow(1 + pow_7, 2.0/3.0) + pow(1 - pow_7, 2.0/3.0))/2;
        Fi_density[i + 7] = (pow(1 + pow_8, 2.0/3.0) + pow(1 - pow_8, 2.0/3.0))/2;
        }
    for (i = size - (size % 8); i < size; i++)
        {
        pow_1 = pow(spin_density[i], 2.0/3.0);
        Fi_density[i] = (pow(1 + pow_1, 2.0/3.0) + pow(1 - pow_1, 2.0/3.0))/2;
        }
    // End of Vectorisation code
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Compute_gradient_thread(T* electron_density, T* gradient_density, unsigned int size_order)
    {
    unsigned int i, j;
    unsigned int side;
    unsigned int size;
    T x_gradient, y_gradient, z_gradient;
    T pixel_lenght;
    
    pixel_lenght = T(this->vector_lenght)/T(size_order);
    side = (2 * size_order + 1);
    size = side * side * side;
    
    for (i = 0; i < side * side; i++) // Optimalization for pipelining with non sequential reading of memory
        {
        if (i >= 1)
            x_gradient = (electron_density[i + 1] - electron_density[i - 1])/(pixel_lenght * 2);
        else
            x_gradient = (electron_density[i + 1] - electron_density[i])/pixel_lenght;
        
        if (i >= side)
            y_gradient = (electron_density[i + side] - electron_density[i - side])/(pixel_lenght * 2);
        else
            y_gradient = (electron_density[i + side] - electron_density[i])/pixel_lenght;
        
        z_gradient = (electron_density[i + (side * side)] - electron_density[i])/pixel_lenght;
        gradient_density[i] = sqrt(((x_gradient * x_gradient) + (y_gradient * y_gradient) + (z_gradient * z_gradient))/3);
        }
    for (i = (side * side); i < size - (side * side); i++)
        {
        gradient_density[i - 1] = sqrt(((x_gradient * x_gradient) + (y_gradient * y_gradient) + (z_gradient * z_gradient))/3);
        x_gradient = (electron_density[i + 1] - electron_density[i - 1])/(pixel_lenght * 2);
        y_gradient = (electron_density[i + side] - electron_density[i - side])/(pixel_lenght * 2);
        z_gradient = (electron_density[i + (side * side)] - electron_density[i - (side * side)])/(pixel_lenght * 2);
        }
    gradient_density[i] = sqrt(((x_gradient * x_gradient) + (y_gradient * y_gradient) + (z_gradient * z_gradient)))/2.00;
    for (i = size - (side * side); i < size ; i++)
        {
        if (size - i > 1)
            x_gradient = (electron_density[i + 1] - electron_density[i - 1])/(pixel_lenght * 2);
        else
            x_gradient = (electron_density[i] - electron_density[i - 1])/pixel_lenght;
        
        if (size - i > side)
            y_gradient = (electron_density[i + side] - electron_density[i - side])/(pixel_lenght * 2);
        else
            y_gradient = (electron_density[i] - electron_density[i - side])/pixel_lenght;
        
        z_gradient = (electron_density[i] - electron_density[i - (side * side)])/(pixel_lenght);
        gradient_density[i] = sqrt(((x_gradient * x_gradient) + (y_gradient * y_gradient) + (z_gradient * z_gradient)));
        }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Compute_gradients(vector<T*> atoms_electron_densities, 
vector<T*> atoms_gradients_densities, unsigned int size_order)
    {
    unsigned int i;
    array<T*, this->max_atoms> electron_densities;
    array<T*, this->max_atoms> gradients_densities;
    
    for (i = 0; i < count_atoms; i++)
        {
        electron_densities[i] = atoms_electron_densities[i];
        gradients_densities[i] = atoms_gradients_densities[i];
        }
    // multithreading code
    #pragma omp parallel
        {
        #pragma omp for
        for (i = 0; i < count_atoms; i++)
            {
            Compute_gradient_thread(electron_densities[i], gradients_densities[i], size_order);
            }
        }
    // End of multithreading code
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Compute_Fi_and_gradients(vector<T*> atoms_electron_densities, vector<T*> atoms_spin_densities,
vector<T*> atoms_Fi_densities, vector<T*> atoms_gradients_densities, unsigned int size_order)
    {
    unsigned int i;
    
    array<T*, this->max_atoms> electron_densities;
    array<T*, this->max_atoms> spin_densities;
    array<T*, this->max_atoms> Fi_densities;
    
    for (i = 0; i < count_atoms; i++)
        {
        electron_densities[i] = atoms_electron_densities[i];
        spin_densities[i] = atoms_spin_densities[i];
        Fi_densities[i] = atoms_Fi_densities[i];
        }
    // multithreading code
    #pragma omp parallel
        {
        #pragma omp for
        for (i = 0; i < count_atoms; i++)
            {
            Compute_Fi_thread(spin_densities[i], electron_densities[i], Fi_densities[i], size_order);
            }
        }
    // End of multithreading code
    Compute_gradients(atoms_electron_densities, atoms_gradients_densities, size_order);
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::PBE_thread(T* atoms_gradients_density, T* atoms_Fi_density, T* electron_density,
unsigned int size_order, T* exchange_energy, T* correlation_energy, T Wigner_Seitz_radius)
    {
    unsigned int i, j;
    unsigned int size;
    T exchange_integral;
    T correlation_integral;
    T n;
    T k_F;
    T s;
    T pre_s;
    T c; // this is not sppeed of light in this context
    T X;
    T e_x_unif;
    T F_x;
    T nl_integral_factor;
    T s_2;
    T fi_2;
    T power;
    
    size = (2 * size_order + 1) * (2 * size_order + 1) * (2 * size_order + 1);
    // Correlation energy computation preparation
    n = (3 * this->Hartree_lenght * this->Hartree_lenght * this->Hartree_lenght)
    /(4 * this->Pi * Wigner_Seitz_radius * Wigner_Seitz_radius * Wigner_Seitz_radius);
    k_F = pow(n * 3 * this->Pi * this->Pi, 1.00/3.00);
    c = pow((3 * this->Pi * this->Pi)/16, 1.00/3.00);
    pre_s = 1.00/(2 * k_F);
    X = (beta/gamma) * c * c * exp(-omega/gamma); // equal 0.72161
    // Exgchange energy computation preparation
    e_x_unif = -(3 * this->e * k_F)/(4 * this->Pi);
    
    exchange_integral = 0;
    correlation_integral = 0;
    for (i = 0; i < size; i++) // Computing corelation_energy
        {
        if (electron_density[i] != 0)
            {
            s = pre_s * abs(atoms_gradients_density[i])/n;
            s_2 = s * s;
            fi_2 = atoms_Fi_density[i] * atoms_Fi_density[i];
            power = pow(X, T(s_2/fi_2));
            nl_integral_factor = log(1 + 1/(power + (power * power)));
            if ((not (isnan(nl_integral_factor))) and (not (isinf(nl_integral_factor))))
                correlation_integral = correlation_integral +
                (gamma * electron_density[i] * fi_2 * atoms_Fi_density[i] * nl_integral_factor);
        
            exchange_integral = exchange_integral + e_x_unif * (1 + kappa - kappa/(1 + (mi * s_2)/kappa)) * electron_density[i];
            }
        }
    exchange_energy[0] = exchange_integral;
    correlation_energy[0] = correlation_integral;
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::PBE_compute(vector<T*> atoms_gradients_densities, vector<T*> atoms_Fi_densities,
vector<T*> electron_densities, unsigned int size_order)
    {
    unsigned int i, j, k;
    unsigned int count_orbitals;
    unsigned int sum_electrons;
    array<unsigned int, this->max_electrons> index;
    // Index of electrons positions forcomputing nuclear atraction integrals
    array<T, this->max_electrons> wavefunction_coefficients;
    array<int, this->max_electrons> bonding;
    array<int, this->max_electrons> spin_paired;
    array<T, this->max_electrons> spins;
    array<T, this->max_electrons> Z;
    array<T, this->max_electrons> exchange_array;
    array<T, this->max_electrons> correlation_array;
    array<T*, this->max_atoms> atoms_gradients_array;
    array<T*, this->max_atoms> atoms_Fi_array;
    array<T*, this->max_atoms> electron_array;
    bool restriction;
    
    restriction = true;
    sum_electrons = this->results.n.size();
    count_orbitals = 0;
    
    if (sum_electrons > this->max_electrons)
        return(-1);
        
    for (i = 0; i < atoms_gradients_densities.size(); i++)
        atoms_gradients_array[i] = atoms_gradients_densities[i];
    
    for (i = 0; i < atoms_Fi_densities.size(); i++)
        atoms_Fi_array[i] = atoms_Fi_densities[i];
        
    for (i = 0; i < electron_densities.size(); i++)
        electron_array[i] = electron_densities[i];
    
    for (i = 0; i < sum_electrons; i++)
        {
        if (this->results.wavefunction_coefficients[i] != 0)
            {
            this->results.wavefunction_coefficients[i] = 1; // Cancelling semi-empirical xpression of correlation energy
            wavefunction_coefficients[i] = 1;
            }
        else
            wavefunction_coefficients[i] = 0;
        
        spins[i] = this->results.spins[i];
        spin_paired[i] = this->results.spin_paired[i];
        bonding[i] = this->results.bonding[i];
        Z[i] = this->results.Z[i];
        }
    // closed-shell Hartree-Fock method optimalization code
    for (i = 0; i < sum_electrons ; i++) // for restricted and unrestricted Hartree-Fock method
        {
        if (spin_paired[i] == -1)
            {
            restriction = false;
            }
        }
    for (i = 0; i < sum_electrons ; i++)
        if (wavefunction_coefficients[i] != 0)
            {
            if (i == 0)
                {
                index[count_orbitals] = i;
                count_orbitals++;
                }
            if (i > 0)
                {
                if ((spins[i] == 0.5) or (bonding[i] >= 0) or restriction == false)
                    {
                    index[count_orbitals] = i;
                    count_orbitals++;
                    }
                }
            }
    i = 0;
    j = 1;
    // multithreading code
    #pragma omp parallel
        {
        #pragma omp for
        for (i = 0; i < count_orbitals; i++)
            {
            PBE_thread(atoms_gradients_array[electron_to_atom_numbers[index[i]] - 1],
            atoms_Fi_array[electron_to_atom_numbers[index[i]] - 1], electron_array[index[i]],
            size_order, &exchange_array[index[i]], &correlation_array[index[i]], Wigner_Seitz_radiuses[Z[index[i]] - 1]);
            }
        }
    // End of multithreading code
    // closed-shell Hartree-Fock method optimalization code
    if (restriction == true) // copying values in restricted basis_set_method
        for (i = 0; i < sum_electrons; i++)
            {
            if (spins[i] == -0.5 and bonding[i] == -1 and spin_paired[i] >= 0)
                {
                exchange_array[i] = exchange_array[spin_paired[i]];
                correlation_array[i] = correlation_array[spin_paired[i]];
                }
            }
    // end closed-shell Hartree-Fock method optimalization code
    for (i = 0; i < sum_electrons; i++) // including the wavefunctions coefficients for linear combination
        {
        exchange_array[i] = exchange_array[i] * wavefunction_coefficients[i] * wavefunction_coefficients[i];
        correlation_array[i] = correlation_array[i] * wavefunction_coefficients[i] * wavefunction_coefficients[i];
        }
    exchange_energies.clear();
    correlation_energies.clear();
    for (i = 0; i < sum_electrons; i++) // copy values into vectors
        {
        exchange_energies.push_back(exchange_array[i]);
        correlation_energies.push_back(-correlation_array[i]* this->Hartree_energy_constant);
        }
    return(0);
    }
// VQE section
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Set_circuit()
    {
    unsigned int i, j;
    unsigned int k, l;
    unsigned int sum_electrons;
    unsigned int previous;
    
    sum_electrons = this->results.n.size();
    for (i = 0; i < count_atoms; i++) // Adding intraaratom entanglements
        {
        if (i < count_atoms -1)
            for (j = index_atoms[i] + 1; j < index_atoms[i + 1]; j++)
                VQE_gates.push_back(make_pair(j - 1, j));
        else
            for (j = index_atoms[i] + 1; j < sum_electrons; j++)
                VQE_gates.push_back(make_pair(j - 1, j));
        }
    for (i = 0; i < count_atoms; i++) // Adding interaratom entanglements
        for (j = i; j < count_atoms; j++)
            {
            if (sqrt(((this->results.x[index_atoms[i]] - this->results.x[index_atoms[j]]) *
                      (this->results.x[index_atoms[i]] - this->results.x[index_atoms[i]])) + 
                     ((this->results.y[index_atoms[i]] - this->results.y[index_atoms[j]]) *
                      (this->results.y[index_atoms[i]] - this->results.y[index_atoms[i]])) +
                     ((this->results.z[index_atoms[i]] - this->results.z[index_atoms[j]]) *
                      (this->results.x[index_atoms[i]] - this->results.z[index_atoms[i]]))) <=
                       VQE_interconnection_distance)
                {
                if (i < count_atoms - 1)
                    k = index_atoms[i + 1] - 1;
                else
                    k = sum_electrons - 1;
                
                if (j < count_atoms - 1)
                    l = index_atoms[j + 1] - 1;
                else
                    l = sum_electrons - 1;
                
                if (k != l)
                    VQE_gates.push_back(make_pair(k, l));
                }
            }
    // Qiskit circuit code
    circuit.set_registers(sum_electrons, sum_electrons);
    for (i = 0; i < sum_electrons; i++) // Apply Hadamard VQE_gates
        circuit.x(i);
    
    for (i = 0; i < VQE_gates.size(); i++) // Apply entanglements
        circuit.cx(VQE_gates[i].first, VQE_gates[i].second); // x-gates for Jordan-Wigner transformation
    
    for (i = 0; i < sum_electrons; i++) // Apply measure gates
        {
        circuit.ry(M_PI/2.00, i);
        circuit.measure(i, i);
        } // End of qiskit circuit code
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Run_circuit(vector<double>& ansatz)
    {
    unsigned int i, j, k;
    unsigned int sum_electrons;
    unsigned int pow;
    double state;
    vector<double> ket;
    vector<double> measured;
    vector<vector<double>> counts;
    // Qiskit circuit code
    QuantumCircuit q_circuit;
    
    sum_electrons = this->results.n.size();
    measured.resize(ansatz.size());
    q_circuit.set_registers(sum_electrons, sum_electrons);
    for (i = 0; i < sum_electrons; i++) // Applying the ansatz
        q_circuit.ry((ansatz[i] - 0.5) * M_PI, i);
    
    q_circuit.add(circuit); // Completing the circuit
    for (i = 0; i < VQE_shots; i++)
        {
        Simulator result(q_circuit);
        ket = result.get_probs(result.qc);
        for (j = 0; j < ansatz.size(); j++)
            {
            state = 0;
            pow = 1;
            for (k = 0; k < j; k++)
                pow = pow * 2;
            for (k = 0; k < ket.size(); k++)
                if (k/pow % 2 == 1)
                    state = state + ket[k];
            
            measured[j] = state;
            }
        counts.push_back(measured);
        if (circuit_string == "")
            circuit_string = result.get_qiskit();
        }
    for (i = 0; i < ansatz.size(); i++)
        {
        state = 0;
        for (j = 0; j < VQE_shots; j++)
            state = state + counts[j][i];
        
        ansatz[i] = state/VQE_shots;
        }
    // End of Qiskit circuit code
    return(0);
    }
// End of VQE section
// TPSS density functional section
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Allocate_gradients_2_densities(unsigned int size_order)
    {
    unsigned int i;
    unsigned int side = (2 * size_order) + 1;
    unsigned int size = side * side * side;
    unsigned int count_electrons = this->results.n.size();
    unsigned int count_orbitals = 0;
    
    array<int, this->max_electrons> bonding;
    array<int, this->max_electrons> spin_paired;
    array<T, this->max_electrons> spins;
    array<T, this->max_electrons> wavefunction_coefficients;
    bool restriction = true;
    T* density_pointer = nullptr;
    
    for (i = 0; i < count_electrons; i++)
        {
        spins[i] = this->results.spins[i];
        spin_paired[i] = this->results.spin_paired[i];
        bonding[i] = this->results.bonding[i];
        wavefunction_coefficients[i] = this->results.wavefunction_coefficients[i];
        }
    // closed-shell Hartree-Fock method optimalization code
    for (i = 0; i < count_electrons ; i++) // for restricted and unrestricted Hartree-Fock method
        {
        if (spin_paired[i] == -1)
            {
            restriction = false;
            }
        }
    // closed-shell Hartree-Fock method optimalization code
    // Allocate kinetic energy densities of orbitals
    for (i = 0; i < count_electrons; i++)
        try {
            if (restriction == false or bonding[i] >=0 or spin_paired[i] > i)
                {
                density_pointer = new T[size];
                electrons_gradients_2_densities.push_back(density_pointer);
                }
            else
                electrons_gradients_2_densities.push_back(electrons_gradients_2_densities[spin_paired[i]]);
            }
        catch (int)
            {
            return(-1);
            }
    // Allocate sums electrons kinetic energies densities of atoms
    for (i = 0; i < count_atoms; i++)
        try {
            density_pointer = new T[size];
            atoms_sum_electrons_gradients_2_densities.push_back(density_pointer);
            }
        catch (int)
            {
            return(-1);
            }
    // Allocate kinetic energy densities of atoms
    for (i = 0; i < count_atoms; i++)
        try {
            density_pointer = new T[size];
            atoms_gradients_2_densities.push_back(density_pointer);
            }
        catch (int)
            {
            return(-1);
            }
    // Allocate spin gradients
    for (i = 0; i < count_atoms; i++)
        try {
            density_pointer = new T[size];
            atoms_spins_gradients.push_back(density_pointer);
            }
        catch (int)
            {
            return(-1);
            }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::TPSS_thread(T* atoms_gradients_density, T* atoms_Fi_density,
T* electron_density, unsigned int size_order, T* exchange_energy, T* correlation_energy, T Wigner_Seitz_radius,
T* atom_electron_density, T* electron_gradients_2_density, T* atom_sum_gradients_2_densities,
T* atom_gradients_2_density, T* atoms_spin_density, T* atoms_spin_gradients_density, T spin)
    { // the same code, as PBE_thread, but extended to use kinetic nergy term
    unsigned int i, j;
    unsigned int size;
    T exchange_integral;
    T correlation_integral;
    T n;
    T k_F;
    T s;
    T pre_s;
    T c; // this is not sppeed of light in this context
    T X;
    T e_x_unif;
    T F_x;
    T nl_integral_factor;
    T s_2;
    T fi_2;
    T power;
    // Changes from PBE density functional
    T X_TPSS;
    T C_TPSS;
    
    T c_TPSS = 1.59096;
    T e_TPSS = 1.537;
    T z_TPSS = 1;
    T q_b_TPSS = 0.4;
    T factor_1_TPSS, factor_2_TPSS, factor_3_TPSS, factor_4_TPSS;
    T factor_5_TPSS, factor_6_TPSS, factor_7_TPSS, factor_8_TPSS;
    
    T pre_Ksi = 1/pow(6 * this->Pi * this->Pi, 2.00/3.00);
    T Dzeta;
    T Ksi;
    T tau_w_tau_2;
    T n_spins_n;
    
    T correlation_increment = 0;
    T exchange_increment = 0;
    
    factor_1_TPSS = 10.00/81.00 + c_TPSS * (z_TPSS * z_TPSS)/((z_TPSS * z_TPSS + 1) * (z_TPSS * z_TPSS + 1));
    factor_2_TPSS = 146.00/2025.00 * q_b_TPSS * q_b_TPSS;
    factor_3_TPSS = -73.00/405.00 * q_b_TPSS * sqrt(0.5);
    factor_4_TPSS = (0.6 * z_TPSS) * (0.6 * z_TPSS);
    factor_5_TPSS = (1.00/kappa) * (10.00/81.00) * (10.00/81.00);
    factor_6_TPSS = 2.00 * sqrt(e_TPSS) * 10.00/81.00 * factor_4_TPSS;
    factor_7_TPSS = e_TPSS * mi;
    factor_8_TPSS = (1 + sqrt(e_TPSS)) * (1 + sqrt(e_TPSS));
    // End changes from PBE density functional
    
    size = (2 * size_order + 1) * (2 * size_order + 1) * (2 * size_order + 1);
    // Correlation energy computation preparation
    n = (3 * this->Hartree_lenght * this->Hartree_lenght * this->Hartree_lenght)
    /(4 * this->Pi * Wigner_Seitz_radius * Wigner_Seitz_radius * Wigner_Seitz_radius);
    k_F = pow(n * 3 * this->Pi * this->Pi, 1.00/3.00);
    c = pow((3 * this->Pi * this->Pi)/16, 1.00/3.00);
    pre_s = 1.00/(2 * k_F);
    X = (beta/gamma) * c * c * exp(-omega/gamma); // equal 0.72161
    // Exgchange energy computation preparation
    e_x_unif = -(3 * this->e * k_F)/(4 * this->Pi);
    
    exchange_integral = 0;
    correlation_integral = 0;
    for (i = 0; i < size; i++) // Computing corelation_energy
        {
        if (electron_density[i] != 0)
            {
            s = pre_s * abs(atoms_gradients_density[i])/n;
            s_2 = s * s;
            fi_2 = atoms_Fi_density[i] * atoms_Fi_density[i];
            if (fi_2 != 0)
                power = pow(X, T(s_2/fi_2));
            nl_integral_factor = log(1 + 1/(power + (power * power)));
            if ((not (isnan(nl_integral_factor))) and (not (isinf(nl_integral_factor))))
                {
                // changes from PBE
                // computing C
                Ksi = pre_Ksi * atoms_spin_gradients_density[i] * pow(electron_density[i], - 1.00/3.00);
                if (atoms_spin_density[i]/atom_electron_density[i] < 1)
                    Dzeta = atoms_spin_density[i]/atom_electron_density[i];
                if (atom_sum_gradients_2_densities[i] > 0)
                    tau_w_tau_2 = (atom_gradients_2_density[i]/atom_sum_gradients_2_densities[i]) *
                    (atom_gradients_2_density[i]/atom_sum_gradients_2_densities[i]);
                n_spins_n = Dzeta * spin + 0.5;
                C_TPSS = 0.53 + (0.87 * Dzeta) + (0.5 * Dzeta * Dzeta) + (2.26 * Dzeta * Dzeta * Dzeta);
                C_TPSS = C_TPSS/pow(1 + Ksi * Ksi * 0.5 * (pow(1 + Ksi, -1.333) + pow(1 - Ksi, -1.333)), 4);
                
                correlation_increment = (gamma * electron_density[i] * fi_2 * atoms_Fi_density[i] * nl_integral_factor)
                * (1 + tau_w_tau_2 * (C_TPSS * (1 - n_spins_n) - n_spins_n));
                if (not isnan(correlation_increment))
                    correlation_integral = correlation_integral + correlation_increment;
                // End changes from PBE density functional
                }
            // Changes from PBE density functional
            // computing X
            X_TPSS = ((factor_1_TPSS * s_2) + factor_2_TPSS + (factor_3_TPSS * sqrt(factor_4_TPSS + s_2 * s_2)) +
            (factor_5_TPSS * s_2 * s_2) + factor_6_TPSS + (factor_7_TPSS * s_2 * s_2 * s_2))/(factor_8_TPSS * s_2 * s_2);
            exchange_increment = e_x_unif * (1 + kappa - kappa)/(1 + (X_TPSS/kappa)) * electron_density[i];
            if (not isnan(exchange_increment))
                exchange_integral = exchange_integral + exchange_increment;
            // End changes from PBE density functional
            }
        }
    exchange_energy[0] = exchange_integral;
    correlation_energy[0] = correlation_integral;
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::TPSS_compute(vector<T*> atoms_gradients_densities,
vector<T*> atoms_Fi_densities, vector<T*> electron_densities, unsigned int size_order)
    {
    unsigned int i, j, k;
    unsigned int count_orbitals;
    unsigned int sum_electrons;
    
    array<unsigned int, this->max_electrons> index;
    // Index of electrons positions forcomputing nuclear atraction integrals
    array<T, this->max_electrons> wavefunction_coefficients;
    array<int, this->max_electrons> bonding;
    array<int, this->max_electrons> spin_paired;
    array<T, this->max_electrons> spins;
    array<unsigned int, this->max_electrons> Z;
    array<T, this->max_electrons> exchange_array;
    array<T, this->max_electrons> correlation_array;
    array<T*, this->max_atoms> atoms_gradients_array;
    array<T*, this->max_atoms> atoms_Fi_array;
    array<T*, this->max_electrons> electron_array;
    
    array<T*, this->max_atoms> atoms_electron_densities_array;
    array<T*, this->max_electrons> electron_2_gradient_density_array;
    array<T*, this->max_atoms> atom_sum_gradients_2_densities_array;
    array<T*, this->max_atoms> atom_gradients_2_density_array;
    array<T*, this->max_atoms> atoms_spin_density_array;
    array<T*, this->max_atoms> atoms_spin_gradients_density;
    
    bool restriction;
    
    restriction = true;
    sum_electrons = this->results.n.size();
    count_orbitals = 0;
    
    if (sum_electrons > this->max_electrons)
        return(-1);
        
    for (i = 0; i < atoms_electron_densities.size(); i++)
        atoms_electron_densities_array[i] = atoms_electron_densities[i];
        
    for (i = 0; i < atoms_gradients_densities.size(); i++)
        atoms_gradients_array[i] = atoms_gradients_densities[i];
    
    for (i = 0; i < atoms_Fi_densities.size(); i++)
        atoms_Fi_array[i] = atoms_Fi_densities[i];
        
    for (i = 0; i < electron_densities.size(); i++)
        electron_array[i] = electron_densities[i];
    
    for (i = 0; i < sum_electrons; i++)
        {
        wavefunction_coefficients[i] = this->results.wavefunction_coefficients[i];
        spins[i] = this->results.spins[i];
        spin_paired[i] = this->results.spin_paired[i];
        bonding[i] = this->results.bonding[i];
        Z[i] = this->results.Z[i];
        
        electron_2_gradient_density_array[i] = electrons_gradients_2_densities[i];
        atom_sum_gradients_2_densities_array[i] = atoms_sum_electrons_gradients_2_densities[i];
        atom_gradients_2_density_array[i] = atoms_gradients_2_densities[i];
        atoms_spin_density_array[i] = atoms_spin_densities[i];
        atoms_spin_gradients_density[i] = atoms_spins_gradients[i];
        }
    // closed-shell Hartree-Fock method optimalization code
    for (i = 0; i < sum_electrons ; i++) // for restricted and unrestricted Hartree-Fock method
        {
        if (spin_paired[i] == -1 or bonding[i] == -1)
            {
            restriction = false;
            }
        }
    for (i = 0; i < sum_electrons ; i++)
        if (wavefunction_coefficients[i] != 0)
            {
            if (i == 0)
                {
                index[count_orbitals] = i;
                count_orbitals++;
                }
            if (i > 0)
                {
                if ((spins[i] == 0.5) or (bonding[i] >= 0) or restriction == false)
                    {
                    index[count_orbitals] = i;
                    count_orbitals++;
                    }
                }
            }
    i = 0;
    j = 1;
    
    // multithreading code
    #pragma omp parallel
        {
        #pragma omp for
        for (i = 0; i < count_orbitals; i++)
            {
            TPSS_thread(atoms_gradients_array[electron_to_atom_numbers[index[i]] - 1],
            atoms_Fi_array[electron_to_atom_numbers[index[i]] - 1], electron_array[index[i]],
            size_order, &exchange_array[index[i]], &correlation_array[index[i]], Wigner_Seitz_radiuses[Z[index[i]] - 1],
            atoms_electron_densities_array[electron_to_atom_numbers[index[i]] - 1],
            electron_2_gradient_density_array[index[i]],
            atom_sum_gradients_2_densities_array[electron_to_atom_numbers[index[i]] - 1],
            atom_gradients_2_density_array[electron_to_atom_numbers[index[i]] - 1],
            atoms_spin_density_array[electron_to_atom_numbers[index[i]] - 1],
            atoms_spin_gradients_density[electron_to_atom_numbers[index[i]] - 1],
            spins[index[i]]);
            }
        }
    // End of multithreading code
    // closed-shell Hartree-Fock method optimalization code
    if (restriction == true) // copying values in restricted basis_set_method
        for (i = 0; i < sum_electrons; i++)
            {
            if (spins[i] == -0.5 and bonding[i] == -1 and spin_paired[i] >= 0)
                {
                exchange_array[i] = exchange_array[spin_paired[i]];
                correlation_array[i] = correlation_array[spin_paired[i]];
                }
            }
    // end closed-shell Hartree-Fock method optimalization code
    for (i = 0; i < sum_electrons; i++) // including the wavefunctions coefficients for linear combination
        {
        exchange_array[i] = exchange_array[i] * wavefunction_coefficients[i] * wavefunction_coefficients[i];
        correlation_array[i] = correlation_array[i] * wavefunction_coefficients[i] * wavefunction_coefficients[i];
        }
    exchange_energies.clear();
    correlation_energies.clear();
    for (i = 0; i < sum_electrons; i++) // copy values into vectors
        {
        exchange_energies.push_back(exchange_array[i]);
        correlation_energies.push_back(-correlation_array[i]* this->Hartree_energy_constant);
        }
    return(0);
    }
// End of TPSS density functional section
template <typename T>
int Slater_basis_set_calculations_DFT<T>::String_to_advanced_parameters(string UI_input, unsigned int size_order,
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
    array<string, this->max_atoms * 3> central_cations_string;
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
    if (this->String_to_list_electrons(UI_input, size_order, extern_coordinates, x_2, y_2, z_2) != 0)
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
T Slater_basis_set_calculations_DFT<T>::Execute_calculation(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool deallocate, vector<T>* values)
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
    
    Hamiltonian = this->Calculate(max_iterations, minimal_fidelity, size_order, allocation_memory, values);
    if (Hamiltonian != -1)
        allocation_memory = false;
    if (deallocate == true) // deallocate atom electron and spin densities
        Clear();
    return(Hamiltonian);
    }
// Density functional theory section
template <typename T>
T Slater_basis_set_calculations_DFT<T>::Execute_PBE(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool deallocate, vector<T>* values)
    {
    unsigned int i, j, k;
    unsigned int matrix_order;
    T exchange_energy;
    T Hamiltonian;
    T electron_energy;
    vector<T> correction_diagonal;
    
    energies = values;
    if (allocation_memory == true)
        {
        if (Create_index_atoms() == -1)
           return(-1);
        if (Detect_center_of_mass() == -1)
           return(-1);
        if (Detect_dipole_moments(&dipoles) == -1)
           return(-1);
        if (central_cation_atoms.atom_numbers.size() > 0)
            Create_crystal_field(&central_cation_atoms);
        }
    matrix_order = this->results.n.size();
    
    for (i = 0; i < matrix_order; i++)
        correction_diagonal.push_back(this->correction_matrix[i * (matrix_order + 1)]);
    
    Hamiltonian = this->Calculate(max_iterations, minimal_fidelity, size_order, allocation_memory, values);
    if (Hamiltonian != -1)
        allocation_memory = false;
    else
        return(-1);
        
    for (i = 0; i < matrix_order; i++)
        pre_PBE_wavefunction_lenght_multipliers.push_back(this->results.wavefunction_lenght_multipliers[i]);
    // Execute calculations
    for (i = 0; i < matrix_order * matrix_order; i++) // avoiding executing PBE for atomic systems
        {
        if (this->resonance_integral_matrix[i] != 0 or this->Helium_correlation_energy == false)
            break;
        
        if (i == matrix_order * matrix_order -1)
            return(Hamiltonian);
        }
    if (allocation_PBE == true)
        if (Allocate_densities(size_order) != 0)
            return(-1);
    allocation_PBE = false;
    // compute electron and spin densities
    Compute_densities(this->results.probabilities, atoms_electron_densities,
    index_atoms, size_order, this->results.spins, this->results.spin_paired,
    this->results.x_range, this->results.y_range, this->results.z_range, false);
    Compute_densities(this->results.probabilities, atoms_spin_densities,
    index_atoms, size_order, this->results.spins, this->results.spin_paired,
    this->results.x_range, this->results.y_range, this->results.z_range, true);
    
    // correct electron and spin densities
    Correct_densities(atoms_electron_densities, this->results.x, this->results.y, this->results.z, size_order);
    Correct_densities(atoms_spin_densities, this->results.x, this->results.y, this->results.z, size_order);
    correlation_energies.clear();
    exchange_energies.clear();
    electron_energies.clear();
    // run PBE method over electron and spin densities
    Compute_Fi_and_gradients(atoms_electron_densities, atoms_spin_densities, atoms_Fi_densities, atoms_gradients_densities,
    size_order);
    PBE_compute(atoms_gradients_densities, atoms_Fi_densities, this->results.probabilities, size_order);
    
    for (i = 0; i < PBE_iterations; i++) 
        {
        // repeat calculation with computed electron correlation energies
        for (j = 0; j < matrix_order; j++)
            this->correction_matrix[j * (matrix_order + 1)] = correlation_energies[j] + exchange_energies[j]
            + correction_diagonal[j]; // PBE density functional not include polar part of bonds
        // Execute HF calculations
        Hamiltonian = this->Calculate(max_iterations, minimal_fidelity, size_order, false, values);
        matrix_order = this->results.n.size();
    
        // compute electron and spin densities
        Compute_densities(this->results.probabilities, atoms_electron_densities,
        index_atoms, size_order, this->results.spins, this->results.spin_paired,
        this->results.x_range, this->results.y_range, this->results.z_range, false);
        Compute_densities(this->results.probabilities, atoms_spin_densities,
        index_atoms, size_order, this->results.spins, this->results.spin_paired,
        this->results.x_range, this->results.y_range, this->results.z_range, true);
        // correct electron and spin densities
        Correct_densities(atoms_electron_densities, this->results.x, this->results.y, this->results.z, size_order);
        Correct_densities(atoms_spin_densities, this->results.x, this->results.y, this->results.z, size_order);
        correlation_energies.clear();
        exchange_energies.clear();
        electron_energies.clear();
        // run PBE method over electron and spin densities
        Compute_Fi_and_gradients(atoms_electron_densities, atoms_spin_densities, atoms_Fi_densities, atoms_gradients_densities,
        size_order);
        PBE_compute(atoms_gradients_densities, atoms_Fi_densities, this->results.probabilities, size_order);
        }
    if (deallocate == true) // deallocate atom electron aand spin densities
        Clear();
        
    return(Hamiltonian);
    }
// end of density functional theory section
// VQE section
template <typename T>
T Slater_basis_set_calculations_DFT<T>::Execute_PBE_VQE(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool deallocate, vector<T>* values)
    {
    unsigned int i, j;
    unsigned int matrix_order;
    bool VQE;
    T VQE_ansatz_sum;
    vector<T> previous_correlation_energies;
    vector<T> previous_exchange_energies;
    T max_PBE_difference = 0;
    
    VQE = true;
    kappa = 0.804; // homogenic electron gass
    mi = (beta * this->Pi * this->Pi)/3.00;
    VQE_Hamiltonian = Execute_PBE(max_iterations, minimal_fidelity, size_order, false, values);
    if (VQE_Hamiltonian == -1)
        return(-1);
    
    matrix_order = this->results.wavefunction_lenght_multipliers.size();
    
    for (i = 0; i < matrix_order * matrix_order; i++) // avoiding executing PBE and VQE for atomic systems
        {
        if (this->resonance_integral_matrix[i] != 0 or this->Helium_correlation_energy == false)
            break;
        
        if (i == matrix_order * matrix_order -1)
            return(VQE_Hamiltonian);
        }
    for (i = 0; i < matrix_order; i++) // Filling first wavefunction_lenghts_multipliers and eigenvectors vectors
        {
        previous_correlation_energies.push_back(correlation_energies[i]);
        previous_exchange_energies.push_back(exchange_energies[i]);
        }
    
    VQE_Wavefunction_lenght_multipliers_1 = this->results.wavefunction_lenght_multipliers;
    VQE_Eigenvectors_1 = values[0];
    
    for (i = 0; i < matrix_order; i++)
        this->results.wavefunction_lenght_multipliers[i] = pre_PBE_wavefunction_lenght_multipliers[i];
    
    kappa = 0.967; // Lieb-Oxford bonding
    mi = 0.235;
    VQE_previous_Hamiltonian = VQE_Hamiltonian;
    VQE_Hamiltonian = Execute_PBE(max_iterations, minimal_fidelity, size_order, false, values);
    if (VQE_Hamiltonian == -1)
        return(-1);
    
    for (i = 0; i < matrix_order; i++) // Filling second wavefunction_lenghts_multipliers and eigenvectors vectors
        {
        VQE_Wavefunction_lenght_multipliers_2.push_back(this->results.wavefunction_lenght_multipliers[i]);
        VQE_Eigenvectors_2.push_back(values->operator[](i));
        if (abs(VQE_Eigenvectors_2[i] - VQE_Eigenvectors_1[i]) > abs(max_PBE_difference))
            max_PBE_difference = VQE_Eigenvectors_2[i] - VQE_Eigenvectors_1[i];
        }
    for (i = 0; i < matrix_order; i++)
        {
        if (VQE_Eigenvectors_1 < VQE_Eigenvectors_2)
            VQE_correlation_energy_sign.push_back(true);
        else
            VQE_correlation_energy_sign.push_back(false);
        }
    Set_circuit();
    if (VQE_gates.size() > VQE_max_connections)
        VQE = false;
    
    for (i = 0; i < matrix_order; i++)
            {
            if (VQE_correlation_energy_sign[i] == true) // Second model - more lower energies
                VQE_ansatz.push_back(0.5 + 0.5 * (VQE_Eigenvectors_2[i] - VQE_Eigenvectors_1[i])/max_PBE_difference);
                
            else // First model - more lower energies
                VQE_ansatz.push_back(0.5 - 0.5 * (VQE_Eigenvectors_1[i] - VQE_Eigenvectors_2[i])/max_PBE_difference);
            }
    VQE_ansatz_sum = 0;
    for (j = 0; j < matrix_order; j++)
         VQE_ansatz_sum = VQE_ansatz_sum + VQE_ansatz[j];
    
    if (VQE == true and (VQE_ansatz_sum > 0.05 * matrix_order and VQE_ansatz_sum < 0.95 * matrix_order))
        {
        for (i = 0; i < VQE_max_iterations; i++)
            {
            Run_circuit(VQE_ansatz); // Computing through quantum circuit
            for (j = 0; j < matrix_order; j++)
                {
                if (VQE_correlation_energy_sign[j] == false) // First model - more higher energies
                    this->results.wavefunction_lenght_multipliers[j] = VQE_Wavefunction_lenght_multipliers_1[j] * (1 - VQE_ansatz[j]) +
                    VQE_Wavefunction_lenght_multipliers_2[j] * VQE_ansatz[j];
                else // Second model - more higher energies
                    this->results.wavefunction_lenght_multipliers[j] = VQE_Wavefunction_lenght_multipliers_2[j] * (1 - VQE_ansatz[j]) +
                    VQE_Wavefunction_lenght_multipliers_1[j] * VQE_ansatz[j];
                }
            // linear combination of homogenous electron gas and Lieb-Oxford bonding
            kappa = 0.804 * (1 - VQE_ansatz_sum/T(matrix_order)) + 0.967 * (VQE_ansatz_sum/T(matrix_order)); 
            mi = (1 - log(2))/(this->Pi * this->Pi) * (1 - VQE_ansatz_sum/T(matrix_order)) + 0.235 * (VQE_ansatz_sum/T(matrix_order));
            for (i = 0; i < matrix_order; i++)
                this->results.wavefunction_lenght_multipliers[i] = pre_PBE_wavefunction_lenght_multipliers[i];
            
            VQE_previous_Hamiltonian = VQE_Hamiltonian;
            VQE_Hamiltonian = Execute_PBE(max_iterations, minimal_fidelity, size_order, false, values);
            if (VQE_Hamiltonian == -1)
                return(-1);
            if ((abs(VQE_Hamiltonian/VQE_previous_Hamiltonian) < (2 - minimal_fidelity)) and
                (abs(VQE_Hamiltonian/VQE_previous_Hamiltonian) > minimal_fidelity) and (VQE_Hamiltonian/VQE_previous_Hamiltonian > 0))
                break;
            
            max_PBE_difference = 0;
            for (i = 0; i < matrix_order; i++) // Calculating new ansatz vector
                {
                if (abs(values->operator[](i) - VQE_Eigenvectors_1[i]) > abs(max_PBE_difference))
                    max_PBE_difference = values->operator[](i) - VQE_Eigenvectors_1[i];
                }
            for (i = 0; i < matrix_order; i++) // Calculating new ansatz vector
                {
                if (VQE_correlation_energy_sign[i] == true) // Second model - more lower energies
                    VQE_ansatz[i] = (0.5 + 0.5 * (values->operator[](i) - VQE_Eigenvectors_1[i])/max_PBE_difference);
                else // First model - more lower energies
                    VQE_ansatz[i] = (0.5 - 0.5 * (VQE_Eigenvectors_1[i] - values->operator[](i))/max_PBE_difference);
                }
            VQE_iterations = VQE_iterations + 1;
            }
        if (deallocate == true) // deallocate atom electron aand spin densities
            Clear();
        return VQE_Hamiltonian;
        }
    else
        {
        if (VQE_previous_Hamiltonian < VQE_Hamiltonian)
            {
            for (i = 0; i < matrix_order; i++) // Calculating new ansatz vector
                {
                this->results.wavefunction_lenght_multipliers[i] = VQE_Wavefunction_lenght_multipliers_1[i];
                
                correlation_energies[i] = previous_correlation_energies[i];
                exchange_energies[i] = previous_exchange_energies[i];
                }
            values[0] = VQE_Eigenvectors_1;
            }
        if (deallocate == true) // deallocate atom electron aand spin densities
            Clear();
        if (VQE_correlation_energy_sign[0] == false)
            return VQE_Hamiltonian;
        else
            return VQE_previous_Hamiltonian;
        }
    }
// end of VQE section
// TPSS section
template <typename T>
T Slater_basis_set_calculations_DFT<T>::Execute_TPSS(unsigned int max_iterations, T minimal_fidelity,
unsigned int size_order, bool deallocate, vector<T>* values)
    {
    unsigned int i, j;
    unsigned int matrix_order = this->results.n.size();
    T PBE_VQE_Hamiltonian;
    T last_Hamiltonian;
    T TPSS_Hamiltonian;
    vector<T> correction_diagonal;
    
    bool higher_density = false;
    
    for (i = 0; i < matrix_order; i++)
        correction_diagonal.push_back(this->correction_matrix[i * (matrix_order + 1)]);
    
    for (i = 0; i < matrix_order; i++)
        if (this->results.reduced_Z[i] > 10)
            higher_density = true;
    
    if (higher_density == true)
        this->Helium_correlation_energy = false;
    
    VQE_Hamiltonian = Execute_PBE_VQE(max_iterations, minimal_fidelity, size_order, false, values);
    if (VQE_Hamiltonian == -1)
        {
        TPSS_error = true;
        return(-1);
        }
    if (allocation_TPSS == true)
        if (Allocate_gradients_2_densities(size_order) == -1)
            {
            TPSS_error = true;
            return(PBE_VQE_Hamiltonian);
            }
        else
            allocation_TPSS = false;
    last_Hamiltonian = PBE_VQE_Hamiltonian;
    
    for (i = 0; i < matrix_order * matrix_order; i++) // avoiding executing PBE and VQE for atomic systems
        {
        if (this->resonance_integral_matrix[i] != 0 or this->Helium_correlation_energy == false)
            break;
        
        if (i == matrix_order * matrix_order -1)
            return(PBE_VQE_Hamiltonian);
        }
    
    for (i = 0; i < TPSS_iterations; i++)
        {
        // Prepare kinetic energies 3D rasters
        Compute_gradients(this->results.Gradients, electrons_gradients_2_densities, size_order);
        Compute_densities(electrons_gradients_2_densities, atoms_sum_electrons_gradients_2_densities,
        index_atoms, size_order, this->results.spins, this->results.spin_paired,
        this->results.x_range, this->results.y_range, this->results.z_range, false);
        Correct_densities(atoms_sum_electrons_gradients_2_densities,
        this->results.x, this->results.y, this->results.z, size_order);
        Compute_gradients(atoms_gradients_densities, atoms_gradients_2_densities, size_order);
        Compute_gradients(atoms_spin_densities, atoms_spins_gradients, size_order);
        // Compute exchange and correlation energies via TPSS
        TPSS_compute(atoms_gradients_densities, atoms_Fi_densities, this->results.probabilities, size_order);
        // compute Hamiltonian and energy levels
        for (j = 0; j < matrix_order; j++)
            this->correction_matrix[j * (matrix_order + 1)] = correlation_energies[j] + exchange_energies[j]
            + correction_diagonal[j]; // PBE density functional not include polar part of bonds
        
        TPSS_Hamiltonian = Execute_calculation(max_iterations, minimal_fidelity, size_order, deallocate, values);
        if (TPSS_Hamiltonian == -1)
            return(-1);
        if (i + 1 < TPSS_iterations)
            {
            // Prepare new atomic densities shared with PBE
            Compute_densities(this->results.probabilities,
            atoms_electron_densities, index_atoms, size_order, this->results.spins, this->results.spin_paired,
            this->results.x_range, this->results.y_range, this->results.z_range, false);
            Compute_densities(this->results.probabilities,
            atoms_spin_densities, index_atoms, size_order, this->results.spins, this->results.spin_paired,
            this->results.x_range, this->results.y_range, this->results.z_range, true);
            // correct electron and spin densities
            Correct_densities(atoms_electron_densities, this->results.x, this->results.y, this->results.z, size_order);
            Correct_densities(atoms_spin_densities, this->results.x, this->results.y, this->results.z, size_order);
            correlation_energies.clear();
            exchange_energies.clear();
            electron_energies.clear();
            // Compute Fi and gradients
            Compute_Fi_and_gradients(atoms_electron_densities, atoms_spin_densities, atoms_Fi_densities,
            atoms_gradients_densities, size_order);
            }
        }
    if (deallocate == true)
        Clear();
    return(TPSS_Hamiltonian);
    }
// End of TPSS section
// Orthonormalizing and Gaussian export section
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Compute_overlap_matrix_reduced()
    {
    unsigned int i, j, k;
    unsigned int min;
    unsigned int max;
    unsigned int matrix_order = this->results.n.size();
    
    // run computing a overlap matrix, if is not completely computed
    if (this->bonded_system == false)
        this->Create_overlap_integral_matrix(this->overlap_integral_matrix, matrix_order, &this->results);
    
    if (this->positive_intraatomic_overlap_matrix == nullptr)
        positive_intraatomic_overlap_matrix = new T[matrix_order * matrix_order];
    
    memset(positive_intraatomic_overlap_matrix, 0, matrix_order * matrix_order);
    // copy the overlaps with electrons on the same atom
    for (i = 0; i < index_atoms.size(); i++)
        {
        min = index_atoms[i];
        if (i + 1 < index_atoms.size())
            max = index_atoms[i + 1];
        else
            max = matrix_order -1;
        // vectorisation code
        for (j = min; j < max; j++)
            {
            for (k = min; k + 7 < max; k++)
                {
                positive_intraatomic_overlap_matrix[j + k * matrix_order] =
                abs(this->overlap_integral_matrix[j + k * matrix_order]);
                positive_intraatomic_overlap_matrix[j + k * matrix_order +1] =
                abs(this->overlap_integral_matrix[j + k * matrix_order +1]);
                positive_intraatomic_overlap_matrix[j + k * matrix_order +2] =
                abs(this->overlap_integral_matrix[j + k * matrix_order +2]);
                positive_intraatomic_overlap_matrix[j + k * matrix_order +3] =
                abs(this->overlap_integral_matrix[j + k * matrix_order +3]);
                positive_intraatomic_overlap_matrix[j + k * matrix_order +4] =
                abs(this->overlap_integral_matrix[j + k * matrix_order +4]);
                positive_intraatomic_overlap_matrix[j + k * matrix_order +5] =
                abs(this->overlap_integral_matrix[j + k * matrix_order +5]);
                positive_intraatomic_overlap_matrix[j + k * matrix_order +6] =
                abs(this->overlap_integral_matrix[j + k * matrix_order +6]);
                positive_intraatomic_overlap_matrix[j + k * matrix_order +7] =
                abs(this->overlap_integral_matrix[j + k * matrix_order +7]);
                }
            for (k = max - (max % 8); k < max; k++)
                {
                positive_intraatomic_overlap_matrix[j + k * matrix_order] =
                abs(this->overlap_integral_matrix[j + k * matrix_order]);
                }
            }
        // end of vectorisation code
        }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Compute_Hamiltonian()
    {
    unsigned int matrix_order = this->results.n.size();
    // recompute Hamiltonian
    this->Calculate_corr_basis_set_matrix(this->basis_set_matrix, this->correction_matrix, this->corr_basis_set_matrix,
    matrix_order);
    return(0);
    }
/*template <typename T>
int Slater_basis_set_calculations_DFT<T>::Create_orthonormalization_matrix()
    {
    unsigned int matrix_order = this->results.n.size();
     
    if (this->orthonormalizing_matrix == nullptr)
        orthonormalizing_matrix = new T[matrix_order * matrix_order];
    // Map the memory locations to Eigen matrices  
    Eigen::Map<Eigen::MatrixXd> A(positive_intraatomic_overlap_matrix, matrix_order, matrix_order);
    Eigen::Map<Eigen::MatrixXd> X(orthonormalizing_matrix, matrix_order, matrix_order);
     
    // Compute the square root of A  
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(A);
    solver.compute(A);
    if (solver.info() != Eigen::Success) {  
        std::cerr << "Error: Failed to compute eigen decomposition." << std::endl;  
        return(-1);
        }
    Eigen::VectorXd eigenvalues = solver.eigenvalues();  
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();
    // Compute the inverse square root of the eigenvalues by canonical orthogonalization
    Eigen::VectorXd invSqrtEigenvalues = Eigen::VectorXd::Zero(eigenvalues.size());  
    for (int i = 0; i < eigenvalues.size(); ++i)
        {  
        if (eigenvalues(i) > 1e-6) {  // Threshold for numerical stability  
            invSqrtEigenvalues(i) = 1.0 / sqrt(eigenvalues(i));  
        }
        }
    X = eigenvectors * invSqrtEigenvalues.asDiagonal();
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Orthonormalize_Hamiltonian()
    {
    unsigned int matrix_order = this->results.n.size();
    T Hamiltonian_magnitude;
    Eigen::VectorXd diagonal;
    Eigen::MatrixXd upper;
     
    if (this->orthonormalizing_matrix == nullptr)
        orthonormalizing_matrix = new T[matrix_order * matrix_order];
    // Map the memory locations to Eigen matrices  
    Eigen::Map<Eigen::MatrixXd> X(orthonormalizing_matrix, matrix_order, matrix_order);
    Eigen::Map<Eigen::MatrixXd> Hamiltonian(this->basis_set_matrix, matrix_order, matrix_order);
     
    Hamiltonian = Hamiltonian * X;
    // solve sum of diagonal and upper diagonal elements
    diagonal = Hamiltonian.diagonal();
    upper = Hamiltonian.triangularView<Eigen::Upper>();
    Hamiltonian_magnitude = diagonal.sum() + upper.sum();
    return(Hamiltonian_magnitude);
    }*/
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Gaussian_quadrature()
    {
    unsigned int i, j;
    unsigned int size;
    unsigned int nodes;
    T exponents_multiplier;
    T exponent;
    T coefficient;
    
    T* wavefunction_lenght_multipliers = this->results.wavefunction_lenght_multipliers.data();
    unsigned int* n = this->results.n.data();
    unsigned int* l = this->results.l.data();
    // Aro 90 basis setas analytical shape for Gaussian Quadrature
    T Aro_90_s1_exponents[] = {35.523221, 6.513144, 1.822143, 0.625955, 0.243077, 0.100112};
    T Aro_90_s1_coefficients[] = {0.009164, 0.049361, 0.168538, 0.370563, 0.416492, 0.130334};
    T Aro_90_s_exponents[] = {0.804524, 0.290399, 0.135201, 0.071840, 0.041270, 0.024899, 0.015522, 0.009899,
    0.006417, 0.004213, 0.002793, 0.001866, 0.001255, 0.000848, 0.000574, 0.000384, 0.000241, 0.000122, 0.000010};
    T Aro_90_p_exponents[] = {0.957502, 0.334869, 0.156031, 0.083822, 0.048802, 0.029849, 0.018858, 0.012192,
    0.008021, 0.005349, 0.003606, 0.002453, 0.001681, 0.001159, 0.000802, 0.000554, 0.000369, 0.000215, 0.000068};
    T Aro_90_d_exponents[] = {0.533897, 0.201903, 0.099729, 0.056230, 0.034095, 0.021582, 0.014046, 0.009324,
    0.006283, 0.004285, 0.002952, 0.002050, 0.001433, 0.001008, 0.000709, 0.000487, 0.000302, 0.000124};
    T Aro_90_f_exponents[] = {0.329197, 0.133383, 0.069204, 0.040519, 0.025300, 0.016391, 0.010876, 0.007342,
    0.005024, 0.003476, 0.002427, 0.001708, 0.001209, 0.000858, 0.000599, 0.000384, 0.000174};
    T Aro_90_g_exponents[] = {0.207044, 0.090784, 0.049515, 0.029994, 0.019180, 0.012651, 0.008515, 0.005820,
    0.004028, 0.002816, 0.001985, 0.001410, 0.001005, 0.000708, 0.000462, 0.000220};
    T Aro_90_h_exponents[] = {0.153191, 0.069181, 0.038612, 0.023758, 0.015343, 0.010183, 0.006885, 0.004724,
    0.003279, 0.002299, 0.001625, 0.001155, 0.000814, 0.000534, 0.000254};
    T Aro_90_i_exponents[] = {0.119421, 0.054816, 0.031022, 0.019236, 0.012456, 0.008269, 0.005587, 0.003829,
    0.002654, 0.001858, 0.001310, 0.000918, 0.000599, 0.000275};
    T Aro_90_k_exponents[] = {0.092378, 0.043392, 0.024899, 0.015512, 0.010037, 0.010037, 0.004474, 0.003054,
    0.002109, 0.001470, 0.001020, 0.000658, 0.000284};
    T Aro_90_l_exponents[] = {0.071875, 0.034537, 0.020008, 0.012454, 0.008016, 0.005271, 0.003525, 0.002389,
    0.001638, 0.001121, 0.000711, 0.000275};
    unsigned int count_nodes[] = {20, 20, 19, 18, 17, 16, 15, 14, 13, 12};
    // Resetting vector with Gaussians
    size = Gaussian_basis.size();
    for (i = 0; i < size; i++) {
        Gaussian_basis[i].clear();
        }
    if (size > 0)
        Gaussian_basis.clear();
    // Gaussian quadrature
    size = this->results.n.size();
    Gaussian_basis.resize(size);
    for (i = 0; i < size; i++)
        if (this->results.wavefunction_constraints[i] == 0) // check, if basis was changed
            {
            nodes = n[i] - l[i];
            exponents_multiplier = n[i] * wavefunction_lenght_multipliers[i];
            for (j = 0; j < 6; j++)
                {
                exponent = Aro_90_s1_exponents[j] * exponents_multiplier;
                coefficient = Aro_90_s1_coefficients[j]/exponents_multiplier;
                Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                }
            if (nodes > count_nodes[l[i]])
                return(-1);
                // check, if basis available
            for (j = 0; j + 1 < nodes; j++)
                {
                switch (l[i]){
                    case 0: {
                        exponent = Aro_90_s_exponents[j] * exponents_multiplier;
                        coefficient = 1/exponents_multiplier;
                        if (j % 2 == 0)
                            coefficient *= -1;
                        Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                        break;
                        }
                    case 1: {
                        exponent = Aro_90_p_exponents[j] * exponents_multiplier;
                        coefficient = 1/exponents_multiplier;
                        if (j % 2 == 0)
                            coefficient *= -1;
                        Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                        break;
                        }
                    case 2: {
                        exponent = Aro_90_d_exponents[j] * exponents_multiplier;
                        coefficient = 1/exponents_multiplier;
                    if (j % 2 == 0)
                        coefficient *= -1;
                        Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                        break;
                        }
                    case 3: {
                        exponent = Aro_90_f_exponents[j] * exponents_multiplier;
                        coefficient = 1/exponents_multiplier;
                        if (j % 2 == 0)
                            coefficient *= -1;
                        Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                        break;
                        }
                    case 4: {
                        exponent = Aro_90_g_exponents[j] * exponents_multiplier;
                        coefficient = 1/exponents_multiplier;
                        if (j % 2 == 0)
                            coefficient *= -1;
                        Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                        break;
                        }
                    case 5: {
                        exponent = Aro_90_g_exponents[j] * exponents_multiplier;
                        coefficient = 1/exponents_multiplier;
                        if (j % 2 == 0)
                            coefficient *= -1;
                        Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                        break;
                        }
                    case 6: {
                        exponent = Aro_90_h_exponents[j] * exponents_multiplier;
                        coefficient = 1/exponents_multiplier;
                        if (j % 2 == 0)
                            coefficient *= -1;
                        Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                        break;
                        }
                    case 7: {
                        exponent = Aro_90_i_exponents[j] * exponents_multiplier;
                        coefficient = 1/exponents_multiplier;
                        if (j % 2 == 0)
                            coefficient *= -1;
                        Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                        break;
                        }
                    case 8: {
                        exponent = Aro_90_k_exponents[j] * exponents_multiplier;
                        coefficient = 1/exponents_multiplier;
                        if (j % 2 == 0)
                            coefficient *= -1;
                        Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                        break;
                        }
                }
            }
        }
    return(0);
    }
// End of orthonormalizing and Gaussian export section

// CI and basis sets creating section
template <typename T>
T Slater_basis_set_calculations_DFT<T>::Execute_Basis_set_creation(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool deallocate, vector<T>* values, unsigned int count_shells,
unsigned int level_correlation_energy, vector<T> correction_energies)
    {
    unsigned int i, j, k, t;
    unsigned int current_atom_size;
    unsigned int n, n_ground;
    unsigned int l, l_ground;
    int m, m_ground;
    T s, s_ground;
    T s_list[] = {-0.5, 0.5};
    T Hamiltonian, Groung_Hamiltoninan;
    
    unsigned int atom_begin;
    unsigned int atom_end;
    unsigned int excited_position, extended_excited_position;
    bool found;
    
    T wavefunction_lenght_multiplier_ground;
    vector<T> nuclear_atraction_integral_matrix_row;
    vector<T> nucleuses_distances_row;
    vector<T> nucleuses_atractions_row;
    vector<T> coulombic_integral_matrix_row;
    vector<T> overlap_integral_matrix_row;
    vector<T> resonance_integral_matrix_row;
    vector<T> kinetic_integral_matrix_row;
    
    unsigned int count_electrons = this->results.n.size();
    // resize vectors
    if (count_shells == 0)
        return(-1);
    
    switch (count_shells)
        {
        case 1:
            {
            count_basis_per_atom = 2;
            break;
            }
        case 2:
            {
            count_basis_per_atom = 10;
            break;
            }
        case 3:
            {
            count_basis_per_atom = 28;
            break;
            }
        case 4:
            {
            count_basis_per_atom = 60;
            break;
            }
        case 5:
            {
            count_basis_per_atom = 110;
            break;
            }
        case 6:
            {
            count_basis_per_atom = 182;
            break;
            }
        case 7:
            {
            count_basis_per_atom = 280;
            break;
            }
        }
    if (Create_index_atoms() == -1)
        return(-1);
    count_basis_overall = count_basis_per_atom * count_atoms;
    
    
    // check, if selected count of basis is sufficient
    for (i = 0; i < count_electrons; i++)
        {
        if (this->results.n[i] > count_shells)
            {
            return(-1);
            }
        this->results.spin_paired[i] = -1;
        }
    
    // resize vectors of basis
    n_extended.resize(count_basis_overall);
    l_extended.resize(count_basis_overall);
    m_extended.resize(count_basis_overall);
    s_extended.resize(count_basis_overall);
    
    Slater_basis_exponents_extended.resize(count_basis_overall);
    Gaussian_basis_extended.resize(count_basis_overall);
    correlation_energies_extended.resize(count_basis_overall);
    
    // Run computation of the ground state
    switch (level_correlation_energy)
        {
        case 0:
            { // semiempirical correlation energy
            Hamiltonian = Execute_calculation(max_iterations, minimal_fidelity, size_order, false, values);
            if (Hamiltonian == -1)
                return(-1);
            break;
            }
        case 1:
            { // PBE density functional
            Hamiltonian = Execute_PBE(max_iterations, minimal_fidelity, size_order, false, values);
            if (Hamiltonian == -1)
                return(-1);
            break;
            }
        case 2:
            { // PBE density functional optimized by variational quantum eigensolver
            Hamiltonian = Execute_PBE_VQE(max_iterations, minimal_fidelity, size_order, false, values);
            if (Hamiltonian == -1)
                return(-1);
            break;
            }
        case 3:
            { // TPSS density functional
            Hamiltonian = Execute_TPSS(max_iterations, minimal_fidelity, size_order, false, values);
            if (Hamiltonian == -1)
                return(-1);
            break;
            }
        case 4:
            { // Externally correlated
            if (correction_energies.size() < count_atoms * count_basis_per_atom)
                return(-1);
            
            // Set correlation energies from vector
            this->Helium_correlation_energy = false;
            for (i = 0; i < count_atoms; i++)
                for (j = 0; j < counts_ground_basis[i]; j++)
                    {
                    this->correction_matrix[(index_atoms[i] + j) * (count_electrons + 1)]
                    = correction_energies[(i * count_basis_per_atom) + j];
                    }
            // copy vector
            for (i = 0; i < count_atoms * count_basis_per_atom; i++)
                correlation_energies_extended.push_back(correction_energies[i]);
            // Compute basis
            Hamiltonian = Execute_calculation(max_iterations, minimal_fidelity, size_order, false, values);
            if (Hamiltonian == -1)
                return(-1);
            break;
            }
        }
    // Export to basis
    for (i = 0; i < count_atoms; i++)
        {
        current_atom_size = this->results.count_electrons[index_atoms[i]];
        for (j = 0; j < current_atom_size; j++)
            {
            Slater_basis_exponents_extended[i * count_basis_per_atom + j] =
            this->results.Z[index_atoms[i] + j] * this->results.wavefunction_lenght_multipliers[index_atoms[i] + j];
            n_extended[i * count_basis_per_atom + j] = this->results.n[index_atoms[i] + j];
            l_extended[i * count_basis_per_atom + j] = this->results.l[index_atoms[i] + j];
            m_extended[i * count_basis_per_atom + j] = this->results.m[index_atoms[i] + j];
            s_extended[i * count_basis_per_atom + j] = this->results.spins[index_atoms[i] + j];
            }
        counts_ground_basis.push_back(current_atom_size);
        }
    Gaussian_quadrature();
    for (i = 0; i < count_atoms; i++)
        {
        current_atom_size = this->results.count_electrons[index_atoms[i]];
        for (j = 0; j < current_atom_size; j++)
            {
            Gaussian_basis_extended[i * count_basis_per_atom + j] = Gaussian_basis[index_atoms[i] + j];
            }
        }
    if (level_correlation_energy > 0 and level_correlation_energy < 3) // If a density functional is used
        {
        for (i = 0; i < count_atoms; i++)
            {
            current_atom_size = this->results.count_electrons[index_atoms[i]];
            for (j = 0; j < current_atom_size; j++)
                {
                correlation_energies_extended[i * count_basis_per_atom + j] =
                correlation_energies[index_atoms[i] + j] + exchange_energies[index_atoms[i] + j];
                }
            }
        }
    Groung_Hamiltoninan = Hamiltonian;
    // Create excitations and run computations of excited states
    
    // Froze all orbitals
    for (i = 0; i < count_electrons; i++)
        this->results.wavefunction_constraints[i] = 1;
    // resize vectors of integral matrices rows for restoring a ground state
    nuclear_atraction_integral_matrix_row.resize(count_electrons);
    nucleuses_distances_row.resize(count_electrons);
    nucleuses_atractions_row.resize(count_electrons);
    coulombic_integral_matrix_row.resize(count_electrons);
    overlap_integral_matrix_row.resize(count_electrons);
    resonance_integral_matrix_row.resize(count_electrons);
    kinetic_integral_matrix_row.resize(count_electrons);
    
    // Create excitations electrons
    for (i = 0; i < count_atoms; i++)
        {
        atom_begin = i * count_basis_per_atom;
        extended_excited_position = atom_begin + counts_ground_basis[i] - 1;
        excited_position = index_atoms[i] + counts_ground_basis[i] - 1;
        atom_end = extended_excited_position;
        // defroze orbital for excitation
        this->results.wavefunction_constraints[excited_position] = 0;
        // copy wavefunction_lenght_multiplies
        wavefunction_lenght_multiplier_ground = this->results.wavefunction_lenght_multipliers[excited_position];
        // copy integral matrixes rows
        for (j = 0; j < count_electrons; j++)
            nuclear_atraction_integral_matrix_row[j] =
            this->nuclear_atraction_integral_matrix[(excited_position * count_electrons)  + j];
        for (j = 0; j < count_electrons; j++)
            nucleuses_distances_row[j] = this->nucleuses_distances[(excited_position * count_electrons) + j];
        for (j = 0; j < count_electrons; j++)
            nucleuses_atractions_row[j] = this->nucleuses_atractions[(excited_position * count_electrons) + j];
        for (j = 0; j < count_electrons; j++)
            coulombic_integral_matrix_row[j] =  this->coulombic_integral_matrix[(excited_position * count_electrons) + j];
        for (j = 0; j < count_electrons; j++)
            overlap_integral_matrix_row[j] = this->overlap_integral_matrix[(excited_position * count_electrons) + j];
        for (j = 0; j < count_electrons; j++)
            resonance_integral_matrix_row[j] = this->resonance_integral_matrix[(excited_position * count_electrons) + j];
        for (j = 0; j < count_electrons; j++)
            kinetic_integral_matrix_row[j] =  this->kinetic_integral_matrix[(excited_position * count_electrons) + j];
        // compute excited states
        for (n = 1; n <= count_shells; n++)
            for (l = 0; l < n; l++)
                for (m = -int(l); m <= int(l); m++)
                    for (j = 0; j < 2; j++)
                        {
                        // check for previous quantum number configuration in atom
                        found = false;
                        for (k = atom_begin; k <= extended_excited_position; k++)
                            if (n_extended[k] == n and l_extended[k] == l and
                                m_extended[k] == m and s_extended[k] == s_list[j])
                                {
                                found = true;
                                break;
                                }
                        
                        if (found == false)
                            {
                            // create excitation;
                            this->Create_excitation(excited_position, n, l, m, s_list[j], true);
                            // recompute basis wavefunctions
                            switch (level_correlation_energy)
                                {
                                case 0:
                                    { // semiempirical correlation energy
                                    this->Helium_correlation_energy = true;
                                    Hamiltonian = Execute_calculation(max_iterations, minimal_fidelity,
                                    size_order, false, values);
                                    if (Hamiltonian == -1)
                                        return(-1);
                                    break;
                                    }
                                case 1:
                                    { // PBE density functional
                                    PBE_thread(atoms_gradients_densities[i], atoms_Fi_densities[i],
                                    this->results.probabilities[excited_position], size_order,
                                    &exchange_energies[excited_position], &correlation_energies[excited_position],
                                    Wigner_Seitz_radiuses[this->results.Z[excited_position] - 1]);

                                    this->correction_matrix[excited_position * (count_electrons + 1)] =
                                    correlation_energies[excited_position] + exchange_energies[excited_position];
                                    // Execute HF calculations
                                    Hamiltonian = Execute_PBE(max_iterations, minimal_fidelity,
                                    size_order, false, values);
                                    if (Hamiltonian == -1)
                                        return(-1);
                                    break;
                                    }
                                case 2:
                                    { // PBE density functional optimized by variational quantum eigensolver
                                    PBE_thread(atoms_gradients_densities[i], atoms_Fi_densities[i],
                                    this->results.probabilities[excited_position], size_order,
                                    &exchange_energies[excited_position], &correlation_energies[excited_position],
                                    Wigner_Seitz_radiuses[this->results.Z[excited_position] - 1]);

                                    this->correction_matrix[excited_position * (count_electrons + 1)] =
                                    correlation_energies[excited_position] + exchange_energies[excited_position];
                                    // Execute HF calculations
                                    Hamiltonian = Execute_PBE(max_iterations, minimal_fidelity,
                                    size_order, false, values);
                                    if (Hamiltonian == -1)
                                        return(-1);
                                    break;
                                    }
                                case 3:
                                    { // TPSS density functional
                                        
                                    TPSS_thread(atoms_gradients_densities[i], atoms_Fi_densities[i],
                                    this->results.probabilities[excited_position], size_order,
                                    &exchange_energies[excited_position], &correlation_energies[excited_position],
                                    Wigner_Seitz_radiuses[this->results.Z[excited_position] - 1],
                                    atoms_electron_densities[i], electrons_gradients_2_densities[excited_position],
                                    atoms_sum_electrons_gradients_2_densities[i], atoms_gradients_2_densities[i],
                                    atoms_spin_densities[i], atoms_spins_gradients[i], s_list[j]);
                                        
                                    this->correction_matrix[excited_position * (count_electrons + 1)] =
                                    correlation_energies[excited_position] + exchange_energies[excited_position];
                                    // Execute HF calculations
                                    Hamiltonian = Execute_PBE(max_iterations, minimal_fidelity,
                                    size_order, false, values);
                                    if (Hamiltonian == -1)
                                        return(-1);
                                    break;
                                    }
                                case 4:
                                    { // Externally correlated
                                    this->Helium_correlation_energy = false;
                                    if (correction_energies.size() < count_atoms * count_basis_per_atom)
                                        return(-1);
            
                                    this->correction_matrix[(index_atoms[i] + excited_position) *
                                    (count_electrons + 1)] = correction_energies[atom_end];
                                                
                                    // Compute basis
                                    Hamiltonian = Execute_calculation(max_iterations, minimal_fidelity, size_order,
                                    false, values);
                                    if (Hamiltonian == -1)
                                        return(-1);
                                    break;
                                    }
                                }
                            // export basis wavefunctions
                            Gaussian_quadrature();
                            // update position for saving basis
                            atom_end++;
                            // save computed basis
                            n_extended[atom_end] = n;
                            l_extended[atom_end] = l;
                            m_extended[atom_end] = m;
                            s_extended[atom_end] = s_list[j];
                            Slater_basis_exponents_extended[atom_end] =
                            this->results.Z[excited_position]
                            * this->results.wavefunction_lenght_multipliers[excited_position];
                            Gaussian_basis_extended[atom_end] = Gaussian_basis[excited_position];
                            }
                        }
        // set back set ground state of quantum numbers and exponent
        this->Create_excitation(excited_position, n_ground, l_ground, m_ground, s_ground, true);
        this->results.wavefunction_lenght_multipliers[excited_position] = wavefunction_lenght_multiplier_ground;
        // froze orbital for excitation
        this->results.wavefunction_constraints[excited_position] = 0;
        // restore integrals of ground state
        for (j = 0; j < count_electrons; j++)
            this->nuclear_atraction_integral_matrix[(excited_position * count_electrons) + j] =
            nuclear_atraction_integral_matrix_row[j];
        for (j = 0; j < count_electrons; j++)
            this->nucleuses_distances[(excited_position * count_electrons)  + j] = nucleuses_distances_row[j];
        for (j = 0; j < count_electrons; j++)
            this->nucleuses_atractions[(excited_position * count_electrons)  + j] = nucleuses_atractions_row[j];
        for (j = 0; j < count_electrons; j++)
            this->coulombic_integral_matrix[(excited_position * count_electrons)  + j] = coulombic_integral_matrix_row[j];
        for (j = 0; j < count_electrons; j++)
            this->overlap_integral_matrix[(excited_position * count_electrons)  + j] = overlap_integral_matrix_row[j];
        for (j = 0; j < count_electrons; j++)
            this->resonance_integral_matrix[(excited_position * count_electrons)  + j] = resonance_integral_matrix_row[j];
        for (j = 0; j < count_electrons; j++)
            this->kinetic_integral_matrix[(excited_position * count_electrons)  + j] = kinetic_integral_matrix_row[j];
        for (j = 0; j < count_electrons; j++)
            {
            this->nuclear_atraction_integral_matrix[(excited_position * j) + count_electrons] =
            this->nuclear_atraction_integral_matrix[(excited_position * count_electrons) + j];
            this->coulombic_integral_matrix[(excited_position * j) + count_electrons] =
            this->coulombic_integral_matrix[(excited_position * count_electrons)  + j];
            this->overlap_integral_matrix[(excited_position * j) + count_electrons] =
            this->overlap_integral_matrix[(excited_position * count_electrons)  + j];
            this->resonance_integral_matrix[(excited_position * j) + count_electrons] =
            this->resonance_integral_matrix[(excited_position * count_electrons)  + j];
            this->kinetic_integral_matrix[(excited_position * j) + count_electrons] =
            this->kinetic_integral_matrix[(excited_position * count_electrons)  + j];
            }
        }
    return(Groung_Hamiltoninan);
    }
template <typename T>
string Slater_basis_set_calculations_DFT<T>::Create_input_from_coordinates(vector<string> species, vector<string> x, vector<string> y, vector<string> z)
    {
    unsigned int i, j;
    unsigned int count_species;
    unsigned int count_layers;
    vector<unsigned int> counts_layers;
    string input = "";
    string ng_1 = "He_";
    string ng_2 = "Ne_";
    string ng_3 = "Ar_";
    string ng_4 = "Kr_";
    string ng_5 = "Xe_";
    string ng_6 = "Rn_";
    string specie;
    size_t found;
    
    string x_string, y_string, z_string;
    
    vector<string> elements = {"H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
    
    count_species = species.size();
    counts_layers.reserve(count_species);
    // Check input vectors
    if (count_species != x.size() or count_species != y.size() or count_species != z.size())
        {
        input = "wrong size of coordinate vector";
        return(input);
        }
    // Generate string of atoms
    for (i = 0; i < count_species; i++)
        {
        specie = species[i];
        for (j = 0; j < 118; j++)
            {
            found = specie.find(elements[j]);
            if (found != string::npos)
                {
                input.append(specie + "_");
                count_layers = 1;
                if (j > 2) {
                    input.append(ng_1);
                    count_layers++;
                    }
                if (j > 10) {
                    input.append(ng_2);
                    count_layers++;
                    }
                if (j > 18) {
                    input.append(ng_3);
                    count_layers++;
                    }
                if (j > 36) {
                    input.append(ng_4);
                    count_layers++;
                    }
                if (j > 54) {
                    input.append(ng_5);
                    count_layers++;
                    }
                if (j > 86) {
                    input.append(ng_6);
                    count_layers++;
                    }
                counts_layers.push_back(count_layers);
                break;
                }
            }
        }
    // Append part of coordinates
    for (i = 0; i < count_species; i++)
        {
        x_string = x[i];
        y_string = y[i];
        z_string = z[i];
        count_layers = counts_layers[i];
        for (j = 0; j < count_layers; j++)
            {
            if (i == 0 and j == 0)
                continue;
            
            input.append("[");
            input.append(x_string);
            input.append("_");
            input.append(y_string);
            input.append("_");
            input.append(z_string);
            input.append("]");
            }
        }
    return(input);
    }
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Set_spins_and_bonds(unsigned int size_order, vector<T>* values)
    {
    unsigned int i, j, k;
    unsigned int matrix_order;
    
    unsigned int count_unpaired = 0;
    vector<T> overlap_rows;
    T overlap_row;
    unsigned int highest_row;
    
    matrix_order = this->results.n.size();
    // Set all electrons bonding
    for (i = 0; i < matrix_order; i++)
        {
        count_unpaired++;
        if (this->results.bonding[i] == -1)
            {
            count_unpaired++;
            this->results.bonding[i] = i;
            }
        }
    // Compute overlap integral matrix
    if (allocation_memory == true)
        Execute_calculation(1, 1, size_order, false, values);
    
    if (this->Create_overlap_integral_matrix(this->overlap_integral_matrix, matrix_order, &this->results) == -1)
        return(-1);
    // Set spins with regards to overlap matrix to minimizing the overlaps
    overlap_rows.resize(matrix_order);
    for (i = 0; i < count_unpaired; i++)
        {
        for (j = 0; j < matrix_order; j++)
            {
            overlap_row = 0;
            for (k = 0; k < matrix_order; k++)
                {
                overlap_row += this->overlap_integral_matrix[j * matrix_order + k];
                }
            overlap_row -= 1;
            overlap_rows[j] = overlap_row;
            }
        highest_row = 0;
        for (j = 0; j < matrix_order; j++)
            if (overlap_rows[j] > overlap_rows[highest_row])
                highest_row = j;
        
        if (overlap_rows[highest_row] > 1)
            {
            // Update spins
            this->results.spins[highest_row] = - this->results.spins[highest_row];
            // Update overlap matrix
            for (j = 0; j < matrix_order; j++)
                {
                if (j =! highest_row)
                    this->overlap_integral_matrix[j * matrix_order + highest_row] =
                    -this->overlap_integral_matrix[j * matrix_order + highest_row];
                }
            for (j = 0; j < matrix_order; j++)
                {
                if (j =! highest_row)
                    this->overlap_integral_matrix[j + highest_row * matrix_order] =
                    -this->overlap_integral_matrix[j + highest_row * matrix_order];
                }
            }
        }
    return(0);
    }
template <typename T>
vector<T> Slater_basis_set_calculations_DFT<T>::Compute_correlation_energies(unsigned int density_matrix_order, T* density_matrix, T* Hamiltonian_matrix)
    {
    vector<T> basis_energy_levels;
    // Back computing the energies to basis
    // Must be performed in the begin and end of CASSCF method and resulting strings must be subtracted
    
    // Transpose density matrix
    
    // Multiply the Hamiltonian by the density matrix
    
    // Diagonalize the resulting matrix
    return(basis_energy_levels);
    }
// end of CI and basis sets creating section

template <typename T>
int Slater_basis_set_calculations_DFT<T>::Calculate_Huckel_Matrix(T* Huckel_matrix, unsigned int* Huckel_matrix_order,
vector<unsigned int> atom_numbers)
    {
    unsigned int i, j;
    unsigned int count_pi_electrons;
    array<unsigned int, this->max_electrons> index_pi_electrons;
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
            if (j != i and electron_to_atom_numbers[index_pi_electrons[i]] != electron_to_atom_numbers[index_pi_electrons[j]])
                Huckel_matrix[(i * count_pi_electrons) + j] = Huckel_matrix[(i * count_pi_electrons) + j]
                + abs(this->resonance_integral_matrix[(index_pi_electrons[i] * matrix_order) +  index_pi_electrons[j]]);
        } /*
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
            */
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
int Slater_basis_set_calculations_DFT<T>::Detect_symetry_information(symetry_axes *symetry_axes, symetry_planes *symetry_planes)
    {
    unsigned int i, j, k, l, m;
    
    array<T, this->max_atoms> x;
    array<T, this->max_atoms> y;
    array<T, this->max_atoms> z;
    array<T, this->max_atoms> Z;
    T x_center; // coordinates of mass center
    T y_center;
    T z_center;
    // axes
    array<T, this->max_atoms> u_x_0;
    array<T, this->max_atoms> u_y_0;
    array<T, this->max_atoms> u_z_0;
    array<T, this->max_atoms> u_x; // (x + y + z) * t + T = 0
    array<T, this->max_atoms> u_y;
    array<T, this->max_atoms> u_z;
    vector<T> u_x_2(this->max_atoms * this->max_atoms);
    vector<T> u_y_2(this->max_atoms * this->max_atoms);
    vector<T> u_z_2(this->max_atoms * this->max_atoms);
    unsigned int count_axes = 0;
    // planes
    array<T, this->max_atoms> a; // a * x + b * y + c * z + d = 0
    array<T, this->max_atoms> b;
    array<T, this->max_atoms> c;
    array<T, this->max_atoms> d;
    // mirroring points
    array<T, this->max_atoms> x_mirroring;
    array<T, this->max_atoms> y_mirroring;
    array<T, this->max_atoms> z_mirroring;
    array<T, this->max_atoms> t;
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
        { // computing mirroring points
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
int Slater_basis_set_calculations_DFT<T>::D3_vector_multiply(T* a, T* b, T* c)
    {c[0] = (a[1] * b[2]) - (a[2] * b[1]); c[1] = (a[2] * b[0]) - (a[0] * b[2]); c[2] = (a[0] * b[1]) - (a[1] * b[0]);
    return(0);}
template <typename T>
int Slater_basis_set_calculations_DFT<T>::Clear()
    {
    unsigned int i, j;
    bool previous_deleted;
    
    Slater_basis_set_calculations<T>::Clear();
    allocation_memory = true;
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
    // Density functional theory section
    allocation_PBE = true;
    
    for (i = 0; i < atoms_spin_densities.size(); i++)
        if (atoms_spin_densities[i] != nullptr)
            delete[] atoms_spin_densities[i];
    for (i = 0; i < atoms_electron_densities.size(); i++)
        if (atoms_electron_densities[i] != nullptr)
            delete[] atoms_electron_densities[i];
    for (i = 0; i < atoms_Fi_densities.size(); i++)
        if (atoms_Fi_densities[i] != nullptr)
            delete[] atoms_Fi_densities[i];
    for (i = 0; i < atoms_gradients_densities.size(); i++)
        if (atoms_gradients_densities[i] != nullptr)
            delete[] atoms_gradients_densities[i];
    
    atoms_electron_densities.clear();
    atoms_spin_densities.clear();
    atoms_Fi_densities.clear();
    atoms_gradients_densities.clear();
    pre_PBE_wavefunction_lenght_multipliers.clear();
    correlation_energies.clear();
    exchange_energies.clear();
    electron_energies.clear();
    kappa = 0.804; // Homogenous electron gass
    mi = (beta * this->Pi * this->Pi)/3.00;
    // end of density functional theory section
    // VQE section
    VQE_Wavefunction_lenght_multipliers_1.clear();
    VQE_Wavefunction_lenght_multipliers_2.clear();
    VQE_Eigenvectors_1.clear();
    VQE_Eigenvectors_2.clear();
    VQE_correlation_energy_sign.clear();
    VQE_ansatz.clear();
    VQE_gates.clear();
    circuit.data.clear();
    T VQE_Hamiltonian = 0;
    T VQE_previous_Hamiltonian = 0;
    VQE_iterations = 0;
    VQE_max_iterations = 10;
    VQE_shots = 1000;
    VQE_max_connections = 8;
    VQE_interconnection_distance = 5;
    Bravyi_Kitaev = false;
    VQE_gates.clear();
    // end of VQE section
    // TPSS section
    allocation_TPSS = true;
    
    for (i = 0; i < electrons_gradients_2_densities.size(); i++)
            {
            previous_deleted = false; // avoid double delete
            for (j = 0; j < i; j++)
                if (electrons_gradients_2_densities[i] == electrons_gradients_2_densities[j])
                    previous_deleted = true;
            
            if (previous_deleted == false) {
                delete[] electrons_gradients_2_densities[i];
                }
            }
    for (i = 0; i < atoms_sum_electrons_gradients_2_densities.size(); i++)
        if (atoms_sum_electrons_gradients_2_densities[i] != nullptr)
            delete[] atoms_sum_electrons_gradients_2_densities[i];
    for (i = 0; i < atoms_gradients_2_densities.size(); i++)
        if (atoms_gradients_2_densities[i] != nullptr)
            delete[] atoms_gradients_2_densities[i];
    for (i = 0; i < atoms_spins_gradients.size(); i++)
        if (atoms_spins_gradients[i] != nullptr)
            delete[] atoms_spins_gradients[i];
    
    electrons_gradients_2_densities.clear();
    atoms_sum_electrons_gradients_2_densities.clear();
    atoms_gradients_2_densities.clear();
    atoms_spins_gradients.clear();
    // End of TPSS section
    // Orthonormalizing and Gaussian export section
    if (positive_intraatomic_overlap_matrix != nullptr)
        delete[] positive_intraatomic_overlap_matrix;
    /*if (orthonormalizing_matrix != nullptr)
        delete[] orthonormalizing_matrix;
    if (orthonormalized_Hamiltonian_matrix != nullptr)
        delete[] orthonormalized_Hamiltonian_matrix;*/
    
    Gaussian_basis.clear();
    /* vector<T*> molecular_orbitals; Backup of pointers to wavefunction Allocated form a previous class instance */
    // End of orthonormalizing and Gaussian export section
    // CI and basis sets creating section
    count_basis_per_atom = 0;
    count_basis_overall = 0;
    n_extended.clear();
    l_extended.clear();
    m_extended.clear();
    s_extended.clear();
    counts_ground_basis.clear();
    
    Slater_basis_exponents_extended.clear();
    Gaussian_basis_extended.clear();
    correlation_energies_extended.clear();
    // End of CI and basis sets creating section
    return(0);
    }
template <typename T>
Slater_basis_set_calculations_DFT<T>::~Slater_basis_set_calculations_DFT(){
    Clear();}
template class Slater_basis_set_calculations_DFT<double>;

/*
Author of this source code Ing. Pavel Florian Ph.D. licensed this source code under the the Apache License:

Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/
*/ 
