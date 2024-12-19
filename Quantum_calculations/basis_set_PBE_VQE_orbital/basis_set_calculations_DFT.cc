#include "basis_set_calculations_DFT.h"
using namespace std;

template <typename T>
T basis_set_calculations_DFT<T>::Create_index_atoms()
    {
    unsigned int i;
    unsigned int sum_electrons = this->results.n.size();
    unsigned int count_bonds;
    T* x = this->results.x.data();
    T* y = this->results.y.data();
    T* z = this->results.z.data();

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
T basis_set_calculations_DFT<T>::Detect_center_of_mass()
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
T basis_set_calculations_DFT<T>::Detect_dipole_moments(dipole_moment *dipole_moment)
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
            next_atom_begin = count_electrons - 1;
        
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
T basis_set_calculations_DFT<T>::Create_crystal_field(central_cations *central_cations)
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
T basis_set_calculations_DFT<T>::Alocate_densities(unsigned int size_order)
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
            atoms_Ksi_densities.push_back(density);
            }
        }
    catch (int) {
        for (i = 0; i < atoms_Ksi_densities.size(); i++)
            delete [] atoms_Ksi_densities[i];
        
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
T basis_set_calculations_DFT<T>::Compute_density_thread(T** densities, T* atom_density, unsigned int begin, unsigned int end, unsigned int atom_size, vector<T> spins, vector<int> spin_paired, bool spin_density)
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
        atom_density[i] = 0;
    
    for (i = begin; i <= end; i++)
        {
        if (spin_density == false or spin_paired[i] == -1 or spin_paired[i] < begin
        or spin_paired[i] > end)
            {
            pointer_to_density = densities[i]; // Unpaired and bonding electrons are included for spin_paired[i] < begin
            if (spin_density == true and spins[i] == -0.5) // For spin density map and negative spins
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
        }
    // End of vectorisation code
    return(0);
    }
template <typename T>
T basis_set_calculations_DFT<T>::Compute_densities(vector<T*> densities, vector<T*> atoms_densities,
vector<unsigned int> index, unsigned int size_order, vector<T> spins, vector<int> spin_paired, bool spin_density)
    { // Compute densities of wavefunctions/probabilities for atoms in 1 thread for each atoms.
    unsigned int i, j;
    unsigned int begin, end;
    unsigned int count_atoms;
    unsigned int size_atoms;
    unsigned int begins_atoms[this->max_atoms];
    unsigned int ends_atoms[this->max_atoms];
    T* pointers_to_atoms[this->max_atoms];
    T* pointers_to_probabilities[this->max_atoms];
    thread t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
    bool t9_flag, t10_flag, t11_flag, t12_flag, t13_flag, t14_flag, t15_flag;

    t9_flag = false;
    t10_flag = false;
    t11_flag = false;
    t12_flag = false;
    t13_flag = false;
    t14_flag = false;
    t15_flag = false;
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
   for (i = 0; (i + 7) < count_atoms; i = i + 8)
           { 
           t1 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0], pointers_to_atoms[i],
           begins_atoms[i], ends_atoms[i], size_atoms, spins, spin_paired, spin_density);
           t2 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0], pointers_to_atoms[i + 1],
           begins_atoms[i + 1], ends_atoms[i + 1], size_atoms, spins, spin_paired, spin_density);
           t3 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0], pointers_to_atoms[i + 2],
           begins_atoms[i + 2], ends_atoms[i + 2], size_atoms, spins, spin_paired, spin_density);
           t4 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0], pointers_to_atoms[i + 3],
           begins_atoms[i + 3], ends_atoms[i + 3], size_atoms, spins, spin_paired, spin_density);
           t5 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0], pointers_to_atoms[i + 4],
           begins_atoms[i + 4], ends_atoms[i + 4], size_atoms, spins, spin_paired, spin_density);
           t6 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0], pointers_to_atoms[i + 5],
           begins_atoms[i + 5], ends_atoms[i + 5], size_atoms, spins, spin_paired, spin_density);
           t7 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0], pointers_to_atoms[i + 6],
           begins_atoms[i + 6], ends_atoms[i + 6], size_atoms, spins, spin_paired, spin_density);
           t8 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0], pointers_to_atoms[i + 7],
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
           t9 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0],
           pointers_to_atoms[count_atoms - 7], begins_atoms[count_atoms - 7], ends_atoms[count_atoms - 7],
           size_atoms, spins, spin_paired, spin_density);
           t9_flag = true;
           }
       if (count_atoms % 8 >= 6)
           {
           t10 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0],
           pointers_to_atoms[count_atoms - 6], begins_atoms[count_atoms - 6], ends_atoms[count_atoms - 6],
           size_atoms, spins, spin_paired, spin_density);
           t10_flag = true;
           }
       if (count_atoms % 8 >= 5)
           {
           t11 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0],
           pointers_to_atoms[count_atoms - 5], begins_atoms[count_atoms - 5], ends_atoms[count_atoms - 5],
           size_atoms, spins, spin_paired, spin_density);
           t11_flag = true;
           }
       if (count_atoms % 8 >= 4)
           {
           t12 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0],
           pointers_to_atoms[count_atoms - 4], begins_atoms[count_atoms - 4], ends_atoms[count_atoms - 4],
           size_atoms, spins, spin_paired, spin_density);
           t12_flag = true;
           }
       if (count_atoms % 8 >= 3)
           {
           t13 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0],
           pointers_to_atoms[count_atoms - 3], begins_atoms[count_atoms - 3], ends_atoms[count_atoms - 3],
           size_atoms, spins, spin_paired, spin_density);
           t13_flag = true;
           }
       if (count_atoms % 8 >= 2)
           {
           t14 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0],
           pointers_to_atoms[count_atoms - 2], begins_atoms[count_atoms - 2], ends_atoms[count_atoms - 2],
           size_atoms, spins, spin_paired, spin_density);
           t14_flag = true;
           }
       if (count_atoms % 8 >= 1)
           {
           t15 = thread(&basis_set_calculations_DFT::Compute_density_thread, this, &pointers_to_probabilities[0],
           pointers_to_atoms[count_atoms - 1], begins_atoms[count_atoms - 1], ends_atoms[count_atoms - 1],
           size_atoms, spins, spin_paired, spin_density);
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
T basis_set_calculations_DFT<T>::Correct_densities(vector<T*> atoms_densities_list, vector<T>  x, vector <T> y, vector <T> z,
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
T basis_set_calculations_DFT<T>::Compute_Ksi_thread(T* spin_density, T* electron_density, T* Ksi_density, unsigned int size_order)
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
        Ksi_density[i] = (pow(1 + pow_1, 2.0/3.0) + pow(1 - pow_1, 2.0/3.0))/2;
        Ksi_density[i + 1] = (pow(1 + pow_2, 2.0/3.0) + pow(1 - pow_2, 2.0/3.0))/2;
        Ksi_density[i + 2] = (pow(1 + pow_3, 2.0/3.0) + pow(1 - pow_3, 2.0/3.0))/2;
        Ksi_density[i + 3] = (pow(1 + pow_4, 2.0/3.0) + pow(1 - pow_4, 2.0/3.0))/2;
        Ksi_density[i + 4] = (pow(1 + pow_5, 2.0/3.0) + pow(1 - pow_5, 2.0/3.0))/2;
        Ksi_density[i + 5] = (pow(1 + pow_6, 2.0/3.0) + pow(1 - pow_6, 2.0/3.0))/2;
        Ksi_density[i + 6] = (pow(1 + pow_7, 2.0/3.0) + pow(1 - pow_7, 2.0/3.0))/2;
        Ksi_density[i + 7] = (pow(1 + pow_8, 2.0/3.0) + pow(1 - pow_8, 2.0/3.0))/2;
        }
    for (i = size - (size % 8); i < size; i++)
        {
        pow_1 = pow(spin_density[i], 2.0/3.0);
        Ksi_density[i] = (pow(1 + pow_1, 2.0/3.0) + pow(1 - pow_1, 2.0/3.0))/2;
        }
    // End of Vectorisation code
    return(0);
    }
template <typename T>
T basis_set_calculations_DFT<T>::Compute_gradient_thread(T* electron_density, T* gradient_density, unsigned int size_order)
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
T basis_set_calculations_DFT<T>::Compute_Ksi_and_gradients(vector<T*> atoms_electron_densities, vector<T*> atoms_spin_densities,
vector<T*> atoms_Ksi_densities, vector<T*> atoms_gradients_densities, unsigned int size_order)
    {
    unsigned int i;
    
    T* electron_densities[this->max_atoms];
    T* spin_densities[this->max_atoms];
    T* Ksi_densities[this->max_atoms];
    T* gradients_densities[this->max_atoms];
    
    
    thread t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45;
    thread t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
    bool t39_flag, t40_flag, t41_flag, t42_flag, t43_flag, t44_flag, t45_flag;
    bool t54_flag, t55_flag, t56_flag, t57_flag, t58_flag, t59_flag, t60_flag;
    
    t39_flag = false;
    t40_flag = false;
    t41_flag = false;
    t42_flag = false;
    t43_flag = false;
    t44_flag = false;
    t45_flag = false;
    
    t54_flag = false;
    t55_flag = false;
    t56_flag = false;
    t57_flag = false;
    t58_flag = false;
    t59_flag = false;
    t60_flag = false;
    
    for (i = 0; i < count_atoms; i++)
        {
        electron_densities[i] = atoms_electron_densities[i];
        spin_densities[i] = atoms_spin_densities[i];
        Ksi_densities[i] = atoms_Ksi_densities[i];
        gradients_densities[i] = atoms_gradients_densities[i];
        }
    
    for (i = 0; (i + 7) < count_atoms; i = i + 8)
        { // multithreading code
        t31 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[i], electron_densities[i],
        Ksi_densities[i], size_order);
        t32 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[i + 1], electron_densities[i + 1],
        Ksi_densities[i + 1], size_order);
        t33 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[i + 2], electron_densities[i + 2],
        Ksi_densities[i + 2], size_order);
        t34 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[i + 3], electron_densities[i + 3],
         Ksi_densities[i + 3], size_order);
        t35 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[i + 4], electron_densities[i + 4],
         Ksi_densities[i + 4], size_order);
        t36 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[i + 5], electron_densities[i + 5],
         Ksi_densities[i + 5], size_order);
        t37 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[i + 6], electron_densities[i + 6],
         Ksi_densities[i + 6], size_order);
        t38 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[i + 7], electron_densities[i + 7],
         Ksi_densities[i + 7], size_order);
        t31.join();
        t32.join();
        t33.join();
        t34.join();
        t35.join();
        t36.join();
        t37.join();
        t38.join();
        }
    if (count_atoms % 8 >= 7)
        {
        t39 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[count_atoms - 7],
        electron_densities[count_atoms - 7], Ksi_densities[count_atoms - 7], size_order);
        t39_flag = true;
        }
    if (count_atoms % 8 >= 6)
        {
        t40 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this,  spin_densities[count_atoms - 6],
        electron_densities[count_atoms - 6], Ksi_densities[count_atoms - 6], size_order);
        t40_flag = true;
        }
    if (count_atoms % 8 >= 5)
        {
        t41 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[count_atoms - 5],
        electron_densities[count_atoms - 5], Ksi_densities[count_atoms - 5], size_order);
        t41_flag = true;
        }
    if (count_atoms % 8 >= 4)
        {
        t42 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[count_atoms - 4],
        electron_densities[count_atoms - 4], Ksi_densities[count_atoms - 4], size_order);
        t42_flag = true;
        }
    if (count_atoms % 8 >= 3)
        {
        t43 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[count_atoms - 3],
        electron_densities[count_atoms - 3], Ksi_densities[count_atoms - 3], size_order);
        t43_flag = true;
        }
    if (count_atoms % 8 >= 2)
        {
        t44 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[count_atoms - 2],
        electron_densities[count_atoms - 2], Ksi_densities[count_atoms - 2], size_order);
        t44_flag = true;
        }
    if (count_atoms % 8 >= 1)
        {
        t45 = thread(&basis_set_calculations_DFT::Compute_Ksi_thread, this, spin_densities[count_atoms - 1],
        electron_densities[count_atoms - 1], Ksi_densities[count_atoms - 1], size_order);
        t45_flag = true;
        }
    if (t39_flag == true)
        {
        t39.join();
        t39_flag = false;
        }
    if (t40_flag == true)
        {
        t40.join();
        t40_flag = false;
        }
    if (t41_flag == true)
        {
        t41.join();
        t41_flag = false;
        }
    if (t42_flag == true)
        {
        t42.join();
        t42_flag = false;
        }
    if (t43_flag == true)
        {
        t43.join();
        t43_flag = false;
        }
    if (t44_flag == true)
        {
        t44.join();
        t44_flag = false;
        }
    if (t45_flag == true)
        {
        t45.join();
        t45_flag = false;
        }
    for (i = 0; (i + 7) < count_atoms; i = i + 8)
        { // multithreading code
        t46 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[i], gradients_densities[i],
        size_order);
        t47 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[i + 1], gradients_densities[i + 1],
        size_order);
        t48 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[i + 2], gradients_densities[i + 2],
        size_order);
        t49 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[i + 3], gradients_densities[i + 3],
        size_order);
        t50 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[i + 4], gradients_densities[i + 4],
        size_order);
        t51 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[i + 5], gradients_densities[i + 5],
        size_order);
        t52 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[i + 6], gradients_densities[i + 6],
        size_order);
        t53 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[i + 7], gradients_densities[i + 7],
        size_order);
        t46.join();
        t47.join();
        t48.join();
        t49.join();
        t50.join();
        t51.join();
        t52.join();
        t53.join();
        }
    if (count_atoms % 8 >= 7)
        {
        t54 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[count_atoms - 7],
        gradients_densities[count_atoms - 7], size_order);
        t54_flag = true;
        }
    if (count_atoms % 8 >= 6)
        {
        t55 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[count_atoms - 6],
        gradients_densities[count_atoms - 6], size_order);
        t55_flag = true;
        }
    if (count_atoms % 8 >= 5)
        {
        t56 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[count_atoms - 5],
        gradients_densities[count_atoms - 5], size_order);
        t56_flag = true;
        }
    if (count_atoms % 8 >= 4)
        {
        t57 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[count_atoms - 4],
        gradients_densities[count_atoms - 4], size_order);
        t57_flag = true;
        }
    if (count_atoms % 8 >= 3)
        {
        t58 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[count_atoms - 3],
        gradients_densities[count_atoms - 3], size_order);
        t58_flag = true;
        }
    if (count_atoms % 8 >= 2)
        {
        t59 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[count_atoms - 2],
        gradients_densities[count_atoms - 2], size_order);
        t59_flag = true;
        }
    if (count_atoms % 8 >= 1)
        {
        t60 = thread(&basis_set_calculations_DFT::Compute_gradient_thread, this, electron_densities[count_atoms - 1],
        gradients_densities[count_atoms - 1], size_order);
        t60_flag = true;
        }
    if (t54_flag == true)
        {
        t54.join();
        t54_flag = false;
        }
    if (t55_flag == true)
        {
        t55.join();
        t55_flag = false;
        }
    if (t56_flag == true)
        {
        t56.join();
        t56_flag = false;
        }
    if (t57_flag == true)
        {
        t57.join();
        t57_flag = false;
        }
    if (t58_flag == true)
        {
        t58.join();
        t58_flag = false;
        }
    if (t59_flag == true)
        {
        t59.join();
        t59_flag = false;
        }
    if (t60_flag == true)
        {
        t60.join();
        t60_flag = false;
        }
    return(0);
    }
template <typename T>
T basis_set_calculations_DFT<T>::PBE_thread(T* atoms_gradients_density, T* atoms_Ksi_density, T* electron_density,
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
            fi_2 = atoms_Ksi_density[i] * atoms_Ksi_density[i];
            power = pow(X, T(s_2/fi_2));
            nl_integral_factor = log(1 + 1/(power + (power * power)));
            if ((not (isnan(nl_integral_factor))) and (not (isinf(nl_integral_factor))))
                correlation_integral = correlation_integral +
                (gamma * electron_density[i] * fi_2 * atoms_Ksi_density[i] * nl_integral_factor);
        
            exchange_integral = exchange_integral + e_x_unif * (1 + kappa - kappa/(1 + (mi * s_2)/kappa)) * electron_density[i];
            }
        }
    exchange_energy[0] = exchange_integral;
    correlation_energy[0] = correlation_integral;
    return(0);
    }
template <typename T>
T basis_set_calculations_DFT<T>::PBE_compute(vector<T*> atoms_gradients_densities, vector<T*> atoms_Ksi_densities,
vector<T*> electron_densities, unsigned int size_order)
    {
    unsigned int i, j, k;
    unsigned int count_orbitals;
    unsigned int sum_electrons;
    unsigned int index[this->max_electrons]; // Index of electrons positions forcomputing nuclear atraction integrals
    T wavefunction_coefficients[this->max_electrons];
    int bonding[this->max_electrons];
    int spin_paired[this->max_electrons];
    T spins[this->max_electrons];
    T Z[this->max_electrons];
    T exchange_array[this->max_electrons];
    T correlation_array[this->max_electrons];
    T* atoms_gradients_array[this->max_atoms];
    T* atoms_Ksi_array[this->max_atoms];
    T* electron_array[this->max_electrons];
    thread t61, t62, t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75;
    bool restriction;
    bool t69_flag, t70_flag, t71_flag, t72_flag, t73_flag, t74_flag, t75_flag;
    
    t69_flag = false;
    t70_flag = false;
    t71_flag = false;
    t72_flag = false;
    t73_flag = false;
    t74_flag = false;
    t75_flag = false;
    restriction = true;
    sum_electrons = this->results.n.size();
    count_orbitals = 0;
    
    if (sum_electrons > this->max_electrons)
        return(-1);
        
    for (i = 0; i < atoms_gradients_densities.size(); i++)
        atoms_gradients_array[i] = atoms_gradients_densities[i];
    
    for (i = 0; i < atoms_Ksi_densities.size(); i++)
        atoms_Ksi_array[i] = atoms_Ksi_densities[i];
        
    for (i = 0; i < electron_densities.size(); i++)
        electron_array[i] = electron_densities[i];
    
    for (i = 0; i < sum_electrons; i++)
        {
        wavefunction_coefficients[i] = this->results.wavefunction_coefficients[i];
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
    for (i = 0; (i + 7) < count_orbitals; i = i + 8)
        {
        t61 = thread(&basis_set_calculations_DFT::PBE_thread, this, atoms_gradients_array[electron_to_atom_numbers[index[i]]],
        atoms_Ksi_array[electron_to_atom_numbers[index[i]]], electron_array[index[i]],
        size_order, &exchange_array[index[i]], &correlation_array[index[i]], Wigner_Seitz_radiuses[Z[index[i]] - 1]);
        t62 = thread(&basis_set_calculations_DFT::PBE_thread, this, atoms_gradients_array[electron_to_atom_numbers[index[i + 1]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[i + 1]] - 1], electron_array[index[i + 1]],
        size_order, &exchange_array[index[i + 1]], &correlation_array[index[i + 1]], Wigner_Seitz_radiuses[Z[index[i + 1]] - 1]);
        t63 = thread(&basis_set_calculations_DFT::PBE_thread, this, atoms_gradients_array[electron_to_atom_numbers[index[i + 2]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[i + 2]] - 1], electron_array[index[i + 2]],
        size_order, &exchange_array[index[i + 2]], &correlation_array[index[i + 2]], Wigner_Seitz_radiuses[Z[index[i + 2]] - 1]);
        t64 = thread(&basis_set_calculations_DFT::PBE_thread, this, atoms_gradients_array[electron_to_atom_numbers[index[i + 3]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[i + 3]] - 1], electron_array[index[i + 3]],
        size_order, &exchange_array[index[i + 3]], &correlation_array[index[i + 3]], Wigner_Seitz_radiuses[Z[index[i + 3]] - 1]);
        t65 = thread(&basis_set_calculations_DFT::PBE_thread, this, atoms_gradients_array[electron_to_atom_numbers[index[i + 4]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[i + 4]] - 1], electron_array[index[i + 4]],
        size_order, &exchange_array[index[i + 4]], &correlation_array[index[i + 4]], Wigner_Seitz_radiuses[Z[index[i + 4]] - 1]);
        t66 = thread(&basis_set_calculations_DFT::PBE_thread, this, atoms_gradients_array[electron_to_atom_numbers[index[i + 5]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[i + 5]] - 1], electron_array[index[i + 5]],
        size_order, &exchange_array[index[i + 5]], &correlation_array[index[i + 5]], Wigner_Seitz_radiuses[Z[index[i + 5]] - 1]);
        t67 = thread(&basis_set_calculations_DFT::PBE_thread, this, atoms_gradients_array[electron_to_atom_numbers[index[i + 6]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[i + 6]] - 1], electron_array[index[i + 6]],
        size_order, &exchange_array[index[i + 6]], &correlation_array[index[i + 6]], Wigner_Seitz_radiuses[Z[index[i + 6]] - 1]);
        t68 = thread(&basis_set_calculations_DFT::PBE_thread, this, atoms_gradients_array[electron_to_atom_numbers[index[i + 7]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[i + 7]] - 1], electron_array[index[i + 7]],
        size_order, &exchange_array[index[i + 7]], &correlation_array[index[i + 7]], Wigner_Seitz_radiuses[Z[index[i + 7]] - 1]);
        
        t61.join();
        t62.join();
        t63.join();
        t64.join();
        t65.join();
        t66.join();
        t67.join();
        t68.join();
        }
    if (count_orbitals % 8 >= 7)
        {
        t69 = thread(&basis_set_calculations_DFT::PBE_thread, this,
        atoms_gradients_array[electron_to_atom_numbers[index[count_orbitals - 7]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[count_orbitals - 7]] - 1],
        electron_array[index[count_orbitals - 7]],
        size_order, &exchange_array[index[count_orbitals - 7]], &correlation_array[index[count_orbitals - 7]],
        Wigner_Seitz_radiuses[Z[index[count_orbitals - 7]] - 1]);
        t69_flag = true;
        }
    if (count_orbitals % 8 >= 6)
        {
        t70 = thread(&basis_set_calculations_DFT::PBE_thread, this,
        atoms_gradients_array[electron_to_atom_numbers[index[count_orbitals - 6]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[count_orbitals - 6]] - 1],
        electron_array[index[count_orbitals - 6]],
        size_order, &exchange_array[index[count_orbitals - 6]], &correlation_array[index[count_orbitals - 6]],
        Wigner_Seitz_radiuses[Z[index[count_orbitals - 6]] - 1]);
        t70_flag = true;
        }
    if (count_orbitals % 8 >= 5)
        {
        t71 = thread(&basis_set_calculations_DFT::PBE_thread, this,
        atoms_gradients_array[electron_to_atom_numbers[index[count_orbitals - 5]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[count_orbitals - 5]] - 1],
        electron_array[index[count_orbitals - 5]],
        size_order, &exchange_array[index[count_orbitals - 5]], &correlation_array[index[count_orbitals - 5]],
        Wigner_Seitz_radiuses[Z[index[count_orbitals - 5]] - 1]);
        t71_flag = true;
        }
    if (count_orbitals % 8 >= 4)
        {
        t72 = thread(&basis_set_calculations_DFT::PBE_thread, this,
        atoms_gradients_array[electron_to_atom_numbers[index[count_orbitals - 4]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[count_orbitals - 4]] - 1],
        electron_array[index[count_orbitals - 4]],
        size_order, &exchange_array[index[count_orbitals - 4]], &correlation_array[index[count_orbitals - 4]],
        Wigner_Seitz_radiuses[Z[index[count_orbitals - 4]] - 1]);
        t72_flag = true;
        }
    if (count_orbitals % 8 >= 3)
        {
        t73 = thread(&basis_set_calculations_DFT::PBE_thread, this,
        atoms_gradients_array[electron_to_atom_numbers[index[count_orbitals - 3]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[count_orbitals - 3]] - 1],
        electron_array[index[count_orbitals - 3]],
        size_order, &exchange_array[index[count_orbitals - 3]], &correlation_array[index[count_orbitals - 3]],
        Wigner_Seitz_radiuses[Z[index[count_orbitals - 3]] - 1]);
        t73_flag = true;
        }
    if (count_orbitals % 8 >= 2)
        {
        t74 = thread(&basis_set_calculations_DFT::PBE_thread, this,
        atoms_gradients_array[electron_to_atom_numbers[index[count_orbitals - 2]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[count_orbitals - 2]] - 1],
        electron_array[index[count_orbitals - 2]],
        size_order, &exchange_array[index[count_orbitals - 2]], &correlation_array[index[count_orbitals - 2]],
        Wigner_Seitz_radiuses[Z[index[count_orbitals - 2]] - 1]);
        t74_flag = true;
        }
   if (count_orbitals % 8 >= 1)
        {
        t75 = thread(&basis_set_calculations_DFT::PBE_thread, this,
        atoms_gradients_array[electron_to_atom_numbers[index[count_orbitals - 1]] - 1],
        atoms_Ksi_array[electron_to_atom_numbers[index[count_orbitals - 1]] - 1],
        electron_array[index[count_orbitals - 1]],
        size_order, &exchange_array[index[count_orbitals - 1]], &correlation_array[index[count_orbitals - 1]],
        Wigner_Seitz_radiuses[Z[index[count_orbitals - 1]] - 1]);
        t75_flag = true;
        }
    if (t69_flag == true)
        {
        t69.join();
        t69_flag = false;
        }
    if (t70_flag == true)
        {
        t70.join();
        t70_flag = false;
        }
    if (t71_flag == true)
        {
        t71.join();
        t71_flag = false;
        }
    if (t72_flag == true)
        {
        t72.join();
        t72_flag = false;
        }
    if (t73_flag == true)
        {
        t73.join();
        t73_flag = false;
        }
    if (t74_flag == true)
        {
        t74.join();
        t74_flag = false;
        }
    if (t75_flag == true)
        {
        t75.join();
        t75_flag = false;
        }
    // end of multithreading code
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
T basis_set_calculations_DFT<T>::Set_circuit()
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
T basis_set_calculations_DFT<T>::Run_circuit(vector<double>& ansatz)
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
template <typename T>
T basis_set_calculations_DFT<T>::String_to_advanced_parameters(string UI_input, unsigned int size_order,
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
T basis_set_calculations_DFT<T>::Execute_calculation(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order,
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
       
    Hamiltonian = this->Calculate(max_iterations, minimal_fidelity, size_order, allocation_memory,
    values, spin_density_vector, spin_values);
    if (dealocate == true) // Dealocate atom electron and spin densities
        Clear();
    return(Hamiltonian);
    }
// Density functional theory section
template <typename T>
T basis_set_calculations_DFT<T>::Execute_PBE(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool dealocate, vector<T>* values, vector<T>* spin_density_vector, vector<T>* spin_values)
    {
    unsigned int i, j, k;
    unsigned int matrix_order;
    T exchange_energy;
    T Hamiltonian;
    T electron_energy;
    vector<T> correction_diagonal;
    vector <T> basis_exchange_energies;
    thread t61, t62;
    
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
    
    if (allocation_memory == true)
        {
        Hamiltonian = this->Calculate(max_iterations, minimal_fidelity, size_order, allocation_memory,
        values, spin_density_vector, spin_values);
        
        for (i = 0; i < matrix_order; i++)
            pre_PBE_wavefunction_lenght_multipliers.push_back(this->results.wavefunction_lenght_multipliers[i]);
        }
    // Execute calculations
    for (i = 0; i < matrix_order * matrix_order; i++) // avoiding executing PBE for atomic systems
        {
        if (this->resonance_integral_matrix[i] != 0)
            break;
        
        if (i == matrix_order * matrix_order -1)
            return(Hamiltonian);
        }
    for (i = 0; i < matrix_order; i++) //
        {
        exchange_energy = 0;
        for (j = 0; j < matrix_order; j++)
            exchange_energy = exchange_energy + this->resonance_integral_matrix[(i * matrix_order) + j];
        
        basis_exchange_energies.push_back(exchange_energy); // Copy HF exchange energies to vector
        }
    if (allocation_memory == true)
        if (Alocate_densities(size_order) == -1)
            return(-1);
    // compute electron and spin densities
    t61 = thread(&basis_set_calculations_DFT::Compute_densities, this, this->results.probabilities, atoms_electron_densities,
    index_atoms, size_order, this->results.spins, this->results.spin_paired, false);
    t62 = thread(&basis_set_calculations_DFT::Compute_densities, this, this->results.probabilities, atoms_spin_densities,
    index_atoms, size_order, this->results.spins, this->results.spin_paired, true);
    t61.join();
    t62.join();
    // correct electron and spin densities
    Correct_densities(atoms_electron_densities, this->results.x, this->results.y, this->results.z, size_order);
    Correct_densities(atoms_spin_densities, this->results.x, this->results.y, this->results.z, size_order);
    correlation_energies.clear();
    exchange_energies.clear();
    electron_energies.clear();
    // run PBE method over electron and spin densities
    Compute_Ksi_and_gradients(atoms_electron_densities, atoms_spin_densities, atoms_Ksi_densities, atoms_gradients_densities,
    size_order);
    PBE_compute(atoms_gradients_densities, atoms_Ksi_densities, this->results.probabilities, size_order);
    
    for (i = 0; i < PBE_iterations; i++) 
        {
        // repeat calculation with computed electron correlation energies
        for (j = 0; j < matrix_order; j++)
            this->correction_matrix[j * (matrix_order + 1)] = correlation_energies[j] + exchange_energies[j]
            + correction_diagonal[j]; // PBE density functional not include polar part of bonds
        // Execute HF calculations
        Hamiltonian = this->Calculate(max_iterations, minimal_fidelity, size_order, false, values, spin_density_vector,  spin_values);
        matrix_order = this->results.n.size();
    
        for (i = 0; i < matrix_order; i++) //
            {
            exchange_energy = 0;
            for (j = 0; j < matrix_order; j++)
                exchange_energy = exchange_energy + this->resonance_integral_matrix[(i * matrix_order) + j];
        
            basis_exchange_energies.push_back(exchange_energy); // Copy HF exchange energies to vector
            }
        // compute electron and spin densities
        t61 = thread(&basis_set_calculations_DFT::Compute_densities, this, this->results.probabilities, atoms_electron_densities,
        index_atoms, size_order, this->results.spins, this->results.spin_paired, false);
        t62 = thread(&basis_set_calculations_DFT::Compute_densities, this, this->results.probabilities, atoms_spin_densities,
        index_atoms, size_order, this->results.spins, this->results.spin_paired, true);
        t61.join();
        t62.join();
        // correct electron and spin densities
        Correct_densities(atoms_electron_densities, this->results.x, this->results.y, this->results.z, size_order);
        Correct_densities(atoms_spin_densities, this->results.x, this->results.y, this->results.z, size_order);
        correlation_energies.clear();
        exchange_energies.clear();
        electron_energies.clear();
        // run PBE method over electron and spin densities
        Compute_Ksi_and_gradients(atoms_electron_densities, atoms_spin_densities, atoms_Ksi_densities, atoms_gradients_densities,
        size_order);
        PBE_compute(atoms_gradients_densities, atoms_Ksi_densities, this->results.probabilities, size_order);
        }
    if (dealocate == true) // Dealocate atom electron aand spin densities
        Clear();
        
    return(Hamiltonian);
    }
// end of density functional theory section
// VQE section
template <typename T>
T basis_set_calculations_DFT<T>::Execute_PBE_VQE(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool dealocate,
vector<T>* values, vector<T>* spin_density_vector, vector<T>* spin_values)
    {
    unsigned int i, j;
    unsigned int matrix_order;
    bool VQE;
    T VQE_ansatz_sum;
    vector<T> previous_spin_density_vector;
    vector<T> previous_spin_values;
    vector<T> previous_correlation_energies;
    vector<T> previous_exchange_energies;
    T max_PBE_difference = 0;
    
    VQE = true;
    allocation_memory = true;
    kappa = 0.804; // homogenic electron gass
    mi = (beta * this->Pi * this->Pi)/3.00;
    VQE_Hamiltonian = Execute_PBE(max_iterations, minimal_fidelity, size_order, false, values, spin_density_vector, spin_values);
    matrix_order = this->results.wavefunction_lenght_multipliers.size();
    
    for (i = 0; i < matrix_order * matrix_order; i++) // avoiding executing PBE and VQE for atomic systems
        {
        if (this->resonance_integral_matrix[i] != 0)
            break;
        
        if (i == matrix_order * matrix_order -1)
            return(VQE_Hamiltonian);
        }
    for (i = 0; i < matrix_order; i++) // Filling first wavefunction_lenghts_multipliers and eigenvectors vectors
        {
        VQE_Wavefunction_lenght_multipliers_1.push_back(this->results.wavefunction_lenght_multipliers[i]);
        VQE_Eigenvectors_1.push_back(values->operator[](i));
        previous_spin_density_vector.push_back(spin_density_vector->operator[](i));
        previous_correlation_energies.push_back(correlation_energies[i]);
        previous_exchange_energies.push_back(exchange_energies[i]);
        }
    for (i = 0; i < spin_values->size(); i++)
        previous_spin_values.push_back(spin_values->operator[](i));
    
    allocation_memory = false;
    for (i = 0; i < matrix_order; i++)
        this->results.wavefunction_lenght_multipliers[i] = pre_PBE_wavefunction_lenght_multipliers[i];
    
    kappa = 0.967; // Lieb-Oxford bonding
    mi = 0.235;
    VQE_previous_Hamiltonian = VQE_Hamiltonian;
    VQE_Hamiltonian = Execute_PBE(max_iterations, minimal_fidelity, size_order, false, values, spin_density_vector, spin_values);
    
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
            VQE_Hamiltonian = Execute_PBE(max_iterations, minimal_fidelity, size_order, false, values, spin_density_vector, spin_values);
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
        if (dealocate == true) // Dealocate atom electron aand spin densities
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
                values->operator[](i) = VQE_Eigenvectors_1[i];
                spin_density_vector->operator[](i) = previous_spin_density_vector[i];
                spin_values->operator[](i) = previous_spin_values[i];
                correlation_energies[i] = previous_correlation_energies[i];
                exchange_energies[i] = previous_exchange_energies[i];
                }
            for (i = 0; i < spin_values->size(); i++)
                spin_values->operator[](i) = previous_spin_values[i];
            }
        if (dealocate == true) // Dealocate atom electron aand spin densities
            Clear();
        if (VQE_correlation_energy_sign[0] == false)
            return VQE_Hamiltonian;
        else
            return VQE_previous_Hamiltonian;
        }
    }
// end of VQE section
// Gaussian export section
template <typename T>
T basis_set_calculations_DFT<T>::Gaussian_quadrature()
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
        {
        nodes = n[i] - l[i];
        exponents_multiplier = n[i] * wavefunction_lenght_multipliers[i];
        for (j = 0; j < 6; j++)
            {
            exponent = Aro_90_s1_exponents[j] * exponents_multiplier;
            coefficient = Aro_90_s1_coefficients[j];
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
                    coefficient = 1;
                    if (j % 2 == 0)
                        coefficient *= -1;
                    Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                    break;
                    }
                case 1: {
                    exponent = Aro_90_p_exponents[j] * exponents_multiplier;
                    coefficient = 1;
                    if (j % 2 == 0)
                        coefficient *= -1;
                    Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                    break;
                    }
                case 2: {
                    exponent = Aro_90_d_exponents[j] * exponents_multiplier;
                    coefficient = 1;
                    if (j % 2 == 0)
                        coefficient *= -1;
                    Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                    break;
                    }
                case 3: {
                    exponent = Aro_90_f_exponents[j] * exponents_multiplier;
                    coefficient = 1;
                    if (j % 2 == 0)
                        coefficient *= -1;
                    Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                    break;
                    }
                case 4: {
                    exponent = Aro_90_g_exponents[j] * exponents_multiplier;
                    coefficient = 1;
                    if (j % 2 == 0)
                        coefficient *= -1;
                    Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                    break;
                    }
                case 5: {
                    exponent = Aro_90_g_exponents[j] * exponents_multiplier;
                    coefficient = 1;
                    if (j % 2 == 0)
                        coefficient *= -1;
                    Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                    break;
                    }
                case 6: {
                    exponent = Aro_90_h_exponents[j] * exponents_multiplier;
                    coefficient = 1;
                    if (j % 2 == 0)
                        coefficient *= -1;
                    Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                    break;
                    }
                case 7: {
                    exponent = Aro_90_i_exponents[j] * exponents_multiplier;
                    coefficient = 1;
                    if (j % 2 == 0)
                        coefficient *= -1;
                    Gaussian_basis[i].push_back(pair<T, T>(exponent, coefficient));
                    break;
                    }
                case 8: {
                    exponent = Aro_90_k_exponents[j] * exponents_multiplier;
                    coefficient = 1;
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
// End of Gaussian export section
template <typename T>
T basis_set_calculations_DFT<T>::Calculate_Huckel_Matrix(T* Huckel_matrix, unsigned int* Huckel_matrix_order,
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
T basis_set_calculations_DFT<T>::Detect_symetry_information(symetry_axes *symetry_axes, symetry_planes *symetry_planes)
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
T basis_set_calculations_DFT<T>::D3_vector_multiply(T* a, T* b, T* c)
    {c[0] = (a[1] * b[2]) - (a[2] * b[1]); c[1] = (a[2] * b[0]) - (a[0] * b[2]); c[2] = (a[0] * b[1]) - (a[1] * b[0]);
    return(0);}
template <typename T>
T basis_set_calculations_DFT<T>::Clear()
    {
    unsigned int i;
    
    basis_set_calculations<T>::Clear();
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
    for (i = 0; i < atoms_spin_densities.size(); i++)
        if (atoms_spin_densities[i] != nullptr)
            delete[] atoms_spin_densities[i];
    for (i = 0; i < atoms_electron_densities.size(); i++)
        if (atoms_electron_densities[i] != nullptr)
            delete[] atoms_electron_densities[i];
    for (i = 0; i < atoms_Ksi_densities.size(); i++)
        if (atoms_Ksi_densities[i] != nullptr)
            delete[] atoms_Ksi_densities[i];
    for (i = 0; i < atoms_gradients_densities.size(); i++)
        if (atoms_gradients_densities[i] != nullptr)
            delete[] atoms_gradients_densities[i];
    
    atoms_electron_densities.clear();
    atoms_spin_densities.clear();
    atoms_Ksi_densities.clear();
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
    return(0);
    }
template <typename T>
basis_set_calculations_DFT<T>::~basis_set_calculations_DFT(){
    Clear();}
template class basis_set_calculations_DFT<double>;



/*
Author of this source code Ing. Pavel Florian Ph.D. licensed this source code under the the Apache License:

Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/  */ 
