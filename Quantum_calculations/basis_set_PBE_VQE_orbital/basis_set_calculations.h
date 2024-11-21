# ifndef BASIS_SET_CALCULATIONS_H
# define BASIS_SET_CALCULATIONS_H
#include <iostream>
#include <string>
#include <cmath>
#include <thread>
#include <vector>
#include <bits/stdc++.h>

using namespace std; // compiler parameters: -pthread -ffast-math -fno-finite-math-only -fmove-loop-invariants -O1

template <typename T>
class basis_set_calculations
{
public:
    vector<T> determinants;
    vector<T> spectra_EPR;
    vector<T> electron_spectra;
    
    const T Phi = (1.00 + sqrt(5))/2.00; // The golden ratio
    const T Pi = 4.00/sqrt(Phi); // Pi constant derived from Phi
    
    T E0 = 8.8541878128E-12; // F·m^-1, m^-3·kg^-1·s^4·A^2
    T h = 6.62607015E-34; // J·s , kg·m^2/s
    T vacuum_permeability = 1.25663706212E-6; // N·A^−2, kg·m·s^-2·A^−2
    T e = 1.602176634E-19; // C
    T me = 9.1093837E-31;
    T mp = 1.67262192E-27;
    T c = 299796000; // m·s^-1
    T Hartree_lenght = (E0 * h * h)/(Pi * me * 1 * e * e);
    T Hartree_energy_constant = me * (e * e /(2 * E0 * h)) * (e * e /(2 * E0 * h));
    T hyperfine_structure_constant = (e * e)/(2 * E0 * h * c);
    
    const unsigned int max_electrons = 1024; // constant for maximum electrons in basis set matrix
    const unsigned int max_atoms = 128; // constant for maximum atoms in basis set matrix
    unsigned int electron_number = 0;
    unsigned int iterations = 0;
    unsigned int determinant_exception_handle = 0;
    unsigned int overlap_integral_exception_handle = 0;
    
    bool Helium_correlation_energy = true;
    bool bonded_system = false;
    bool partial_overlap_matrix = true;
    
    T vector_lenght = 5.00; // constant for dimension of wavefunction vectors in Hartree lenghts
    T nucleus_repulsive_energy;
    T relative_permitivity = 1;
    T* nuclear_atraction_integral_matrix = nullptr;
    T* kinetic_integral_matrix = nullptr;
    T* nucleuses_distances = nullptr;
    T* nucleuses_atractions = nullptr;
    T* coulombic_integral_matrix = nullptr;
    T* overlap_integral_matrix = nullptr;
    T* overlap_efective_lenght_integral_matrix = nullptr;
    T* resonance_integral_matrix = nullptr;
    T* basis_set_matrix = nullptr;
    T* correction_matrix = nullptr;
    T* corr_basis_set_matrix = nullptr;
    T* spin_density_matrix = nullptr;
    
    vector <unsigned int> index_atoms;
    vector <unsigned int> electron_to_atom_numbers;
    vector <T> additive_atom_overlaps;
    
    
    struct atom_orbitals {
    vector<unsigned int> n;
    vector<unsigned int> l;
    vector<int> m;
    vector<T> s;
    vector<T> wavefunction_lenght_multipliers;
    vector<T> wavefunction_coefficients;
    vector <int> bonding;
    vector <bool> paired_with_previous;
    int charge;
    unsigned int reduced_Z;
    unsigned int Z;
    T electronegativity;
    } atoms;
    
    struct atom_wavefunctions {
    vector<T*> lenghts;
    vector<T*> wavefunctions;
    vector<T*> probabilities;
    vector<T*> Laplacians;
    vector<unsigned int> lenght_orders;
    vector<T> wavefunction_coefficients;
    vector<T> wavefunction_lenght_multipliers; // coefficient for conversion bohr radius to pixels
    vector<T> efective_radius_base; // elecron wavefunction radius fase for inicialization of calculation 
    vector<T> spins;
    vector<int> spin_paired;
    vector<int> bonding;
    vector <int> antibonding;
    vector <int> pi_bonding;
    vector <unsigned int> wavefunction_constraints;
    vector<unsigned int> n;
    vector<unsigned int> l;
    vector<int> m;
    vector<int> charge;
    vector<unsigned int> count_electrons;
    vector<T> reduced_Z;
    vector<unsigned int> Z;
    vector<unsigned int> electron_numbers;
    vector<T> x;
    vector<T> y;
    vector<T> z;
    } results;
    
    struct small_atom_wavefunctions {
    // This structure is for generating small wavefunction vectors for integration lenght of exchange integration
    vector<unsigned int> electron_numbers;
    vector<unsigned int> lenght_orders;
    vector<T*> lenghts;
    vector<T*> relative_lenghts;
    vector<T*> wavefunctions;
    vector<T*> probabilities;
    vector<T> x;
    vector<T> y;
    vector<T> z;
    vector<unsigned int> n;
    vector<unsigned int> l;
    vector<int> m;
    vector<unsigned int> Z;
    } small_results;
public:
// Section 1 - solving basis set matrixes - inherited from Huckel_calculations
inline T Determinant(unsigned int order, T* pointer, T* buffer, T* denominator, T* temp1, T* temp2);
T basis_set_Determinant_set(unsigned int order, T* pointer,unsigned int count, T min, T step, T* output_values);
T basis_set_Determinant_solver(unsigned int order, T* pointer);
private:
// Section 2 - generating the wavefunctions
T Wavefunction_lenghts_generate(T* lenghts, unsigned int lenght_order);

T Wavefunction_1s_generate(T* wavefunction, T* lenghts, int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_2s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_3s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_4s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
    
T Wavefunction_2px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_2pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_2py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_3px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_3pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_3py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_4px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_4pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_4py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_5px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_6px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_7px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
    
T Wavefunction_3dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_3dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_3dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_3dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_3dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_4dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_4dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_4dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_4dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_4dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_5dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_6dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_7dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
    
T Wavefunction_4fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
 T Wavefunction_4fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_4fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_4fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_4fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_4fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
 T Wavefunction_4fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
 
T Wavefunction_5fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_6fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_7fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
    
T Wavefunction_normalize(T* wavefunction_pointer, T normalisation_constant, unsigned int size);
T Orbitals_to_wavefunctions(unsigned int n, unsigned int l, int m, unsigned int lenght_order, T* wavefunction, T* lenghts,
unsigned int Z, T multiplier);
T Get_relative_Hartree_length(unsigned int Z, unsigned int n);
// Section 3 - mathematical operations for wavefunctions, probabilities densities and integrals
unsigned int Wavefunction_multiply(T* wavefunction_1, T* wavefunction_2, T* probabilities,
unsigned int lenght_order, T d_x, T d_y, T d_z);
T Wavefunction_multiply(T* wavefunction_1, T* wavefunction_2, T* probabilities, unsigned int lenght_order);
T Wavefunction_relative_lenghts_generate(T* lenghts, unsigned int lenght_order);
T Wavefunction_square(T* wavefunction_1, T* probabilities, unsigned int lenght_order);
T Probabilities_lenght(T* Probabilities, unsigned int lenght_order, int x, int y, int z);
T Probabilities_thread(T* Probabilities, unsigned int lenght_order, T* lenght);
T Probabilities_lenght(T* probabilities, unsigned int lenght_order);
T Laplacian_thread(T* Laplacian_1, T* wavefunction_2, unsigned int lenght_order);
T Integral_overlap(T* Wavefunction_1, T* Wavefunction_2, T* result, unsigned int lenght_order, T x, T y, T z);
T Integrate_Integral_overlap(T* wavefunction_1, T* wavefunction_2, T* result, unsigned int lenght_order, T x, T y, T z);
T Integral_coulombic(T radius_1, T radius_2, T distance, T* result, bool spin_bonded);
T Integrate_Integral_coulombic(T* density_1, T* density_2, T* result, unsigned int lenght_order, T x, T y, T z);
T Integral_nucleus_atraction(T probabilities_lenght, T multiplier, T* result, T* lenght, unsigned int Z);
T Integrate_Integral_nucleus_atraction(T* probabilities,
T* result, T* lenght, unsigned int lenght_order, T lenght_x, T lenght_y, T lenght_z, unsigned int Z);
T Integral_kinetic(T* Laplacian_1, T* Laplacian_2, T* result, unsigned int lenght_order, T d_x, T d_y, T d_z);
T Rydberg_energy(T Z, unsigned int n);
T Orbital_moment_energy(int m, T* B0);
T Spin_moment_energy(T s, T B0);
T Orbital_magnetic_field(T potential_energy, T radius, int l);
// Section 4 - generating lists of electrons
T Quantum_numbers_to_orbitals(unsigned int n, unsigned int l, int fulness, atom_orbitals* atom_orbitals_PTR);
T Atoms_to_valence_orbitals(string atom, atom_orbitals* atom_orbitals_PTR);
T Create_atomic_wavefunctions(atom_orbitals *atom_orbitals_PTR, atom_wavefunctions *atom_wavefunctions,
unsigned int size_order,T x,T y,T z);
T Create_bond_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions_1, atom_wavefunctions *atom_wavefunctions_2,
unsigned int count, T electronegativity_1, T electronegativity_2, T x_diference, T y_diference, T z_diference);
T Sum_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions_1, atom_wavefunctions *atom_wavefunctions_2);
// section 5 - generating matrices of integrals and Fock matrices
T Create_nuclear_atraction_integral_matrix(T* matrix, T* nucleuses, unsigned int order, atom_wavefunctions *atom_wavefunctions);
T Create_coulombic_integral_matrix(T* matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions,
small_atom_wavefunctions *small_atom_wavefunctions);
T Create_overlap_integral_matrix(T* matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions);
T Calculate_resonance_integral_matrix(T* overlap_matrix, T* overlap_efective_lenght_integral_matrix, T* resonance_integral_matrix,
unsigned int order, atom_wavefunctions *atom_wavefunctions, small_atom_wavefunctions *small_atom_wavefunctions);
T Calculate_kinetic_integral_matrix(T* matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions);
T Calculate_basis_set_matrix(T* nuclear_atraction_integral_matrix, T* coulombic_integral_matrix, T* resonance_integral_matrix,
T* kinetic_integral_matrix, T* basis_set_matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions);
T Calculate_corr_basis_set_matrix(T* basis_set_matrix, T* correction_matrix, T* corr_basis_set_matrix, unsigned int order);
T Calculate_spin_density_matrix(T* overlap_integral_matrix, T* spin_density_matrix, unsigned int order,
atom_wavefunctions *atom_wavefunctions);
T Solve_basis_set_matrix(T* basis_set_matrix, T* overlap_integral_matrix, unsigned int order,
vector<T>* values, atom_wavefunctions *atom_wavefunctions);
T Solve_spin_density_matrix(T* spin_density_matrix, unsigned int order, vector<T>* energetic_levels);
public:
// section 6 - completing basis set method and user interface handling
T Atom_orbitals_generate(string UI_input, atom_orbitals *atom_orbitals_PTR);
T Nucleus_repulsive_energy(atom_wavefunctions *atom_wavefunctions);
T Generate_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions, small_atom_wavefunctions *small_atom_wavefunctions,
unsigned int size_order, bool alocate, bool compute_densities);
T String_to_list_electrons(string UI_input, unsigned int size_order,
bool extern_coordinates, vector<T>* x_2, vector<T>* y_2, vector<T>* z_2);
T Calculate(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool alocate,
vector<T>* values, vector<T>* spin_density_vector,  vector<T>* spin_values);
T Clear();
~basis_set_calculations();};

#include "basis_set_calculations.cc"
# endif

/* BASIS_SET_CALCULATIONS_H

Author of this source code Ing. Pavel Florian Ph.D. licensed this source code under the the Apache License:
Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/
*/
