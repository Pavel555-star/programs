# ifndef Slater_basis_set_calculations_H
# define Slater_basis_set_calculations_H
#include <iostream>
#include <string>
#include <cmath>
#include <omp.h>
#include <array> 
#include <vector>
#include <bits/stdc++.h>
using namespace std; // compiler parameters: -fopenmp -O1 -ffast-math -fno-finite-math-only -fmove-loop-invariants 

template <typename T>
class Slater_basis_set_calculations
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
    
    static const unsigned int max_electrons = 1024; // constant for maximum electrons in basis set matrix
    static const unsigned int max_atoms = 512; // constant for maximum atoms in basis set matrix
    unsigned int electron_number = 0;
    unsigned int iterations = 0;
    unsigned int determinant_exception_handle = 0;
    unsigned int overlap_integral_exception_handle = 0;
    T small_wavefunction_size_factor = 1.5;
    T wavefunction_integration_threshold = 1.0E-15;
    
    bool Helium_correlation_energy = true;
    bool bonded_system = false;
    bool s1_memory_optimization = true; // if is true, memory for small_results grid is not allocated for s1 system
    
    T vector_lenght = 5.00; // constant for dimension of wavefunction vectors in Hartree lenghts
    T nucleus_repulsive_energy;
    T relative_permitivity = 1;
    T* nuclear_atraction_integral_matrix = nullptr;
    T* kinetic_integral_matrix = nullptr;
    T* nucleuses_distances = nullptr;
    T* nucleuses_atractions = nullptr;
    T* coulombic_integral_matrix = nullptr;
    T* overlap_integral_matrix = nullptr;
    T* overlap_effective_lenght_integral_matrix = nullptr;
    T* resonance_integral_matrix = nullptr;
    T* basis_set_matrix = nullptr;
    T* correction_matrix = nullptr;
    T* corr_basis_set_matrix = nullptr;
    T* spin_density_matrix = nullptr;
    
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
    vector<T*> Gradients;
    vector<unsigned int> lenght_orders;
    vector<unsigned int> x_range;
    vector<unsigned int> y_range;
    vector<unsigned int> z_range;
    vector<T> wavefunction_coefficients;
    vector<T> wavefunction_lenght_multipliers; // coefficient for conversion bohr radius to pixels
    vector<T> effective_radius_base; // elecron wavefunction radius fase for inicialization of calculation 
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
    vector<unsigned int> reduced_Z;
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
    vector<unsigned int> x_range;
    vector<unsigned int> y_range;
    vector<unsigned int> z_range;
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
T Determinant(unsigned int order, T* pointer, T* buffer, T* denominator, T* temp1, T* temp2);
int basis_set_Determinant_set(unsigned int order, T* pointer,unsigned int count, T min, T step, T* output_values);
int basis_set_Determinant_solver(unsigned int order, T* pointer);
protected:
// Section 2 - generating the wavefunctions
int Wavefunction_lenghts_generate(T* lenghts, unsigned int lenght_order);

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

T Wavefunction_5gz4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5gz3y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5gz3x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5gz2xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5gz2_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5gz_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5gz_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5gxy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_5gx4_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_6gz4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6gz3y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6gz3x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6gz2xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6gz2_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6gz_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6gz_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6gxy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6gx4_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_7gz4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7gz3y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7gz3x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7gz2xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7gz2_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7gz_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7gz_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7gxy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7gx4_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_6hz5_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6hz4y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6hz4x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6hz3xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6hz3_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6hz2_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6hz2_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6hz_x3y_xy3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6hz_x4_x2y2_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6hy_x4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_6hx_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_7hz5_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7hz4y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7hz4x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7hz3xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7hz3_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7hz2_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7hz2_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7hz_x3y_xy3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7hz_x4_x2y2_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7hy_x4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7hx_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

T Wavefunction_7iz6_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7iz5y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7iz5x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7iz4xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7iz4_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7iy__x2_y2__z3_zr2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7ix__x2_y2__z3_zr2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7i_xy__x2_y2__z2_r2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7i_x4_x2y2_y4__z2_r2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7i_yz__x4_x2y2_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7i_xz__x4_x2y2_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7i_yx5_x3y3_xy5_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);
T Wavefunction_7i_x6__x4y2_x2y4_y6_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order);

int Wavefunction_normalize(T* wavefunction_pointer, T normalisation_constant, unsigned int size);
int Orbitals_to_wavefunctions(unsigned int n, unsigned int l, int m, unsigned int lenght_order, T* wavefunction,
T* lenghts, unsigned int Z, T multiplier, unsigned int* x_range, unsigned int* y_range, unsigned int* z_range);
T Get_relative_Hartree_length(unsigned int Z, unsigned int n);
int Wavefunction_range_detect(T* wavefunction_pointer, unsigned int lenght_order, unsigned int* x_range,
unsigned int* y_range, unsigned int* z_range);
// Section 3 - mathematical operations for wavefunctions, probabilities densities and integrals
int Wavefunction_multiply(T* wavefunction_1, T* wavefunction_2, T* probabilities,
unsigned int lenght_order, T d_x, T d_y, T d_z);
int Wavefunction_multiply(T* wavefunction_1, T* wavefunction_2, T* probabilities, unsigned int lenght_order);
int Wavefunction_relative_lenghts_generate(T* reverse_relative_lenghts, unsigned int lenght_order);
int Wavefunction_square(T* wavefunction_1, T* probabilities, unsigned int lenght_order);
T Probabilities_lenght(T* Probabilities, unsigned int lenght_order,
unsigned int x_range, unsigned int y_range, unsigned int z_range, int x, int y, int z);
T Probabilities_lenght(T* probabilities, unsigned int lenght_order,
unsigned int x_range, unsigned int y_range, unsigned int z_range);
int Probabilities_thread(T* Probabilities, unsigned int lenght_order,
unsigned int x_range, unsigned int y_range, unsigned int z_range, T* lenght);

int Gradient_thread(T* Gradient_1, T* wavefunction_2, unsigned int lenght_order);
int Integral_overlap(T* Wavefunction_1, T* Wavefunction_2, T* result, unsigned int lenght_order, unsigned int x_range_1, unsigned int x_range_2, unsigned int y_range_1, unsigned int y_range_2, unsigned int z_range_1, unsigned int z_range_2,
T x, T y, T z);
int Integrate_Integral_overlap(T* wavefunction_1, T* wavefunction_2, T* result,
unsigned int lenght_order, T x, T y, T z);
int Integral_coulombic(T radius_1, T radius_2, T distance, T* result, bool spin_bonded);
int Integrate_Integral_coulombic(T* density_1, T* density_2, T* result, unsigned int lenght_order, T x, T y, T z,
unsigned int x_range_1, unsigned int x_range_2, unsigned int y_range_1, unsigned int y_range_2,
unsigned int z_range_1, unsigned int z_range_2);
int Integral_nucleus_atraction(T probabilities_lenght, T multiplier, T* result, T* lenght, unsigned int Z);
int Integrate_Integral_nucleus_atraction(T* probabilities,
T* result, T* lenght, unsigned int lenght_order, unsigned int x_range, unsigned int y_range, unsigned int z_range,
T lenght_x, T lenght_y, T lenght_z, unsigned int Z);
int Integral_kinetic(T* Gradient_1, T* Gradient_2, T* result,
unsigned int lenght_order, T d_x, T d_y, T d_z, unsigned int x_range_1, unsigned int x_range_2,
unsigned int y_range_1, unsigned int y_range_2, unsigned int z_range_1, unsigned int z_range_2);
T Rydberg_energy(unsigned int Z, unsigned int n);
T Orbital_moment_energy(int m, T* B0);
T Spin_moment_energy(T s, T B0);
T Orbital_magnetic_field(T potential_energy, T radius, int l);
// Section 4 - generating lists of electrons
int Quantum_numbers_to_orbitals(unsigned int n, unsigned int l, int fulness, atom_orbitals* atom_orbitals_PTR);
int Atoms_to_valence_orbitals(string atom, atom_orbitals* atom_orbitals_PTR);
int Create_atomic_wavefunctions(atom_orbitals *atom_orbitals_PTR, atom_wavefunctions *atom_wavefunctions,
unsigned int size_order,T x,T y,T z);
int Create_bond_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions_1, atom_wavefunctions *atom_wavefunctions_2,
unsigned int count, T electronegativity_1, T electronegativity_2, T x_diference, T y_diference, T z_diference);
int Sum_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions_1, atom_wavefunctions *atom_wavefunctions_2);
// section 5 - generating matrices of integrals and Fock matrices
int Create_nuclear_atraction_integral_matrix(T* matrix, T* nucleuses, unsigned int order, atom_wavefunctions *atom_wavefunctions);
int Create_coulombic_integral_matrix(T* matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions,
small_atom_wavefunctions *small_atom_wavefunctions);
int Create_overlap_integral_matrix(T* matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions);
int Calculate_resonance_integral_matrix(T* overlap_matrix, T* overlap_effective_lenght_integral_matrix, T* resonance_integral_matrix,
unsigned int order, atom_wavefunctions *atom_wavefunctions, small_atom_wavefunctions *small_atom_wavefunctions);
int Calculate_kinetic_integral_matrix(T* matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions);
int Calculate_basis_set_matrix(T* nuclear_atraction_integral_matrix, T* coulombic_integral_matrix, T* resonance_integral_matrix,
T* kinetic_integral_matrix, T* basis_set_matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions);
int Calculate_corr_basis_set_matrix(T* basis_set_matrix, T* correction_matrix, T* corr_basis_set_matrix,
unsigned int order);
T Solve_basis_set_matrix(T* basis_set_matrix, T* overlap_integral_matrix, unsigned int order,
vector<T>* values, atom_wavefunctions *atom_wavefunctions);
int Generate_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions,
small_atom_wavefunctions *small_atom_wavefunctions, unsigned int size_order, bool alocate, bool compute_densities);
public:
// section 6 - completing basis set method and user interface handling
int Atom_orbitals_generate(string UI_input, atom_orbitals *atom_orbitals_PTR);
T Nucleus_repulsive_energy(atom_wavefunctions *atom_wavefunctions);

int String_to_list_electrons(string UI_input, unsigned int size_order,
bool extern_coordinates, vector<T>* x_2, vector<T>* y_2, vector<T>* z_2);
T Calculate(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool alocate,
vector<T>* values);
int Create_excitation(unsigned int electron_number,
unsigned int n, unsigned int l, int m, T spin, bool generate);
int Clear();
~Slater_basis_set_calculations();};
#include "Slater_basis_set_calculations.cc"
# endif
/* Slater_basis_set_calculations_H
Author of this source code Ing. Pavel Florian Ph.D. licensed this source code under the the Apache License:
Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/
*/
