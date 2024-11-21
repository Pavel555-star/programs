#ifndef WAVEFUNCTION_CALCULATIONS_H
#define WAVEFUNCTION_CALCULATIONS_H
#include <iostream>
#include <string>
#include <cmath>
#include <thread> // compiler parameters: -pthread -ffast-math -fno-finite-math-only
#include <vector>
#include <bits/stdc++.h>
#include "basis_set_calculations.h"
// VQE section, Qiskit circuit code
#include "MicroQiskitCpp.h" // End of Qiskit circuit code,  End of VQE section
using namespace std;

template <typename T>
class basis_set_calculations_DFT: public basis_set_calculations<T>
{
public:
    unsigned int count_atoms = 0;
    bool allocation_memory = true;
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
    // Density functional theory section
    vector <T*> atoms_electron_densities;
    vector <T*> atoms_spin_densities;
    vector <T*> atoms_Ksi_densities;
    vector <T*> atoms_gradients_densities;
    vector <T> pre_PBE_wavefunction_lenght_multipliers;
    vector <T> exchange_energies;
    vector <T> correlation_energies;
    vector <T> electron_energies;
    unsigned int PBE_iterations = 2; // Number of iterations of including the correlation energy to main calculations
    array<T, 118> Wigner_Seitz_radiuses = {
    1.7794259739938E-10, 1.85089784636951E-10, 1.73099989649705E-10, 9.88433935284313E-11, 8.82013859640198E-11, 8.06137530478938E-11,
    1.11120028851942E-10, 9.74549532475575E-11, 8.58972752524277E-11, 1.04004524925882E-10, 2.11090198754033E-10, 1.40427965161849E-10,
    1.09675312651626E-10, 1.06087905940566E-10, 9.69777732449294E-11, 1.00743770215178E-10, 1.08670994417806E-10, 1.03498633350404E-10,
    2.59121991588759E-10, 1.72357333264976E-10, 1.25744420684479E-10, 1.01697433727408E-10, 8.7612312122663E-11, 7.82983859721842E-11,
    7.55431254547781E-11, 7.05893047394917E-11, 6.62965314287656E-11, 6.39025227414756E-11, 6.34893657365994E-11, 6.71073572830079E-11,
    7.10989839420626E-11, 7.28135314589623E-11, 7.01634901522086E-11, 7.40677258652833E-11, 8.43441641421331E-11, 9.14217630814674E-11,
    2.80605537309278E-10, 1.87314512234621E-10, 1.37938049287963E-10, 1.11476200596967E-10, 9.50487816563437E-11, 8.50936086726222E-11,
    7.9584183900951E-11, 7.38757352170976E-11, 7.14589146593838E-11, 7.05101739568561E-11, 7.18053389233462E-11, 7.54140277163509E-11,
    7.8216443813395E-11, 7.69680510485447E-11, 7.82995640650464E-11, 7.96948671471751E-11, 8.43115896785351E-11, 9.76751547507002E-11,
    3.01007150707118E-10, 1.97873631293056E-10, 1.43775855407312E-10, 1.27016850712164E-10, 1.18133648359178E-10, 1.10744578030526E-10,
    1.04158031655681E-10, 9.96605195244992E-11, 1.08303483550119E-10, 9.23782018863133E-11, 8.85894221803847E-11, 8.56353479157359E-11,
    8.2992288759032E-11, 8.05125469613091E-11, 7.82201170636002E-11, 8.53072299237178E-11, 7.45425334560279E-11, 6.65737100155775E-11,
    6.09110618935683E-11, 5.73994958045938E-11, 5.5076074987394E-11, 5.33205985200425E-11, 5.27420012606784E-11, 5.31425865053357E-11,
    5.44770252171394E-11, 6.08907845425743E-11, 6.32464394290588E-11, 6.36973266212036E-11, 6.63291682577297E-11, 6.64627306585628E-11,
    7.44554338152102E-11, 8.54735529189076E-11, 3.37861036869392E-10, 2.07627303324919E-10, 1.43876019309543E-10, 1.25138071717085E-10,
    1.0597965475137E-10, 9.38455672697644E-11, 8.68563503438371E-11, 8.47896675065719E-11, 9.18938610003413E-11, 8.9423555928356E-11,
    8.44225102834228E-11, 8.39530470913391E-11, 9.54076703541614E-11, 9.08051059645699E-11, 8.71255215815781E-11, 8.6399809048525E-11,
    7.51181602534702E-11, 6.32892809067158E-11, 6.37124243526582E-11, 6.06956641748638E-11, 5.74879645457555E-11, 5.51352496590779E-11,
    5.02749365992974E-11, 5.55851465462774E-11, 5.36653161612468E-11, 6.76855780795347E-11, 7.05765446808478E-11, 5.70629246965641E-11,
    6.62931350236222E-11, 6.69320102285978E-11, 8.01313714554136E-11, 7.965320751222E-11};
    T beta = 0.066725; // scalable argument
    T gamma = (1 - log(2))/(this->Pi * this->Pi); // equal 0.031091, scalable argument to 0.025
    T kappa = 0.804;
    T mi = (beta * this->Pi * this->Pi)/3.00;
    T omega = 0.046644; // scalable argument
    // End of density functional theory section
    // VQE section
    vector<T> VQE_Wavefunction_lenght_multipliers_1;
    vector<T> VQE_Wavefunction_lenght_multipliers_2;
    vector<T> VQE_Eigenvectors_1;
    vector<T> VQE_Eigenvectors_2;
    vector<bool> VQE_correlation_energy_sign;
    vector<double> VQE_ansatz;
    T VQE_Hamiltonian;
    T VQE_previous_Hamiltonian;
    unsigned int VQE_iterations = 0;
    unsigned int VQE_max_iterations = 2;
    unsigned int VQE_shots = 1; // Simulator settings
    unsigned int VQE_max_connections = 16;
    T VQE_interconnection_distance = 5;
    bool Bravyi_Kitaev = false;
    vector<pair<unsigned int, unsigned int>> VQE_gates; // Qiskit circuit code
    string circuit_string;
    QuantumCircuit circuit;
    // End of Qiskit circuit code,  End of VQE section
private:
    T Detect_center_of_mass();
    T Detect_dipole_moments(dipole_moment *dipole_moment);
    T Create_crystal_field(central_cations *central_cations);
    // Density functional theory section
    T Alocate_densities(unsigned int size_order);
    T Compute_density_thread(T** densities, T* atom_density, unsigned int begin, unsigned int end, unsigned int atom_size,
    vector<T> spins, vector<int> spin_paired, bool spin_density);
    T Compute_densities(vector<T*> densities, vector<T*> atoms_densities, vector<unsigned int> index,
    unsigned int size_order, vector<T> spins, vector<int> spin_paired, bool spin_density);
    T Correct_densities(vector<T*> atoms_densities_list, vector<T>  x, vector <T> y, vector <T> z, unsigned int size_order);
    T Compute_Ksi_thread(T* spin_density, T* electron_density, T* Ksi_density, unsigned int size_order);
    T Compute_gradient_thread(T* electron_density, T* gradient_density, unsigned int size_order);
    T Compute_Ksi_and_gradients(vector<T*> atoms_electron_densities, vector<T*> atoms_spin_densities,
    vector<T*> atoms_Ksi_densities, vector<T*> atoms_gradients_densities, unsigned int size_order);
    T PBE_thread(T* atoms_gradients_density, T* atoms_Ksi_density, T* electron_density,
    unsigned int size_order, T* exchange_energy, T* correlation_energy, T Wigner_Seitz_radius);
    T PBE_compute(vector<T*> atoms_gradients_densities, vector<T*> atoms_Ksi_densities, vector<T*> electron_densities,
    unsigned int size_order);
    // End of density functional theory section
    // VQE section
    T Set_circuit();
    T Run_circuit(vector<double>& ansatz);
    // End of VQE section
public:
    T String_to_advanced_parameters(string UI_input, unsigned int size_order,
    bool extern_coordinates, vector<T>* x_2, vector<T>* y_2, vector<T>* z_2);
    T Execute_calculation(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool dealocate,
    vector<T>* values, vector<T>* spin_density_vector,  vector<T>* spin_values);
    // Density functional theory section
    T Execute_PBE(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool dealocate,
    vector<T>* values, vector<T>* spin_density_vector, vector<T>* spin_values);
    // End of density functional theory section
    // VQE section
    T Execute_PBE_VQE(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool dealocate,
    vector<T>* values, vector<T>* spin_density_vector, vector<T>* spin_values);
    // End of VQE section
    T Calculate_Huckel_Matrix(T* Huckel_matrix, unsigned int* Huckel_matrix_order, vector<unsigned int> atom_numbers);
    T Detect_symetry_information(symetry_axes *symetry_axes, symetry_planes *symetry_planes);
    T D3_vector_multiply(T* a, T* b, T* c);
    T Clear();
    ~basis_set_calculations_DFT();
};
#include "basis_set_calculations_DFT.cc"
#endif // WAVEFUNCTION_CALCULATIONS_H
/*
Author of this source code Ing. Pavel Florian Ph.D. licensed this source code under the the Apache License:
Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/  */ 
