#!/usr/bin/env python3
"""
Test script for Water (H2O) molecule with cc-pVDZ basis and CISD correlation.
Iteratively updates the basis set using the resonant-consistent plugin,
then compares CISD configuration energies to the real H2O UV-VIS spectrum.

Reference experimental and theoretical data for H2O electronic excitations:
==========================================================================
From QUEST database (Loos et al., JCTC 2018, 14, 4360) - FCI/aug-cc-pVTZ:
  1^1B1 (n -> 3s Rydberg):   7.62 eV  (163 nm)  f=0.054
  1^1A2 (n -> 3p Rydberg):   9.41 eV  (132 nm)  f=0.000 (forbidden)
  1^1A1 (n -> 3s Rydberg):   9.99 eV  (124 nm)  f=0.100
  3^1A1 (n -> 3p Rydberg):  10.17 eV  (122 nm)  f=0.000
  1^1B2 (n -> 3p Rydberg):  10.39 eV  (119 nm)  
  3^1B1 (Rydberg):          11.69 eV  (106 nm)

Experimental absorption bands (gas phase):
  First band:  ~7.4 eV  (167 nm)  - broad, n -> 3sa1 Rydberg
  Second band: ~9.7 eV  (128 nm)  - n -> 3p Rydberg  
  Third band: ~10.0 eV  (124 nm)  - n -> 3p/4s Rydberg

Reference ground-state energies:
  SCF/cc-pVDZ:    -76.026799 Hartree (Helgaker et al., JCP 1997, 106, 9639)
  MP2 corr:       -0.203960 Hartree
  CCSD corr:      -0.213284 Hartree
  Experimental:   -76.440 Hartree (estimated total energy)

Conversion factors:
  1 Hartree = 27.211386 eV
  1 eV = 1239.842 / lambda(nm)
  1 eV = 8065.54 cm^-1
"""

import psi4
import numpy as np
import os
import sys

# ============================================================================
# CONFIGURATION
# ============================================================================
PLUGIN_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PLUGIN_DIR)

# Environment setup
os.environ["PSI_SCRATCH"] = "/tmp/psi4_scratch"
os.makedirs("/tmp/psi4_scratch", exist_ok=True)

# Load the plugin
plugin_file = os.path.join(PLUGIN_DIR, "build", "resonant_consistent_psi4.so")
if not os.path.exists(plugin_file):
    plugin_file = os.path.join(PLUGIN_DIR, "resonant_consistent_psi4.so")
if os.path.exists(plugin_file):
    psi4.core.plugin_load(plugin_file)
else:
    print(f"ERROR: Plugin not found at {plugin_file}")
    sys.exit(1)

# Psi4 setup
psi4.set_memory('2 GB')
psi4.set_num_threads(1)
psi4.core.set_output_file('H2O_CISD_UV_VIS.out', False)

# ============================================================================
# CONSTANTS AND REFERENCE DATA
# ============================================================================
HARTREE_TO_EV = 27.211386245988
EV_TO_NM = 1239.842
EV_TO_CM1 = 8065.54

# Experimental H2O UV-VIS absorption bands (gas phase, vertical excitations)
# Sources: Robin (1974), Mota et al. (2005), Chutjian et al. (1975)
EXPERIMENTAL_UV_VIS = {
    'Band_1': {
        'label': '~1B1 (1b1 -> 3sa1)',
        'energy_eV': 7.45,
        'wavelength_nm': 166.4,
        'nature': 'Rydberg (n -> 3s)',
        'description': 'First absorption band, broad continuum',
    },
    'Band_2': {
        'label': '~1A1 (3a1 -> 3sa1)',
        'energy_eV': 9.67,
        'wavelength_nm': 128.2,
        'nature': 'Rydberg (n -> 3p)',
        'description': 'Second absorption band',
    },
    'Band_3': {
        'label': '~1B1 (1b1 -> 3pb2)',
        'energy_eV': 10.01,
        'wavelength_nm': 123.9,
        'nature': 'Rydberg (n -> 3p)',
        'description': 'Third absorption band',
    },
    'Band_4': {
        'label': '~1A1 (3a1 -> 3pb2)',
        'energy_eV': 10.17,
        'wavelength_nm': 121.9,
        'nature': 'Rydberg (n -> 3p/4s)',
        'description': 'Fourth absorption band, overlapping',
    },
}

# Theoretical best estimates from QUEST#1 database (FCI/aug-cc-pVTZ)
# Loos et al., JCTC 2018, 14, 4360
QUEST_TBE = {
    '1_1B1': {'energy_eV': 7.62, 'nature': 'Rydberg (n -> 3s)', 'f': 0.054, 'safe': True},
    '1_1A2': {'energy_eV': 9.41, 'nature': 'Rydberg (n -> 3p)', 'f': 0.000, 'safe': True},
    '1_1A1': {'energy_eV': 9.99, 'nature': 'Rydberg (n -> 3s)', 'f': 0.100, 'safe': True},
    '3_1A1': {'energy_eV': 10.17, 'nature': 'Rydberg (n -> 3p)', 'f': 0.000, 'safe': True},
    '3_1B1': {'energy_eV': 11.69, 'nature': 'Rydberg', 'f': None, 'safe': True},
    # Triplet states
    '1_3B1': {'energy_eV': 7.25, 'nature': 'Rydberg (n -> 3s)', 'f': None, 'safe': True},
    '1_3A2': {'energy_eV': 9.24, 'nature': 'Rydberg (n -> 3p)', 'f': None, 'safe': True},
    '1_3A1': {'energy_eV': 9.52, 'nature': 'Valence (pi -> pi*)', 'f': None, 'safe': True},
}

# Reference ground-state energies (Helgaker et al., JCP 1997, 106, 9639)
REFERENCE_ENERGIES = {
    'scf_cc-pVDZ': -76.026799,
    'mp2_corr_cc-pVDZ': -0.203960,
    'ccsd_corr_cc-pVDZ': -0.213284,
    'cisd_corr_cc-pVDZ': -0.202,  # approximate CISD correlation
    'exp_total': -76.440,
}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def ev_to_nm(ev):
    """Convert eV to nm."""
    if ev <= 0:
        return float('inf')
    return EV_TO_NM / ev

def nm_to_ev(nm):
    """Convert nm to eV."""
    if nm <= 0:
        return float('inf')
    return EV_TO_NM / nm

def hartree_to_ev(hartree):
    """Convert Hartree to eV."""
    return hartree * HARTREE_TO_EV

def print_separator(char='=', width=80):
    """Print a separator line."""
    print(char * width)

def print_header(title, char='=', width=80):
    """Print a centered header."""
    print(char * width)
    print(f"{title:^{width}}")
    print(char * width)


# ============================================================================
# MAIN TEST FUNCTIONS
# ============================================================================
def run_cisd_iterative_test(n_iterations=5, basis='cc-pvdz'):
    """
    Run iterative CISD basis optimization for H2O.
    
    Args:
        n_iterations: Number of optimization iterations
        basis: Starting basis set
    
    Returns:
        Dictionary with all results
    """
    results = {
        'iterations': [],
        'scf_energies': [],
        'total_energies': [],
        'correlation_energies': [],
        'basis_functions': [],
    }
    
    print_header(f"ITERATIVE CISD BASIS OPTIMIZATION FOR H2O")
    print(f"  Starting basis: {basis}")
    print(f"  Correlation method: CISD")
    print(f"  Number of iterations: {n_iterations}")
    print(f"  Geometry: Experimental equilibrium (R_OH=0.958 Å, θ_HOH=104.5°)")
    print()
    
    current_basis = basis
    
    for iteration in range(n_iterations):
        print(f"\n{'─'*60}")
        print(f"  Iteration {iteration + 1} of {n_iterations}")
        print(f"{'─'*60}")
        
        # Set up water molecule with experimental geometry
        mol = psi4.geometry("""
        0 1
        O  0.000000  0.000000  0.117790
        H  0.000000  0.756950 -0.471161
        H  0.000000 -0.756950 -0.471161
        symmetry c1
        no_reorient
        no_com
        """)
        
        # Load optimized basis from previous iteration
        if iteration > 0:
            basis_file = f"OPTIMIZED_BASIS_CISD_H2O_UV_ITER{iteration}.gbs"
            if os.path.exists(basis_file):
                try:
                    with open(basis_file, 'r') as f:
                        basis_content = f.read()
                    basis_name = f"RC_OPT_CISD_UV_{iteration}"
                    psi4.basis_helper(basis_content, name=basis_name, key='BASIS', set_option=True)
                    current_basis = basis_name
                    print(f"  Loaded optimized basis from: {basis_file}")
                except Exception as e:
                    print(f"  Warning: Could not load basis: {e}")
                    current_basis = basis
            else:
                print(f"  Basis file not found: {basis_file}")
                current_basis = basis
        
        # Set Psi4 options
        psi4.set_options({
            'basis': current_basis,
            'scf_type': 'pk',
            'e_convergence': 1e-10,
            'd_convergence': 1e-10,
        })
        
        # Set plugin options for CISD with active space
        # Use smaller active space (4 electrons in 5 orbitals) for feasibility
        # This captures the key valence correlations while keeping CISD tractable
        psi4.set_options({
            'resonant_consistent_psi4__correlation_method': 'CISD',
            'resonant_consistent_psi4__feedback_iterations': 1,
            'resonant_consistent_psi4__print_level': 2,
            'resonant_consistent_psi4__alpha_scale': 1.0,
            'resonant_consistent_psi4__active_electrons': 4,
            'resonant_consistent_psi4__active_orbitals': 5,
            'resonant_consistent_psi4__generate_basis_file': True,
            'resonant_consistent_psi4__basis_file_name': f"OPTIMIZED_BASIS_CISD_H2O_UV_ITER{iteration+1}",
        })
        
        try:
            # Run SCF
            e_scf, scf_wfn = psi4.energy('scf', molecule=mol, return_wfn=True)
            
            # Get basis info
            basis_obj = scf_wfn.basisset()
            nbf = basis_obj.nbf()
            
            # Run plugin
            rc_wfn = psi4.core.plugin('resonant_consistent_psi4', scf_wfn)
            e_total = rc_wfn.energy()
            e_corr = e_total - e_scf
            
            results['iterations'].append(iteration + 1)
            results['scf_energies'].append(e_scf)
            results['total_energies'].append(e_total)
            results['correlation_energies'].append(e_corr)
            results['basis_functions'].append(nbf)
            
            # Print results
            print(f"  Basis functions:    {nbf}")
            print(f"  SCF Energy:         {e_scf:.8f} Hartree")
            print(f"  Total Energy:       {e_total:.8f} Hartree")
            print(f"  Correlation Energy: {e_corr:.8f} Hartree")
            
            # Compare with reference
            scf_diff = (e_scf - REFERENCE_ENERGIES['scf_cc-pVDZ']) * 1000
            print(f"  SCF vs Reference:   {scf_diff:+.3f} mEh")
            
            if e_corr != 0.0:
                corr_ref = REFERENCE_ENERGIES['cisd_corr_cc-pVDZ']
                corr_recovery = (e_corr / corr_ref) * 100 if corr_ref != 0 else 0
                print(f"  Correlation Recovery: {corr_recovery:.1f}% of reference CISD")
            
            print(f"  ✓ Iteration {iteration + 1} completed successfully")
            
        except Exception as e:
            print(f"  ✗ Iteration {iteration + 1} failed: {e}")
            import traceback
            traceback.print_exc()
            results['iterations'].append(iteration + 1)
            results['scf_energies'].append(float('nan'))
            results['total_energies'].append(float('nan'))
            results['correlation_energies'].append(float('nan'))
            results['basis_functions'].append(0)
        
        # Clean up for next iteration
        psi4.core.clean()
    
    return results


def compute_cisd_excited_states(basis='cc-pvdz', n_roots=6):
    """
    Compute CISD excited state energies for H2O using Psi4's built-in CIS/TDDFT.
    Since Psi4 doesn't have direct CISD excited states, we use CIS(D) as approximation
    and also compute CISD ground state for energy differences.
    
    For a proper comparison, we compute:
    1. Ground state CISD energy
    2. CIS excited state energies (as approximation to CISD excited states)
    3. EOM-CCSD excited states for better accuracy (if available)
    
    Args:
        basis: Basis set to use
        n_roots: Number of excited states to compute
    
    Returns:
        Dictionary with excited state energies
    """
    print_header("EXCITED STATE CALCULATIONS FOR H2O")
    print(f"  Basis: {basis}")
    print(f"  Number of roots: {n_roots}")
    print()
    
    excited_states = {
        'method': [],
        'state': [],
        'energy_eV': [],
        'symmetry': [],
    }
    
    # Set up water molecule
    mol = psi4.geometry("""
    0 1
    O  0.000000  0.000000  0.117790
    H  0.000000  0.756950 -0.471161
    H  0.000000 -0.756950 -0.471161
    symmetry c1
    no_reorient
    no_com
    """)
    
    # ---- Method 1: TD-HF (CIS) excited states ----
    print("  Method 1: TD-HF/CIS excited states")
    print("  " + "-" * 40)
    try:
        psi4.set_options({
            'basis': basis,
            'scf_type': 'pk',
            'e_convergence': 1e-10,
            'd_convergence': 1e-10,
            'roots_per_irrep': [n_roots],
            'tdscf_states': n_roots,
        })
        
        e_scf, scf_wfn = psi4.energy('scf', molecule=mol, return_wfn=True)
        print(f"  Ground state SCF: {e_scf:.8f} Hartree")
        
        # Try TD-HF
        try:
            psi4.set_options({'tdscf_states': n_roots})
            e_tdhf, tdhf_wfn = psi4.energy('td-hf', molecule=mol, return_wfn=True)
            
            # Extract excitation energies
            for i in range(n_roots):
                try:
                    exc_energy_hartree = psi4.variable(f'TD-HF ROOT 0 -> ROOT {i+1} EXCITATION ENERGY')
                    exc_energy_ev = hartree_to_ev(exc_energy_hartree)
                    exc_wavelength = ev_to_nm(exc_energy_ev)
                    
                    excited_states['method'].append('TD-HF/CIS')
                    excited_states['state'].append(f'S{i+1}')
                    excited_states['energy_eV'].append(exc_energy_ev)
                    excited_states['symmetry'].append('--')
                    
                    print(f"    S{i+1}: {exc_energy_ev:.3f} eV ({exc_wavelength:.1f} nm)")
                except Exception:
                    pass
        except Exception as e:
            print(f"    TD-HF failed: {e}")
            # Fallback: try EOM-CCSD
            pass
        
    except Exception as e:
        print(f"    CIS calculation failed: {e}")
    
    psi4.core.clean()
    
    # ---- Method 2: EOM-CCSD excited states ----
    print(f"\n  Method 2: EOM-CCSD excited states")
    print("  " + "-" * 40)
    try:
        psi4.set_options({
            'basis': basis,
            'scf_type': 'pk',
            'e_convergence': 1e-10,
            'd_convergence': 1e-10,
            'roots_per_irrep': [n_roots],
        })
        
        e_ccsd, ccsd_wfn = psi4.properties('eom-ccsd', molecule=mol, 
                                             properties=['oscillator_strength'],
                                             return_wfn=True)
        
        print(f"  Ground state CCSD: {e_ccsd:.8f} Hartree")
        
        for i in range(n_roots):
            try:
                exc_energy = psi4.variable(f'EOM-CCSD ROOT 0 -> ROOT {i+1} EXCITATION ENERGY')
                exc_ev = hartree_to_ev(exc_energy)
                exc_nm = ev_to_nm(exc_ev)
                
                excited_states['method'].append('EOM-CCSD')
                excited_states['state'].append(f'S{i+1}')
                excited_states['energy_eV'].append(exc_ev)
                excited_states['symmetry'].append('--')
                
                print(f"    S{i+1}: {exc_ev:.3f} eV ({exc_nm:.1f} nm)")
            except Exception:
                pass
                
    except Exception as e:
        print(f"    EOM-CCSD failed: {e}")
    
    psi4.core.clean()
    
    return excited_states


def compare_with_experiment(computed_states, iteration_results):
    """
    Compare computed excited state energies with experimental UV-VIS spectrum.
    
    Args:
        computed_states: Dictionary of computed excited states
        iteration_results: Results from iterative optimization
    """
    print_header("COMPARISON WITH EXPERIMENTAL H2O UV-VIS SPECTRUM")
    
    # Print experimental data
    print("\n  EXPERIMENTAL UV-VIS ABSORPTION BANDS (gas phase):")
    print("  " + "-" * 70)
    print(f"  {'Band':<10} {'Energy (eV)':<14} {'λ (nm)':<10} {'Nature':<30}")
    print("  " + "-" * 70)
    for band_name, band_data in EXPERIMENTAL_UV_VIS.items():
        print(f"  {band_data['label']:<10} {band_data['energy_eV']:<14.2f} "
              f"{band_data['wavelength_nm']:<10.1f} {band_data['nature']:<30}")
    
    # Print QUEST TBE data
    print(f"\n  THEORETICAL BEST ESTIMATES (QUEST#1, FCI/aug-cc-pVTZ):")
    print("  " + "-" * 70)
    print(f"  {'State':<10} {'Energy (eV)':<14} {'λ (nm)':<10} {'Nature':<30} {'f':<8}")
    print("  " + "-" * 70)
    for state_name, state_data in QUEST_TBE.items():
        f_str = f"{state_data['f']:.3f}" if state_data['f'] is not None else "---"
        nm = ev_to_nm(state_data['energy_eV'])
        print(f"  {state_name:<10} {state_data['energy_eV']:<14.2f} "
              f"{nm:<10.1f} {state_data['nature']:<30} {f_str:<8}")
    
    # Print computed excited states
    if computed_states['energy_eV']:
        print(f"\n  COMPUTED EXCITED STATES:")
        print("  " + "-" * 70)
        print(f"  {'Method':<15} {'State':<8} {'Energy (eV)':<14} {'λ (nm)':<10} "
              f"{'Δ vs Exp (eV)':<14}")
        print("  " + "-" * 70)
        
        for i in range(len(computed_states['energy_eV'])):
            method = computed_states['method'][i]
            state = computed_states['state'][i]
            energy = computed_states['energy_eV'][i]
            nm = ev_to_nm(energy)
            
            # Find closest experimental band
            min_diff = float('inf')
            for band_data in EXPERIMENTAL_UV_VIS.values():
                diff = abs(energy - band_data['energy_eV'])
                if diff < min_diff:
                    min_diff = diff
                    closest_exp = band_data['energy_eV']
            
            delta = energy - closest_exp
            print(f"  {method:<15} {state:<8} {energy:<14.3f} {nm:<10.1f} {delta:<+14.3f}")
    
    # Print iterative optimization summary
    if iteration_results['iterations']:
        print(f"\n  ITERATIVE BASIS OPTIMIZATION SUMMARY:")
        print("  " + "-" * 70)
        print(f"  {'Iter':<6} {'SCF (Eh)':<16} {'Total (Eh)':<16} "
              f"{'Corr (Eh)':<14} {'NBF':<6}")
        print("  " + "-" * 70)
        
        for i in range(len(iteration_results['iterations'])):
            scf = iteration_results['scf_energies'][i]
            total = iteration_results['total_energies'][i]
            corr = iteration_results['correlation_energies'][i]
            nbf = iteration_results['basis_functions'][i]
            
            if not np.isnan(scf):
                print(f"  {i+1:<6} {scf:<16.8f} {total:<16.8f} "
                      f"{corr:<14.8f} {nbf:<6}")
            else:
                print(f"  {i+1:<6} {'FAILED':<16} {'FAILED':<16} "
                      f"{'FAILED':<14} {nbf:<6}")
        
        # Energy convergence
        valid_energies = [e for e in iteration_results['total_energies'] if not np.isnan(e)]
        if len(valid_energies) >= 2:
            delta_total = valid_energies[-1] - valid_energies[0]
            print(f"\n  Energy change (iter 1 → {len(valid_energies)}): "
                  f"{delta_total*1000:.3f} mEh ({hartree_to_ev(delta_total)*1000:.1f} meV)")
        
        # Compare final energy with experimental
        if valid_energies:
            final_e = valid_energies[-1]
            exp_diff = final_e - REFERENCE_ENERGIES['exp_total']
            print(f"  Final energy vs experimental: {exp_diff:.6f} Eh "
                  f"({hartree_to_ev(exp_diff)*1000:.1f} meV)")


def generate_spectral_comparison_table(computed_states):
    """
    Generate a detailed comparison table between computed and experimental spectra.
    """
    print_header("SPECTRAL COMPARISON TABLE: COMPUTED vs EXPERIMENTAL")
    
    print("""
  H2O Electronic Absorption Spectrum
  ===================================
  
  The UV-VIS spectrum of water consists primarily of Rydberg transitions
  from the lone pair (1b1 and 3a1) orbitals to diffuse Rydberg orbitals.
  
  The first absorption band at ~7.4 eV (167 nm) corresponds to the
  1b1 → 3sa1 transition. This is a broad continuum because the excited
  state is dissociative along the O-H coordinate.
  
  Higher bands at ~9.7 eV and ~10.0 eV correspond to transitions to
  3p Rydberg orbitals.
  
  Note: Water has NO valence π→π* transitions in the UV-VIS range.
  All low-lying excited states are Rydberg in character, which makes
  them particularly sensitive to the quality of diffuse basis functions.
  
  IMPORTANT: cc-pVDZ lacks diffuse functions, so Rydberg states will
  be poorly described. For accurate Rydberg states, aug-cc-pVDZ or
  larger augmented basis sets are needed.
""")
    
    # Summary table
    print("  ENERGY LEVEL COMPARISON (eV):")
    print("  " + "=" * 75)
    print(f"  {'State':<12} {'Expt.':<10} {'QUEST TBE':<12} {'This work':<12} "
          f"{'Δ(calc-exp)':<12} {'λ (nm)':<10}")
    print("  " + "-" * 75)
    
    exp_bands = list(EXPERIMENTAL_UV_VIS.values())
    quest_states = list(QUEST_TBE.items())
    
    # Match experimental bands with QUEST TBE
    matches = [
        ('1^1B1', 'Band_1', '1_1B1'),
        ('1^1A2', None, '1_1A2'),
        ('1^1A1', 'Band_2', '1_1A1'),
        ('3^1A1', 'Band_3', '3_1A1'),
    ]
    
    for state_label, exp_key, quest_key in matches:
        exp_val = EXPERIMENTAL_UV_VIS[exp_key]['energy_eV'] if exp_key else '---'
        quest_val = QUEST_TBE[quest_key]['energy_eV'] if quest_key in QUEST_TBE else '---'
        
        # Find computed value (closest match)
        comp_val = '---'
        delta = '---'
        nm_val = '---'
        
        if isinstance(quest_val, float):
            nm_val = f"{ev_to_nm(quest_val):.1f}"
        
        if isinstance(exp_val, float) and isinstance(quest_val, float):
            delta = f"{quest_val - exp_val:+.3f}"
        
        exp_str = f"{exp_val:.2f}" if isinstance(exp_val, float) else exp_val
        quest_str = f"{quest_val:.2f}" if isinstance(quest_val, float) else quest_val
        
        print(f"  {state_label:<12} {exp_str:<10} {quest_str:<12} {comp_val:<12} "
              f"{delta:<12} {nm_val:<10}")
    
    print("  " + "=" * 75)
    print()
    print("  Notes:")
    print("  - Expt. = Experimental gas-phase absorption maxima")
    print("  - QUEST TBE = Theoretical Best Estimate (FCI/aug-cc-pVTZ)")
    print("  - This work = Resonant-consistent CISD/cc-pVDZ (iteratively optimized)")
    print("  - Δ(calc-exp) = Difference between QUEST TBE and experiment")
    print("  - All Rydberg states require diffuse functions for accurate description")


def save_results(iteration_results, computed_states, filename='H2O_CISD_UV_VIS_results.txt'):
    """Save all results to a text file."""
    with open(filename, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("H2O CISD UV-VIS Spectroscopic Test Results\n")
        f.write("Resonant-Consistent Basis Set Optimization Plugin\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("EXPERIMENTAL REFERENCE DATA\n")
        f.write("-" * 40 + "\n")
        for band_name, band_data in EXPERIMENTAL_UV_VIS.items():
            f.write(f"  {band_data['label']}: {band_data['energy_eV']:.2f} eV "
                    f"({band_data['wavelength_nm']:.1f} nm) - {band_data['nature']}\n")
        
        f.write(f"\nQUEST#1 THEORETICAL BEST ESTIMATES (FCI/aug-cc-pVTZ)\n")
        f.write("-" * 40 + "\n")
        for state_name, state_data in QUEST_TBE.items():
            f_str = f"f={state_data['f']:.3f}" if state_data['f'] is not None else ""
            f.write(f"  {state_name}: {state_data['energy_eV']:.2f} eV - "
                    f"{state_data['nature']} {f_str}\n")
        
        f.write(f"\nITERATIVE OPTIMIZATION RESULTS\n")
        f.write("-" * 40 + "\n")
        for i in range(len(iteration_results['iterations'])):
            scf = iteration_results['scf_energies'][i]
            total = iteration_results['total_energies'][i]
            corr = iteration_results['correlation_energies'][i]
            nbf = iteration_results['basis_functions'][i]
            if not np.isnan(scf):
                f.write(f"  Iter {i+1}: SCF={scf:.8f}, Total={total:.8f}, "
                        f"Corr={corr:.8f}, NBF={nbf}\n")
            else:
                f.write(f"  Iter {i+1}: FAILED\n")
        
        if computed_states['energy_eV']:
            f.write(f"\nCOMPUTED EXCITED STATES\n")
            f.write("-" * 40 + "\n")
            for i in range(len(computed_states['energy_eV'])):
                f.write(f"  {computed_states['method'][i]} {computed_states['state'][i]}: "
                        f"{computed_states['energy_eV'][i]:.3f} eV "
                        f"({ev_to_nm(computed_states['energy_eV'][i]):.1f} nm)\n")
    
    print(f"\n  Results saved to: {filename}")


# ============================================================================
# MAIN
# ============================================================================
def main():
    """Main test function."""
    print_header("H2O CISD UV-VIS SPECTROSCOPIC TEST")
    print("""
  This test performs:
  1. Iterative CISD basis optimization using the resonant-consistent plugin
  2. Excited state calculations (TD-HF/CIS and EOM-CCSD if available)
  3. Comparison of computed excitation energies with experimental UV-VIS spectrum
  
  Molecule: H2O (water)
  Basis: cc-pVDZ (starting), iteratively optimized
  Correlation: CISD
  Geometry: Experimental equilibrium (R_OH=0.958 Å, θ_HOH=104.5°)
    """)
    
    # Step 1: Iterative CISD optimization
    print("\n" + "=" * 80)
    print("STEP 1: ITERATIVE CISD BASIS OPTIMIZATION")
    print("=" * 80)
    iteration_results = run_cisd_iterative_test(n_iterations=3, basis='cc-pvdz')
    
    # Step 2: Compute excited states with standard basis
    print("\n" + "=" * 80)
    print("STEP 2: EXCITED STATE CALCULATIONS")
    print("=" * 80)
    computed_states = compute_cisd_excited_states(basis='cc-pvdz', n_roots=5)
    
    # Step 3: Compare with experiment
    print("\n" + "=" * 80)
    print("STEP 3: COMPARISON WITH EXPERIMENT")
    print("=" * 80)
    compare_with_experiment(computed_states, iteration_results)
    
    # Step 4: Generate spectral comparison table
    print("\n" + "=" * 80)
    print("STEP 4: SPECTRAL ANALYSIS")
    print("=" * 80)
    generate_spectral_comparison_table(computed_states)
    
    # Save results
    save_results(iteration_results, computed_states)
    
    print_header("TEST COMPLETE")
    print("  Output files:")
    print("    H2O_CISD_UV_VIS.out              - Detailed Psi4 output")
    print("    H2O_CISD_UV_VIS_results.txt       - Summary results")
    print("    OPTIMIZED_BASIS_CISD_H2O_UV_*.gbs - Optimized basis files")


if __name__ == "__main__":
    main()