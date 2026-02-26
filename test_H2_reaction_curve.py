#!/usr/bin/env python3
"""
Test script for H2 reaction curve with iterative basis updating.
Tests the potential energy surface of H2 at various bond distances.
"""

import psi4
import numpy as np
import os
import sys
import json

# Add plugin directory to path
plugin_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, plugin_dir)

# Set up environment
os.environ["PSI_SCRATCH"] = "/tmp/psi4_scratch"
os.makedirs("/tmp/psi4_scratch", exist_ok=True)

# Load the plugin
plugin_file = os.path.join(plugin_dir, "build", "resonant_consistent_psi4.so")
if os.path.exists(plugin_file):
    psi4.core.plugin_load(plugin_file)
else:
    plugin_file = os.path.join(plugin_dir, "resonant_consistent_psi4.so")
    if os.path.exists(plugin_file):
        psi4.core.plugin_load(plugin_file)
    else:
        print(f"ERROR: Plugin not found")
        sys.exit(1)

# Set up Psi4
psi4.set_options({'scf_type': 'pk', 'e_convergence': 1e-8, 'd_convergence': 1e-8})

def get_h2_molecule(bond_length):
    """Create H2 molecule with specified bond length (in Angstrom)."""
    return psi4.geometry(f"""
    0 1
    H  0.000000  0.000000  0.000000
    H  0.000000  0.000000  {bond_length}
    symmetry c1
    """)

def run_single_point(mol, basis_label, method, active_electrons, active_orbitals, 
                     generate_basis=True, basis_file_name="H2_opt"):
    """Run a single point calculation."""
    psi4.set_options({
        'CORRELATION_METHOD': method.upper(),
        'ACTIVE_ELECTRONS': active_electrons,
        'ACTIVE_ORBITALS': active_orbitals,
        'FEEDBACK_ITERATIONS': 10,
        'GENERATE_BASIS_FILE': generate_basis,
        'BASIS_FILE_NAME': basis_file_name
    })
    
    psi4.set_options({'basis': basis_label})
    
    # Run SCF
    e_scf, scf_wfn = psi4.energy('scf', molecule=mol, return_wfn=True)
    
    # Run plugin
    rc_wfn = psi4.core.plugin('resonant_consistent_psi4', scf_wfn)
    e_total = rc_wfn.energy()
    
    return e_scf, e_total

def run_h2_reaction_curve(basis_name, method, bond_lengths, n_iterations=5):
    """
    Run H2 reaction curve with iterative basis updating.
    
    For each bond length:
    1. Start with initial basis
    2. Run n_iterations of basis optimization
    3. Save the optimized basis for that geometry
    
    Args:
        basis_name: Initial basis set name or path
        method: Correlation method
        bond_lengths: List of bond lengths to test (in Angstrom)
        n_iterations: Number of basis optimization iterations per geometry
    
    Returns:
        Dictionary with results for each bond length
    """
    results = {
        'bond_lengths': bond_lengths,
        'energies': {},  # {bond_length: [energies per iteration]}
        'final_energies': [],
        'scf_energies': [],
        'method': method,
        'basis': basis_name
    }
    
    is_custom_basis = os.path.exists(basis_name) or basis_name.endswith('.gbs')
    basis_short = basis_name.split('/')[-1].replace('.gbs', '').replace('-', '')
    
    # Determine active space for H2
    if "PVDZ" in basis_name.upper():
        active_electrons, active_orbitals = 2, 10
    else:
        active_electrons, active_orbitals = 2, 2
    
    for r_idx, r in enumerate(bond_lengths):
        print(f"\n{'='*60}")
        print(f"Bond Length: {r:.4f} Angstrom ({r_idx+1}/{len(bond_lengths)})")
        print(f"{'='*60}")
        
        results['energies'][r] = []
        
        # File prefix for this geometry
        
        for iteration in range(n_iterations):
            print(f"\n--- Iteration {iteration + 1} of {n_iterations} ---")
            
            # File prefix for this iteration (includes iteration number to avoid overwriting)
            file_prefix = f"H2_{basis_short}_{method.upper()}_R{r:.3f}_iter{iteration}".replace('.', 'p')
            
            # Determine basis to use
            if iteration == 0:
                if is_custom_basis:
                    with open(basis_name, 'r') as f:
                        basis_content = f.read()
                    basis_label = f"INIT_BASIS_{r_idx}"
                    psi4.basis_helper(basis_content, name=basis_label, key='BASIS', set_option=True)
                else:
                    basis_label = basis_name
                print(f"  Basis: {basis_name} (initial)")
            else:
                # Look for the optimized basis file from previous iteration
                # The plugin may create uppercase or lowercase filenames
                prev_file_prefix = f"H2_{basis_short}_{method.upper()}_R{r:.3f}_iter{iteration-1}".replace('.', 'p')
                opt_basis_file_lower = f"{prev_file_prefix}.gbs"
                opt_basis_file_upper = f"{prev_file_prefix.upper()}.gbs"
                
                # Try lowercase first, then uppercase
                opt_basis_file = None
                if os.path.exists(opt_basis_file_lower):
                    opt_basis_file = opt_basis_file_lower
                elif os.path.exists(opt_basis_file_upper):
                    opt_basis_file = opt_basis_file_upper
                
                if opt_basis_file:
                    try:
                        with open(opt_basis_file, 'r') as f:
                            basis_content = f.read()
                        basis_label = f"OPT_{r_idx}_{iteration}"
                        psi4.basis_helper(basis_content, name=basis_label, key='BASIS', set_option=True)
                        print(f"  Basis: {opt_basis_file} (optimized from iteration {iteration})")
                    except Exception as e:
                        print(f"  Warning: Could not load optimized basis: {e}")
                        if is_custom_basis:
                            with open(basis_name, 'r') as f:
                                basis_content = f.read()
                            basis_label = f"INIT_BASIS_{r_idx}_{iteration}"
                            psi4.basis_helper(basis_content, name=basis_label, key='BASIS', set_option=True)
                        else:
                            basis_label = basis_name
                else:
                    if is_custom_basis:
                        with open(basis_name, 'r') as f:
                            basis_content = f.read()
                        basis_label = f"INIT_BASIS_{r_idx}_{iteration}"
                        psi4.basis_helper(basis_content, name=basis_label, key='BASIS', set_option=True)
                    else:
                        basis_label = basis_name
                    print(f"  Basis: {basis_name} (no optimized file from iteration {iteration})")
            
            try:
                mol = get_h2_molecule(r)
                e_scf, e_total = run_single_point(
                    mol, basis_label, method, 
                    active_electrons, active_orbitals,
                    generate_basis=True, basis_file_name=file_prefix
                )
                
                results['energies'][r].append(e_total)
                print(f"  SCF Energy: {e_scf:.8f} Hartree")
                print(f"  Total Energy: {e_total:.8f} Hartree")
                print(f"  ✓ Completed")
                
            except Exception as e:
                print(f"  ✗ Failed: {e}")
                import traceback
                traceback.print_exc()
                results['energies'][r].append(float('nan'))
        
        # Store final energy for this geometry
        if results['energies'][r]:
            results['final_energies'].append(results['energies'][r][-1])
        else:
            results['final_energies'].append(float('nan'))
    
    return results

def print_reaction_curve_table(results):
    """Print reaction curve results in a table."""
    print(f"\n{'='*80}")
    print(f"H2 REACTION CURVE RESULTS")
    print(f"Method: {results['method'].upper()}, Basis: {results['basis']}")
    print(f"{'='*80}")
    
    print(f"\n{'Bond Length (A)':<16}", end="")
    n_iter = max(len(results['energies'][r]) for r in results['bond_lengths'])
    for i in range(n_iter):
        print(f"{'Iter '+str(i+1):>14}", end="")
    print()
    
    for r in results['bond_lengths']:
        print(f"{r:<16.4f}", end="")
        for e in results['energies'][r]:
            if np.isnan(e):
                print(f"{'N/A':>14}", end="")
            else:
                print(f"{e:>14.8f}", end="")
        print()
    
    print(f"\nFinal Energies:")
    print(f"{'Bond Length (A)':<16}{'Energy (Ha)':>16}")
    for r, e in zip(results['bond_lengths'], results['final_energies']):
        if np.isnan(e):
            print(f"{r:<16.4f}{'N/A':>16}")
        else:
            print(f"{r:<16.4f}{e:>16.8f}")

def save_results(results, filename):
    """Save results to JSON file."""
    # Convert numpy types to Python types
    save_data = {
        'bond_lengths': list(results['bond_lengths']),
        'final_energies': [float(e) if not np.isnan(e) else None for e in results['final_energies']],
        'method': results['method'],
        'basis': results['basis'],
        'energies_per_iteration': {
            str(r): [float(e) if not np.isnan(e) else None for e in energies]
            for r, energies in results['energies'].items()
        }
    }
    
    with open(filename, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"Results saved to {filename}")

def main():
    """Main test function."""
    # Configuration
    bond_lengths = np.linspace(0.5, 3.0, 11)  # 0.5 to 3.0 Angstrom
    methods = ["mp2", "dmrg", "ccsd", "cisd"]  # Test with MP2, DMRG, CCSD and CISD
    n_iterations = 3  # Fewer iterations for speed
    
    # Basis set - use cc-pVDZ as requested
    basis_name = "cc-pVDZ"
    
    psi4.core.set_output_file('H2_reaction_curve_output.log', False)
    
    print(f"{'='*80}")
    print("H2 REACTION CURVE TEST WITH ITERATIVE BASIS UPDATING")
    print(f"Basis: {basis_name}")
    print(f"Bond lengths: {bond_lengths[0]:.2f} to {bond_lengths[-1]:.2f} Angstrom")
    print(f"Iterations per geometry: {n_iterations}")
    print(f"{'='*80}")
    
    all_results = {}
    
    for method in methods:
        print(f"\n{'='*80}")
        print(f"TESTING {method.upper()} METHOD")
        print(f"{'='*80}")
        
        try:
            results = run_h2_reaction_curve(
                basis_name, method, bond_lengths, n_iterations
            )
            all_results[method] = results
            
            print_reaction_curve_table(results)
            save_results(results, f"H2_reaction_curve_{method}_{basis_name.replace('-','')}.json")
            
        except Exception as e:
            print(f"Error testing {method}: {e}")
            import traceback
            traceback.print_exc()
    
    # Print comparison
    print(f"\n{'='*80}")
    print("COMPARISON OF METHODS")
    print(f"{'='*80}")
    
    print(f"\n{'Bond Length (A)':<16}", end="")
    for method in methods:
        print(f"{method.upper():>14}", end="")
    print()
    
    for i, r in enumerate(bond_lengths):
        print(f"{r:<16.4f}", end="")
        for method in methods:
            if method in all_results:
                e = all_results[method]['final_energies'][i]
                if np.isnan(e):
                    print(f"{'N/A':>14}", end="")
                else:
                    print(f"{e:>14.8f}", end="")
            else:
                print(f"{'N/A':>14}", end="")
        print()
    
    print(f"\n{'='*80}")
    print("TEST COMPLETE")
    print(f"{'='*80}")
    
    return all_results

if __name__ == "__main__":
    results = main()
