from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.algorithms import GroundStateEigensolver, NumPyMinimumEigensolverFactory
from qiskit_nature.second_q.mappers import JordanWignerMapper, QubitConverter

driver = PySCFDriver(
    atom="H 0 0 0; H 0 0 0.735",
    basis="sto3g",
    charge=0,
    spin=0,
    unit=DistanceUnit.ANGSTROM,
)

problem = driver.run()
hamiltonian = problem.hamiltonian
coefficients = hamiltonian.electronic_integrals
second_q_op = hamiltonian.second_q_op()

solver = GroundStateEigensolver(
    QubitConverter(JordanWignerMapper()),
    NumPyMinimumEigensolverFactory(),
)

result = solver.solve(problem) # ElectronicStructureResult
print(result) # items hartree_fock_energy, nuclear_repulsion_energy, computed_energies, electronic_energies, total_energies, ...
