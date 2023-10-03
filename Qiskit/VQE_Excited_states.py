from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper, QubitConverter
from qiskit.algorithms.optimizers import SLSQP
from qiskit.primitives import Estimator
from qiskit_nature.second_q.algorithms import GroundStateEigensolver, QEOM, VQEUCCFactory
from qiskit_nature.second_q.circuit.library import UCCSD

driver = PySCFDriver(
    atom="H 0 0 0; H 0 0 0.735",
    basis="sto3g",
    charge=0,
    spin=0,
    unit=DistanceUnit.ANGSTROM,
)

es_problem = driver.run()
converter = QubitConverter(JordanWignerMapper())

estimator = Estimator()
# This first part sets the ground state solver
# see more about this part in the ground state calculation tutorial
solver = VQEUCCFactory(estimator, UCCSD(), SLSQP())
gse = GroundStateEigensolver(converter, solver)
# The qEOM algorithm is simply instantiated with the chosen ground state solver and Estimator primitive
qeom_excited_states_solver = QEOM(gse, estimator, "sd")

qeom_results = qeom_excited_states_solver.solve(es_problem) # ElectronicStructureResult
# items hartree_fock_energy, nuclear_repulsion_energy, computed_energies, electronic_energies, total_energies, ...
print(qeom_results)
