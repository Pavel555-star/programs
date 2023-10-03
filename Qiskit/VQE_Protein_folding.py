from qiskit_research.protein_folding.interactions.random_interaction import (RandomInteraction,)
from qiskit_research.protein_folding.interactions.miyazawa_jernigan_interaction import (MiyazawaJerniganInteraction,)
from qiskit_research.protein_folding.peptide.peptide import Peptide
from qiskit_research.protein_folding.protein_folding_problem import (ProteinFoldingProblem,)
from qiskit_research.protein_folding.penalty_parameters import PenaltyParameters
from qiskit.utils import algorithm_globals, QuantumInstance
from qiskit.circuit.library import RealAmplitudes
from qiskit.algorithms.optimizers import COBYLA
from qiskit.algorithms import NumPyMinimumEigensolver, VQE
from qiskit.opflow import PauliExpectation, CVaRExpectation
from qiskit import execute, Aer
import matplotlib.pyplot as plt

algorithm_globals.random_seed = 23
random_interaction = RandomInteraction()
mj_interaction = MiyazawaJerniganInteraction()
penalty_back = 10
penalty_chiral = 10
penalty_1 = 10
penalty_terms = PenaltyParameters(penalty_chiral, penalty_back, penalty_1)
peptide = Peptide("APRLR", ["", "", "F", "Y", ""])
protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)
qubit_op = protein_folding_problem.qubit_op()
print(qubit_op)
# set classical optimizerвведите аминокислотную последовательность
optimizer = COBYLA(maxiter=50)
# set variational ansatz
ansatz = RealAmplitudes(reps=1)
# set the backend
backend_name = "aer_simulator"
backend = QuantumInstance(
    Aer.get_backend(backend_name),
    shots=8192,
    seed_transpiler=algorithm_globals.random_seed,
    seed_simulator=algorithm_globals.random_seed,
)
counts = []
values = []

def store_intermediate_result(eval_count, parameters, mean, std):
    counts.append(eval_count)
    values.append(mean)

# initialize CVaR_alpha objective with alpha = 0.1
cvar_exp = CVaRExpectation(0.1, PauliExpectation())
# initialize VQE using CVaR
vqe = VQE(
    expectation=cvar_exp,
    optimizer=optimizer,
    ansatz=ansatz,
    quantum_instance=backend,
    callback=store_intermediate_result,
)

raw_result = vqe.compute_minimum_eigenvalue(qubit_op)
print(raw_result)

result = protein_folding_problem.interpret(raw_result=raw_result)
print("The bitstring representing the shape of the protein during optimization is: ",result.turn_sequence,)
print("The expanded expression is:", result.get_result_binary_vector())
print(f"The folded protein's main sequence of turns is: {result.protein_shape_decoder.main_turns}")
print(f"and the side turn sequences are: {result.protein_shape_decoder.side_turns}")

fig = plt.figure()
plt.plot(counts, values)
plt.ylabel("Conformation Energy")
plt.xlabel("VQE Iterations")
fig.add_axes([0.44, 0.51, 0.44, 0.32])
plt.plot(counts[40:], values[40:])
plt.ylabel("Conformation Energy")
plt.xlabel("VQE Iterations")
plt.show()

fig_2 = result.get_figure(title="Protein Structure", ticks=False, grid=True)
fig_2.get_axes()[0].view_init(10, 60)
