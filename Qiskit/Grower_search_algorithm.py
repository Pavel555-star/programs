from qiskit import *
import matplotlib.pyplot as plt
import numpy as np
# define the oracle circuit
oracle = QuantumCircuit(2, name='oracle')
oracle.cz(0,1) # -> cz-matrix ((1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1))
oracle.to_gate()

backend = Aer.get_backend('statevector_simulator')



reflect = QuantumCircuit(2, name = 'reflection')
reflect.h([0,1])  # hadamard matrix - ((1,1), (1,-1))
reflect.z([0,1])  # z-matrix ((0,1), (1,0))
reflect.cz(0,1)   # cz-matrix ((1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1))
reflect.h([0,1])  # hadamard matrix - ((1,1), (1,-1))
reflect.to_gate() 

backend = Aer.get_backend('qasm_simulator')

grover_circ = QuantumCircuit(2, 2)
grover_circ.h([0,1]) # -> hadamard matrix - ((1,1), (1,-1))
grover_circ.append(oracle, [0,1])
grover_circ.append(reflect, [0,1])
grover_circ.measure([0,1],[0,1])


job = execute(grover_circ,backend, shots = 1)
result = job.result()
print(grover_circ.draw())
print(result.get_counts())
