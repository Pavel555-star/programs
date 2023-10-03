import numpy as np
from qiskit import QuantumCircuit, Aer, execute


backend = Aer.get_backend('qasm_simulator')
n_count = 8
a = 7

def c_amod15(a, exponent):
    m = QuantumCircuit(4)
    for iterations in range (exponent):
        m.swap(2,3)
        m.swap(1,2)
        m.swap(0,1)
        for q in range (4):
            m.x(q)
        m = m.to_gate()
        m.name = "%i ^ %i mod 15" %(a, exponent)
        c_m = m.control()
        return(c_m)
    

def qft_cross(n):
    d = QuantumCircuit(n)
    for qubit in range(n//2):
        d.swap(qubit, n-qubit-1)
    for i in range(n):
       for j in range(i):
           d.cp(-np.pi/float(2*(i-j)), j, i) # cu1 matrix ((1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,e^il))
       d.h(i)

    d.name = "QFT Dagger"
    return d

qc = QuantumCircuit(n_count + 4, n_count)
for q in range(n_count):   
    qc.h(q)

qc.x(3 + n_count)

for q in range(n_count):
    qc.append(c_amod15(a, 2**q), [q] + [i + n_count for i in range(4)])
qc.append(qft_cross(n_count), range(n_count))
qc.measure(range(n_count), range(n_count))
print(qc.draw("text"))
results = execute(qc, backend, shots = 32).result()
counts = results.get_counts()
print(counts)
