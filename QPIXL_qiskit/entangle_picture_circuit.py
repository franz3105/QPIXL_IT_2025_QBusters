import qiskit as qk


def entangle_picture_circuit(circ_list, nq_list):

    # Create a Quantum Circuit with 2 qubits
    ctrl_circ = qk.QuantumCircuit(1)
    circ_list = [qk.QuantumCircuit(nq) for nq in range(nq_list)]

    for i_circ, circ in enumerate(circ_list):
        # Create a Quantum Circuit with 2 qubits
        


    
    # Apply Hadamard gate to the first qubit
    circuit.h(0)

    # Apply CNOT gate to entangle the qubits
    circuit.cx(0, 1)

    # Draw the circuit
    circuit.draw('mpl')

    return circuit