# Example of implementing iterative phase estimation
# QuID Hackathon
import math

import projectq
from projectq.meta import Control
from projectq.ops import All, H, Measure, R, X

if __name__ == '__main__':
    # Let's choose an easy unitary U = R
    # R is the phase shift gate, see the docs and implementation
    # https://projectq.readthedocs.io/en/latest/projectq.ops.html#projectq.ops.R
    #
    # R(angle) = [[1, 0], [0, cmath.exp(1j * angle)]]
    #
    # We choose angle = 2pi*a where a = 0.101 in binary notation, then:
    #
    # R(2pi*a) |0> = 1 * |0>
    # R(2pi*a) |1> = exp(2pi*j*a) |1>
    #
    # Our goal is to determine `a` given U
    a1 = 1
    a2 = 0
    a3 = 1
    a = a1 * 1./2 + a2 * 1./4 + a3 * 1/8.
    U = R(2 * math.pi * a)

    eng = projectq.MainEngine()
    qubit = eng.allocate_qubit()
    X | qubit  # prepares qubit in eigenstate |1>

    # TODO implement phase estimation using one ancilla qubit:
    ancilla_qb = eng.allocate_qubit()



    # Before qubits go out of scope, measure them
    All(Measure) | qubit + ancilla_qb
