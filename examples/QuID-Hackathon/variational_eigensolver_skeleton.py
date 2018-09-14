"""
Implementation of a variational eigensolver using ProjectQ.

The example shown here is from the paper "Scalable Quantum Simulation of 
Molecular Energies" by P.J.J. O'Malley et al. arXiv:1512.06860

We use the variational quantum eigensolver to find the ground state energy
of H2 for one specific bond length of 0.75 Angstrom.

Eq. 2 of the paper shows the functional to minimize and Eq. 3 shows the
coupled cluster ansatz for the trial wavefunction (using the unitary coupled
cluster approach). The Hamiltonian is given in Eq. 1. The coefficients can be
found in Table 1. Note that both the ansatz and the Hamiltonian can be
calculated using FermiLib which is a library for Quantum Chemistry on top of
ProjectQ.
"""

import projectq
from projectq.ops import All, Measure, QubitOperator, TimeEvolution, X

import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar


def get_hamiltonian(bond_distance):
    """ 
    Returns the Hamiltonian of the system for a secific bond distance.

    Uses the data in Table 1 of the paper. One could use FermiLib to
    calculate and the Hamiltonian for an arbitrary molecule and bond distance.

    Args:
        bond_radius(float): Only bond distances of paper allowed, i.e.
                            0.2 <= bond_radius <= 2.85 in steps of 0.05
    Returns:
        QubitOperator which represents the Hamiltonian
    """
    with open("hamiltonian.txt", 'r') as file:
        found = False
        for line in file:
            r, g0, g1, g2, g3, g4, g5 = [float(x) for x in line.split()]
            if r == bond_distance:
                found = True
                break
        if not found:
            raise Exception("Hamiltonian for bond_distance ", bond_distance,
                            " not found in data file.")
    hamiltonian = g0*QubitOperator(()) # = identity
    hamiltonian += g1*QubitOperator("Z0")
    hamiltonian += g2*QubitOperator("Z1")
    hamiltonian += g3*QubitOperator("Z0 Z1")
    hamiltonian += g4*QubitOperator("X0 X1")
    hamiltonian += g5*QubitOperator("Y0 Y1")
    return hamiltonian


def energy(theta, bond_distance):
    hamiltonian = get_hamiltonian(bond_distance)
    # create a ProjectQ compiler with a simulator as a backend
    
    # TODO: implement algorithm to calculate energy

    return energy


if __name__ == '__main__':
    bond_distances = [x/100. for x in range(20, 285, 5)]
    minimal_energies = []
    for bond_distance in bond_distances:
        minimum = minimize_scalar(lambda theta: energy(theta, bond_distance))
        minimal_energies.append(minimum.fun)

    # print result
    plt.xlabel("Parameter $\\theta$")
    plt.ylabel("Expectation value $E(\\theta)$")
    plt.title("Expectation value $E = \left<\psi(\\theta) | H | "
             +"\psi(\\theta)\\right>$ for different $\\theta$")
    plt.plot(bond_distances, minimal_energies, "b.-",
             label="Ground-state energy")
    plt.show()
