# This file is a lab aimed at solving the Schrodinger equation
# for arbitrary potentials, in 3D.
#
# Some methods to consider:
#
# - Density Functional theory
# - Finite Element Method
# - Diffusion Monte carlo
# - Hartree Fock
# - DVR
# - Spectral methods
# Kohn-Sham / plane-wave basis

# Article comparing plane waves to atomic orbitals:
# https://davidbowler.github.io/AtomisticSimulations/blog/basis-sets-in-dft
# Has some notable detractions about plane waves. Good info in general
# comparing the 2.

# https://chem.libretexts.org/Courses/University_of_California_Davis/UCD_Chem_107B%3A_Physical_Chemistry_for_Life_Scientists/Chapters/4%3A_Quantum_Theory/4.10%3A_The_Schr%C3%B6dinger_Wave_Equation_for_the_Hydrogen_Atom

# Good resource overall?: http://davidbowler.github.io/AtomisticSimulations/index.html


# Basis sets:
# - psi values (Is this called real-space?)
# - plane waves
# - gaussians
# - atomic orbitals

from dataclasses import dataclass;

import numpy as np
import matplotlib.pyplot as plt

from numpy import sqrt, pi, cos

A_0 = 1.
Z_H = 1.
K_C = 1.
Q_PROT = 1.
Q_ELEC = -1.
# M_ELEC = 5.45e-4
M_ELEC = 1.  # todo: Which?
ħ = 1

@dataclass
class Vec3:
    x: float
    y: float
    z: float

    def __sub__(self, other):
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)


def h_V(posit_nuc: Vec3, posit_sample: Vec3) -> float:
    """Analytic solution for n=1, s orbital"""
    diff = posit_sample - posit_nuc
    r = sqrt(diff.x**2 + diff.y**2 + diff.z**2)

    # todo: Why negative? Seems we need it.
    return -K_C * Q_PROT / r


def h_wf_100(posit_nuc: Vec3, posit_sample: Vec3) -> float:
    """Analytic solution for n=1, s orbital"""
    diff = posit_sample - posit_nuc
    r = sqrt(diff.x**2 + diff.y**2 + diff.z**2)

    ρ = Z_H * r / A_0
    return 1. / sqrt(pi) * (Z_H / A_0)**(3. / 2.) * np.exp(-ρ)
    # return 1. / sqrt(pi) * 1./ A_0**(3. / 2.) * np.exp(-ρ)


def h_wf_200(posit_nuc: Vec3, posit_sample: Vec3) -> float:
    diff = posit_sample - posit_nuc
    r = sqrt(diff.x**2 + diff.y**2 + diff.z**2)

    ρ = Z_H * r / A_0
    return 1. / sqrt(32.*pi) * (Z_H / A_0)**(3. / 2.) * (2. - ρ) * np.exp(-ρ/2.)


def h_wf_210(posit_nuc: Vec3, posit_sample: Vec3, theta: float) -> float:
    diff = posit_sample - posit_nuc
    r = sqrt(diff.x**2 + diff.y**2 + diff.z**2)

    ρ = Z_H * r / A_0
    return 1. / sqrt(32.*pi) * (Z_H / A_0)**(3. / 2.) * ρ * np.exp(-ρ/2.) * cos(theta)



def main():
    # Let's stick to 3D. And try the H WF as a test that it
    # works, to validate our method.

    posit_nuc = Vec3(0., 0., 0.)

    # Schrod eq for H: 
    # V for hydrogen: K_C * Q_PROT / r

    # psi(r)'' = (E - V(r)) * 2*m/ħ**2 * psi(r)
    # psi(r) = (E - V(R))^-1 * ħ**2/2m * psi(r)''

    E = -1/2 # n = 1

    N = 10_000

	# As a test iter, keep y and z = 0

    V_vals = np.zeros(N)
    psi = np.zeros(N)
    psi_pp_expected = np.zeros(N)
    psi_pp_measured = np.zeros(N)
    x_vals = np.linspace(-10, 10, N)

    # Used for calculating numerical psi''.
    # Smaller is more precise. Applies to dx, dy, and dz
    h = 0.00001  # aka dx?

    wf = h_wf_100

    for i, x in enumerate(x_vals):
        posit_sample = Vec3(x, 0., 0.)

        psi[i] = wf(posit_nuc, posit_sample)

        # todo: Experimenting with negative sign here. Why?
        V = h_V(posit_nuc, posit_sample)
        V_vals[i] = V

        psi_pp_expected[i] = (E - V) * -2. * M_ELEC / ħ**2 * psi[i]

        # Calculate psi'' based on a numerical derivative of psi
        # in 3D.

        # This crude approach uses a 3d-block style grid. Maybe works for now.
        # todo: Better, more spherical way of working this??
        # This is currently a plus-style grid.
        # todo: Also - reduce your approach to simplified operations by
        # writing it out, or ref a prev doc that probably exists that you
        # wrote it out on.
        # We're 

        # todo: You need to subtract and scale by an appropriate dt on your adding
        # todo and subtracting from x, y etc.

        # todo: second deriv is perhaps avg of the above? Validate this.
        # found in prev notes. Test on 1D, 2D, 3D.

        x_prev = Vec3(posit_sample.x - h, posit_sample.y, posit_sample.z)
        x_next = Vec3(posit_sample.x + h, posit_sample.y, posit_sample.z)
        y_prev = Vec3(posit_sample.x, posit_sample.y - h, posit_sample.z)
        y_next = Vec3(posit_sample.x, posit_sample.y + h, posit_sample.z)
        z_prev = Vec3(posit_sample.x, posit_sample.y, posit_sample.z - h)
        z_next = Vec3(posit_sample.x, posit_sample.y, posit_sample.z + h)

        psi_x_prev = wf(posit_nuc, x_prev)
        psi_x_next = wf(posit_nuc, x_next)
        psi_y_prev = wf(posit_nuc, y_prev)
        psi_y_next = wf(posit_nuc, y_next)
        psi_z_prev = wf(posit_nuc, z_prev)
        psi_z_next = wf(posit_nuc, z_next)

        # Could combine these, but this split-by-axis is easier to understand.

        # todo: Is this how you compensate for h??. We removed the 1/4 term and added h^2.
        
        # Initialize using the second derivative temr
        psi_pp_measured[i] = 0.
        psi_pp_measured[i] += (psi_x_prev + psi_x_next - 2. * psi[i])
        psi_pp_measured[i] += (psi_y_prev + psi_y_next - 2. * psi[i])
        psi_pp_measured[i] += (psi_z_prev + psi_z_next - 2. * psi[i])
        psi_pp_measured[i] /= h**2

		# psi_pp_measured[i] = 0.25 * psi_x_prev2 + psi_x_next2 + psi_y_prev2 + psi_y_next2 + \
		# 	psi_z_prev2 + psi_z_next2 - 6. * psi[i]

    scaler = 1  # to make the plots show together.

    plt.plot(x_vals, V_vals)
    plt.plot(x_vals, psi)
    plt.plot(x_vals, psi_pp_expected * scaler)  # Scaler for visualization
    plt.plot(x_vals, psi_pp_measured * scaler)

    plt.ylim(-0.5, 0.8)

    plt.show()



main()