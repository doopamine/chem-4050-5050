import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid

#Constants
kB = 1.380649e-23  #Boltzmann constant in J/K
NA = 6.02214076e23  #Avogadro's number in 1/mol
eV_to_J = 1.602176634e-19
angstrom_to_m = 1e-10

#Changeable parameters
sigma = 3.4 * angstrom_to_m       #3.4 Å -> meters
epsilon = 0.01 * eV_to_J          #0.01 eV -> Joules
lambda_factor = 1.5

#Integration grid
r_min = 1e-3
r_max = 5 * sigma
r_grid = np.linspace(r_min, r_max, 1000)

#Definitions
def hard_sphere(r, sigma):
    if isinstance(r, (list, np.ndarray)):
        r = r[0]
    return 1000 if r < sigma else 0.0

def square_well(r, epsilon, sigma, lambda_factor):
    if isinstance(r, (list, np.ndarray)):
        r = r[0]
    if r < sigma:
        return 1000
    elif sigma <= r < lambda_factor * sigma:
        return -epsilon
    else:
        return 0.0

def lennard_jones(r, epsilon, sigma):
    if isinstance(r, (list, np.ndarray)):
        r = r[0]
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

#B2V integrand
def integrand(r, T, potential_func):
    U = np.vectorize(potential_func)(r)
    return r**2 * (np.exp(-U / (kB * T)) - 1)

#Temperature range
T_range = np.arange(100, 801, 50)  # From 100 K to 800 K

#Store results
B2V_hs_list = []
B2V_sw_list = []
B2V_lj_list = []

#Loop over temperatures
for T in T_range:
    B2V_hs = trapezoid(integrand(r_grid, T, lambda r: hard_sphere(r, sigma)), r_grid)
    B2V_sw = trapezoid(integrand(r_grid, T, lambda r: square_well(r, epsilon, sigma, lambda_factor)), r_grid)
    B2V_lj = trapezoid(integrand(r_grid, T, lambda r: lennard_jones(r, epsilon, sigma)), r_grid)

    #Convert to molar units
    factor = 2 * np.pi * NA
    B2V_hs_list.append(B2V_hs * factor)
    B2V_sw_list.append(B2V_sw * factor)
    B2V_lj_list.append(B2V_lj * factor)

#Plotting
plt.figure(figsize=(10, 6))
plt.plot(T_range, B2V_hs_list, label='Hard Sphere', marker='o')
plt.plot(T_range, B2V_sw_list, label='Square Well', marker='s')
plt.plot(T_range, B2V_lj_list, label='Lennard-Jones', marker='^')
plt.axhline(0, color='gray', linestyle='--', label='B2V = 0')

plt.xlabel('Temperature (K)')
plt.ylabel(r'$B_2^V$ (m$^3$/mol)')
plt.title(r'Second Virial Coefficient $B_2^V$ vs Temperature')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()