import math
import numpy as np
from scipy.optimize import minimize

#First, we define the Lennard-Jones potential function
def lennard_jones(r, epsilon=0.01, sigma=3.4):
    """
    Calculates the Lennard-Jones potential between two particles. (Specifically Argon)

    Parameters:
    r (float): Distance between the two particles.
    epsilon (float): The minimum energy value. (depth of the potential well)
    sigma (float): Distance at which the potential is zero.

    Returns:
    float: The Lennard-Jones potential at distance r.
    """
    if isinstance(r, (list, np.ndarray)):
        r = r[0]
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

#We need a function for the total energy of the system
def total_energy(coords, r12=1.0, epsilon=0.01, sigma=3.4):
    #Atom 1: (0, 0)
    #Atom 2: (r12, 0)
    #Atom 3: (x3, y3) from coords
    x3, y3 = coords
    
    #Distances
    r13 = np.sqrt(x3**2 + y3**2)
    r23 = np.sqrt((x3 - r12)**2 + y3**2)
    
    #Sum of LJ potentials between each pair
    V_total = (
        lennard_jones(r12, epsilon, sigma) +
        lennard_jones(r13, epsilon, sigma) +
        lennard_jones(r23, epsilon, sigma)
    )

    return V_total

#Initial guess for (x3, y3)
initial_guess = [0.5, 1.0]

#Minimize total energy
result = minimize(total_energy, initial_guess)

#Extracting the optimal coordinates
x3, y3 = result.x
r12 = 1.0
r13 = np.sqrt(x3**2 + y3**2)
r23 = np.sqrt((x3 - r12)**2 + y3**2)

print("Optimal position of third atom (x3, y3):", result.x)
print("Minimum total potential energy:", result.fun)

#Print optimal distances
print()
print(f"Optimal distances (in same units as input):")
print(f"r12 (Atom 1-2): {r12:.6f}")
print(f"r13 (Atom 1-3): {r13:.6f}")
print(f"r23 (Atom 2-3): {r23:.6f}")

#Compute angles using Law of Cosines
#Generalize angle calculation so any input order of sides works
def angle(a, b, c):
    #Angle opposite side a, given sides a, b, c
    return math.degrees(math.acos((b**2 + c**2 - a**2) / (2 * b * c)))

#Angles at each atom
angle1 = angle(r23, r12, r13)  #at Atom 1
angle2 = angle(r13, r12, r23)  #at Atom 2
angle3 = angle(r12, r13, r23)  #at Atom 3

#Printing results
print()
print(f"Optimal angles (in degrees):")
print(f"Angle at Atom 1: {angle1:.2f}")
print(f"Angle at Atom 2: {angle2:.2f}")
print(f"Angle at Atom 3: {angle3:.2f}")

#Comment on geometry of angles
#1e-2 gives enough of a buffer to cover most floating point issues
if abs(r12 - r13) < 1e-2 and abs(r12 - r23) < 1e-2:
    print("The atoms form an approximately equilateral triangle.")
else:
    print("The atoms do not form an equilateral triangle.")