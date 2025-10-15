from scipy.optimize import minimize
import math
import numpy as np

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
def total_energy(coords, epsilon=0.01, sigma=3.4):
    # coords = [x2, y2, x3, y3]
    x2, y2, x3, y3 = coords
    # Atom 1 at (0,0,0), Atom 2 at (x2,y2,0), Atom 3 at (x3,y3,0)
    r12 = np.hypot(x2, y2)
    r13 = np.hypot(x3, y3)
    r23 = np.hypot(x3 - x2, y3 - y2)
    return lennard_jones(r12, epsilon, sigma) + lennard_jones(r13, epsilon, sigma) + lennard_jones(r23, epsilon, sigma)


#Initial guess for the distance = 4.0 Angstroms for atom 2 and (2.0, 3.5) for atom 3
initial_guess = [4.0, 0.0, 2.0, 3.5]

#We use the minimize function from scipy to find the minimum of the Lennard-Jones potential
result = minimize(lennard_jones, initial_guess, bounds=[(0, None), (None, None), (0, None), (None, None)])
x2, y2, x3, y3 = result.x

#Building a dictionary to represent the molecule
molecule = {
    'Ar1': np.array([0.0, 0.0, 0.0]),
    'Ar2': np.array([x2, y2, 0.0]),
    'Ar3': np.array([x3, y3, 0.0])
}

#Function to print the molecule in XYZ format
def print_xyz(molecule, comment="Argon cluster"):
    print(len(molecule))
    print(comment)
    for label, coords in molecule.items():
        #Use 'Ar' for argon, format to 6 decimal places
        print(f"Ar {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}")

#Results of the optimization
print(f"Minimum at r = {result.x[0]:.4f} Angstroms")
print(f"Minimum potential = {result.fun:.6f} eV")

#Now we can plot the Lennard-Jones potential to visualize it
import matplotlib.pyplot as plt

#Generating values for r from 3 to 6 Angstroms
r_values = np.linspace(3, 6, 100)
lj_values = [lennard_jones([r]) for r in r_values]
plt.plot(r_values, lj_values, label="Lennard-Jones Potential")

#Marking the equilibrium distance (minimum point) on the graph in red
plt.scatter(result.x[0], result.fun, color='red', label="Equilibrium Point")

#Formatting for the plot
plt.title("Lennard-Jones Potential")
plt.xlabel("Distance (Å)")
plt.ylabel("Potential Energy (eV)")
plt.legend()
plt.grid(True)
plt.show()

#To print bond lengths and bond angles, I'll be reusing some code from previous assignments

#Bond length function
def compute_bond_length(coords1, coords2):
    return math.sqrt((coords1[0] - coords2[0])**2 + (coords1[1] - coords2[1])**2 + (coords1[2] - coords2[2])**2)
#Input sorting function
def coord_sorter(input_str):
    coords = input_str.split()
    if len(coords) % 4 != 0:
        print("Error: Each atom should be listed with a unique ID and its x, y, z coordinates.")
        return None
    molecule = {}
    for x in range(0, len(coords), 4):
        label = coords[x]
        coord_list = np.array(list(map(float, coords[x+1:x+4])))
        molecule[label] = coord_list
    return molecule
#Function to calculate and print all bond lengths in the molecule
def calculate_all_bond_lengths(molecule):
    atom_labels = list(molecule.keys())
    n = len(atom_labels)
    for i in range(n):
        for j in range(i + 1, n):
            atom1 = atom_labels[i]
            atom2 = atom_labels[j]
            bond_length = compute_bond_length(molecule[atom1], molecule[atom2])
            print(f"Bond length between {atom1} and {atom2}: {bond_length:.2f} Angstroms")

#Computing bond angle code from previous assignment
#Bond angle function
def compute_bond_angle(coord1, coord2, coord3):
    vector1 = coord1 - coord2
    vector2 = coord3 - coord2
    dot_product = np.dot(vector1, vector2)
    norm1 = np.linalg.norm(vector1)
    norm2 = np.linalg.norm(vector2)
    if norm1 == 0 or norm2 == 0:
        return None
    cos_angle = dot_product / (norm1 * norm2)
    angle_rad = np.arccos(cos_angle)
    angle_deg = np.degrees(angle_rad)
    return angle_deg
#Function to calculate and print all bond angles in the molecule
def calculate_all_bond_angles(molecule):
    atom_labels = list(molecule.keys())
    n = len(atom_labels)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                atom1 = atom_labels[i]
                atom2 = atom_labels[j]
                atom3 = atom_labels[k]
                bond_angle = compute_bond_angle(molecule[atom1], molecule[atom2], molecule[atom3])
                if bond_angle is None:
                    print(f"Cannot calculate bond angle for {atom1}, {atom2}, {atom3}. (Vector is zero-length)")
                elif abs(bond_angle - 90) < 1e-2:
                    print("The bond angle is right.")
                elif abs(bond_angle - 180) < 1e-2:
                    print("The bond angle is linear.")
                elif abs(bond_angle) < 1e-2:
                    print("The bond angle is invalid.")
                elif bond_angle < 90:
                    print("The bond angle is acute.")
                elif bond_angle > 180:
                    print("The bond angle is too large.")
                else:
                    print("The bond angle is obtuse.")
                print(f"Bond angle between {atom1}, {atom2}, and {atom3}: {bond_angle:.2f} degrees")
             
                print()

if molecule:
    print()
    calculate_all_bond_lengths(molecule)
    print()
    calculate_all_bond_angles(molecule)
    print_xyz(molecule, comment="Argon cluster")