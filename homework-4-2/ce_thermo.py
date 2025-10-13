import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.constants import k as kB, eV

def partition_function(energies, degeneracies, T):
    # Computes partition function Z for any given values for energies and degeneracies
    beta = 1 / (kB * T) # in 1/J
    return np.sum(degeneracies * np.exp(-np.array(energies) * eV * beta))

def thermodynamic_properties(energies, degeneracies, T):
    # Solves variables U (J), F (J), and S (J/K) for given values
    beta = 1 / (kB * T) # in 1/J
    Z = partition_function(energies, degeneracies, T) 
    E_exp = np.sum(degeneracies * np.array(energies) * eV * np.exp(-np.array(energies) * eV * beta)) 
    U = E_exp / Z 
    F = -kB * T * np.log(Z) 
    S = (U - F) / T 
    return U, F, S

def compute_all_cases(T_range=np.linspace(300, 2000, 100)):
    # Compute Ce3+ for all 3 cases
    results = []

    # Case definitions were given 

    # Case 1: Isolated Ce3+
    energies1 = [0]
    degeneracies1 = [14]

    # Case 2: With SOC
    energies2 = [0, 0.28]
    degeneracies2 = [6, 8]

    # Case 3: With SOC + CFS
    energies3 = [0.00, 0.07, 0.12, 0.13, 0.14]
    degeneracies3 = [4, 2, 2, 4, 2]

    # Prepare data storage
    data = {"T": [], "U1": [], "U2": [], "U3": [], "F1": [], "F2": [], "F3": [], "S1": [], "S2": [], "S3": []}

    # Calculate properties for each temperature and append data
    for T in T_range:
        for i, (E, g, label) in enumerate(zip(
            [energies1, energies2, energies3],
            [degeneracies1, degeneracies2, degeneracies3],
            ["1", "2", "3"])):
            U, F, S = thermodynamic_properties(E, g, T)
            data[f"U{label}"].append(U)
            data[f"F{label}"].append(F)
            data[f"S{label}"].append(S)
        data["T"].append(T)

    # Save as a CSV
    with open("ce_thermo.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["T (K)", "U1", "U2", "U3", "F1", "F2", "F3", "S1", "S2", "S3"])
        for i in range(len(T_range)):
            writer.writerow([
                data["T"][i],
                data["U1"][i], data["U2"][i], data["U3"][i],
                data["F1"][i], data["F2"][i], data["F3"][i],
                data["S1"][i], data["S2"][i], data["S3"][i]
            ])

# Plot 1: Internal Energy vs Temperature
    plt.figure(figsize=(8, 5))
    plt.plot(T_range, np.array(data["U1"]) / eV, label="U (Isolated)")
    plt.plot(T_range, np.array(data["U2"]) / eV, label="U (SOC)")
    plt.plot(T_range, np.array(data["U3"]) / eV, label="U (SOC+CFS)")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Internal Energy (eV)")
    plt.title("Ce³⁺ Internal Energy vs Temperature")
    plt.legend()
    plt.tight_layout()
    plt.savefig("ce_internal_energy.png", dpi=300)

# Plot 2: Entropy vs Temperature
    plt.figure(figsize=(8, 5))
    plt.plot(T_range, np.array(data["S1"]), label="S (Isolated)")
    plt.plot(T_range, np.array(data["S2"]), label="S (SOC)")
    plt.plot(T_range, np.array(data["S3"]), label="S (SOC+CFS)")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Entropy (J/K)")
    plt.title("Ce³⁺ Entropy vs Temperature")
    plt.legend()
    plt.tight_layout()
    plt.savefig("ce_entropy.png", dpi=300)

# Plot 3: Free Energy (F2 vs Temperature)
    plt.figure(figsize=(8, 5))
    plt.plot(T_range, np.array(data["F2"]) / eV, lw=2, color="blue")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Free Energy F₂ (eV)")
    plt.title("Free Energy vs Temperature — Ce³⁺ with Spin-Orbit Coupling")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("ce_free_energy_soc.png", dpi=300)

    plt.show()

if __name__ == "__main__":
    compute_all_cases()