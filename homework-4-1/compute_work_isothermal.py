import numpy as np
from scipy.integrate import trapezoid
import csv

# The entirety of this code is defining and calling the function for computing work done during an isothermal process
def compute_work_isothermal(Vi=0.1, Vf_max=0.3, n=1.0, R=8.314, T=300.0, n_points=100):
    
    # Generate volume values and compute work for each final volume
    Vf_values = np.linspace(Vi, Vf_max, n_points)
    works = []

    for Vf in Vf_values:
        V = np.linspace(Vi, Vf, 500) # Volume array for integration
        P = n * R * T / V # Pressure array for integration
        w = -trapezoid(P, V) # Work done during the process
        works.append(w) # Store the work value

    #Save results and export to a CSV file
    with open("work_isothermal.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Vf (m^3)", "Work (J)"])
        for v, w in zip(Vf_values, works):
            writer.writerow([v, w])

    return Vf_values, works

# Running the code for isothermal
if __name__ == "__main__":
    compute_work_isothermal()