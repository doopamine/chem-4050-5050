import numpy as np
from scipy.integrate import trapezoid
import csv

def compute_work_adiabatic(Vi=0.1, Vf_max=0.3, n=1.0, R=8.314, T=300.0, gamma=1.4, n_points=100):
    # Constant from initial condition: P_i * V_i^y = constant
    # Using known values for n, R, T to simplify P_i = nRT/V_i directly.
    Pi = n * R * T / Vi
    constant = Pi * Vi ** gamma

    Vf_values = np.linspace(Vi, Vf_max, n_points)
    works = []

    for Vf in Vf_values:
        V = np.linspace(Vi, Vf, 500)
        P = constant / (V ** gamma)
        w = -trapezoid(P, V)
        works.append(w)

    # Save results and export to CSV file
    with open("work_adiabatic.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Vf (m^3)", "Work (J)"])
        for v, w in zip(Vf_values, works):
            writer.writerow([v, w])

    return Vf_values, works

# Main execution
if __name__ == "__main__":
    compute_work_adiabatic()