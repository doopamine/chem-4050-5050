"""
This code takes the input of a CSV file containing boiling points (TB),
enthalpies of vaporization (Hv), and substance class labels. It performs an
ordinary least squares (OLS) regression to find the relationship between TB and Hv,
calculates confidence intervals for the slope and intercept, and plots the data.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats

# OLS functions taken from Lecture 7
def ols_slope(x, y):
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    numerator = np.sum((x - x_mean) * (y - y_mean))
    denominator = np.sum((x - x_mean) ** 2)
    return numerator / denominator

def ols_intercept(x, y):
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    slope = ols_slope(x, y)
    return y_mean - slope * x_mean

def ols(x, y):
    slope = ols_slope(x, y)
    intercept = ols_intercept(x, y)
    return slope, intercept

# Loading data from CSV
# Had some issues with relative paths, so using absolute path
script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, "trouton.csv")
data = pd.read_csv(file_path)

# Extracting data from columns in the CSV
TB = data["T_B (K)"].values          # Boiling points in K
Hv = data["H_v (kcal/mol)"].values          # Enthalpies of vaporization in kJ/mol
Class = data["Class"].values    # Substance class labels (text)

# Perform OLS regression on the data from the CSV
a, b = ols(TB, Hv)              # slope and intercept
Hv_pred = a * TB + b            # predicted values

# 95% confidence intervals
n = len(TB) # number of data points
dof = n - 2 # degrees of freedom
residuals = Hv - Hv_pred 
s_err = np.sqrt(np.sum(residuals**2) / dof) # standard error of the estimate

t_val = stats.t.ppf(0.975, dof) # t-value for 95% CI
x_mean = np.mean(TB) 
Sxx = np.sum((TB - x_mean)**2) # sum of squares of x deviations

s_a = s_err / np.sqrt(Sxx) # standard error of the slope
s_b = s_err * np.sqrt(1/n + x_mean**2 / Sxx) # standard error of the intercept

a_CI = (a - t_val*s_a, a + t_val*s_a) # confidence interval for slope
b_CI = (b - t_val*s_b, b + t_val*s_b) # confidence interval for intercept

# Convert slope to J/mol·K
a_JmolK = a * 1000
s_a_JmolK = s_a * 1000

# Plot
plt.figure(figsize=(8,6))

# Encode text labels into numbers for coloring
class_labels, class_ids = np.unique(Class, return_inverse=True)

scatter = plt.scatter(TB, Hv, c=class_ids, cmap="viridis", label="Data")
plt.plot(TB, Hv_pred, color="red", label="Fit")

eq_text = (
    f"Hv = {a:.3f}·TB + {b:.2f}\n"
    f"Slope (ΔSv) = {a_JmolK:.1f} ± {s_a_JmolK:.1f} J/mol·K\n"
    f"Intercept = {b:.2f} ± {s_b:.2f} kJ/mol"
)

plt.text(0.05, 0.9, eq_text, transform=plt.gca().transAxes, fontsize=10,
         bbox=dict(facecolor="white", alpha=0.7))

plt.xlabel("Boiling Point TB [K]")
plt.ylabel("Enthalpy of Vaporization Hv [kJ/mol]")
plt.title("Trouton’s Rule")

# Add legend with original labels
handles, _ = scatter.legend_elements()
plt.legend(handles, class_labels, title="Class")

# Save plot
os.makedirs("homework-3-1", exist_ok=True)
plt.savefig("homework-3-1/troutons_rule.png", dpi=300)
plt.show()
plt.close()
