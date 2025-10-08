import numpy as np
import matplotlib.pyplot as plt
from compute_work_isothermal import compute_work_isothermal
from compute_work_adiabatic import compute_work_adiabatic

# Compute both
Vf_iso, w_iso = compute_work_isothermal()
Vf_adi, w_adi = compute_work_adiabatic()

# Plot
plt.figure(figsize=(7, 5))
plt.plot(Vf_iso, w_iso, label="Isothermal", lw=2)
plt.plot(Vf_adi, w_adi, label="Adiabatic", lw=2)
plt.xlabel("Final Volume Vf (m³)")
plt.ylabel("Work Done on Gas (J)")
plt.title("Work Done in Isothermal vs Adiabatic Expansion")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("work_comparison.png", dpi=300)
plt.show()