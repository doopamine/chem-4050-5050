import numpy as np
import matplotlib.pyplot as plt

# We define a constant seed for reproducibility
np.random.seed(42)
# We define the half-length of the cube in which we will sample points
L = 20

# Section 1: Define the Hydrogen 2p Orbital Function
# We define the wavefunction for any point (x,y,z) in space
# Important to note that the Bohr radius a_0 will be set to 1 for simplicity sake
def psi_2p_z(x,y,z):
    psi = (1/(4*np.sqrt(2*np.pi))) * np.sqrt(x**2+y**2+z**2) * z/np.sqrt(x**2+y**2+z**2) * np.exp(-np.sqrt(x**2+y**2+z**2)/2)
    return psi

# Section 2: Compute the Overlap Integral Using Random Sampling
# Compute overlap integral S(R) by random uniform sampling
# We will sample N points uniformly in the cube of side length 2L centered at the origin
def random_overlap(R, N, L):
    # We sample random points in cube [-L, L]
    x = np.random.uniform(-L, L, N)
    y = np.random.uniform(-L, L, N)
    z = np.random.uniform(-L, L, N)

    psi1 = psi_2p_z(x, y, z + R / 2) # Wavefunction centered at +R/2
    psi2 = psi_2p_z(x, y, z - R / 2) # Wavefunction centered at -R/2

    integrand = psi1 * psi2  # The integrand of the overlap integral
    avg = np.mean(integrand) # Average value of the integrand
    V = (2 * L) ** 3 # Volume of the cube
    return V * avg

# Section 3: Improve Efficiency Using Importance Sampling
# We will sample points according to a distribution similar to |ψ|^2
def importance_sampling(N):
    sigma = 1 # Sigma = a0, which is 1 in our units
    x = np.random.normal(0, sigma, N)
    y = np.random.normal(0, sigma, N)
    z = np.random.normal(0, sigma, N)
    return x, y, z

# We can compute the probability density function for Gaussian sampling
def prob_density(x, y, z, sigma=1):
    coefficients = 1 / ((2 * np.pi * sigma**2) ** 1.5)
    return coefficients * np.exp(-(x**2 + y**2 + z**2) / (2 * sigma**2))

# We can compute the overlap integral S(R) by using importance sampling
def importance_overlap(R, N):
    x, y, z = importance_sampling(N) # Sample points according to the Gaussian distribution
    psi1 = psi_2p_z(x, y, z + R / 2) 
    psi2 = psi_2p_z(x, y, z - R / 2)
    g_values = prob_density(x, y, z) # Evaluate the probability density at the sampled points
    integrand = psi1 * psi2 / g_values
    return np.mean(integrand)

# Section 4: Plot the Overlap Integral as a Function of Separation Distance
# We can plot S(R) for R from 0 to 10
def plot_overlap(R=2.0, L=20.0):
    N_values = [10**k for k in range(2, 8)]  # up to 1e7 samples (I went to 1e8 but it took too long)
    s_random, s_importance = [], [] 
    for N in N_values:
        s_random.append(random_overlap(R, N, L)) # Random sampling
        s_importance.append(importance_overlap(R, N)) # Importance sampling
    # Graphing the results
    plt.figure()
    plt.loglog(N_values, np.abs(s_random), "o-", label="Random sampling")
    plt.loglog(N_values, np.abs(s_importance), "s-", label="Importance sampling")
    plt.xlabel("Number of samples N")
    plt.ylabel("|S(R)|")
    plt.title(f"Convergence of Overlap Integral S(R={R})")
    plt.legend()
    plt.grid(True, which="both")
    plt.savefig("convergence_overlap.png", dpi=300)
    plt.show()
    return N_values, s_random, s_importance

def overlap_vs_R():
    R_values = np.arange(0.5, 20.5, 0.5) # R from 0.5 to 20 in steps of 0.5
    N = 10**6 # This number of samples takes a bit longer, but provides a better graph 
    S_R = [] 
    for R in R_values: 
        S_R.append(importance_overlap(R, N)) # Importance sampling
    # Graphing the results again but against R
    plt.figure()
    plt.plot(R_values, S_R, "o-")
    plt.xlabel("Separation distance R (a0)")
    plt.ylabel("Overlap integral S(R)")
    plt.title("Overlap Integral vs. Separation Distance")
    plt.grid(True)
    plt.savefig("overlap_vs_R.png", dpi=300)
    plt.show()
    return R_values, S_R

if __name__ == "__main__":
    # Running the convergence test at R=2
    plot_overlap(R=2.0, L=L)
    # Compute and plot S(R)
    overlap_vs_R()

# Using a timing function:
# 1e7 samples = 28.75 seconds
# 1e6 samples = 20.36 seconds
# 1e5 samples = 19.46 seconds

# Seperately changing N:
# N = 10**6 = 16.21 seconds
# Changing N from 10**6 to 10**7 increases time by ~12 seconds, but the graph looks significantly better.
# N = 10**8 took over a minute and didn't provide a single output, so I reverted it back to 10**7. 