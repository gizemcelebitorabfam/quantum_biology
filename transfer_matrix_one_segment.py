import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# Define the potential function coefficients (should be for x in Angstroms)
p = [-0.0032480562835404563, 0.13387641022560062, -1.3976143856453171,
     5.300165914337288, -5.501800261100386, -4.892603997626288, 4.046960362958583]

# Constants for unit conversion
hbar = 1.0545718e-34  # J*s
m = 6.492e-26  # kg (mass of potassium ion)
eV_to_J = 1.60218e-19  # conversion factor from eV to Joules

# Particle energy (example)
E = 0.0178  # eV
E = E * eV_to_J  # convert to Joules

# Discretize the range (in Angstroms, so we use x in Angstroms)
x_start = 2.5  # Angstrom
x_end = 5  # Angstrom
N = 5000
x = np.linspace(x_start, x_end, N)  # x in Angstroms

# Convert x to meters for calculations later
x_meters = x * 1e-10

# Compute the segment size dx (in meters)
dx = (x_meters[-1] - x_meters[0]) / (N - 1)

# Define the potential in eV (before converting to Joules)
V = 0.0103643*(p[0]*x**6 + p[1]*x**5 + p[2]*x**4 + p[3]*x**3 + p[4]*x**2 + p[5]*x + p[6])

# Apply Gaussian smoothing to the potential (in eV)
V_smoothed = gaussian_filter1d(V, sigma=2)  # Apply some smoothing

# Plot the potential in eV (as a sixth-order polynomial)
plt.plot(x, V, label="Original Potential (eV)")
plt.plot(x, V_smoothed, label="Smoothed Potential (eV)")
plt.title("Potential V(x) in eV")
plt.xlabel("x (Angstroms)")
plt.ylabel("V (eV)")
plt.legend()
plt.show()

# Convert the potential to Joules for later calculations
V = V * eV_to_J  # convert to Joules
V_smoothed = V_smoothed * eV_to_J  # also convert the smoothed potential

# Now use V_smoothed for the transfer matrix calculations in Joules...
# Rest of the code follows as before

# Initialize the total transfer matrix
M_total = np.identity(2, dtype=complex)

def transfer_matrix_segment(V_segment, dx):
    if E == V_segment:
        print("E equals V_segment; returning approximate identity matrix.")
        return np.array([[1, dx], [0, 1]], dtype=complex)

    if E > V_segment:
        k2 = np.sqrt(2 * m * (E - V_segment)) / hbar
        print(f"k2 (real): {k2}")

        if k2 * dx > 20:
            exp_k2_dx = 0.5 * np.exp(k2 * dx)
            exp_neg_k2_dx = 0.5 / exp_k2_dx
            print(f"Exp terms (approximated): exp_k2_dx = {exp_k2_dx}, exp_neg_k2_dx = {exp_neg_k2_dx}")
        else:
            exp_k2_dx = np.exp(k2 * dx)
            exp_neg_k2_dx = np.exp(-k2 * dx)
            print(f"Exp terms: exp_k2_dx = {exp_k2_dx}, exp_neg_k2_dx = {exp_neg_k2_dx}")

        sinh_k2_dx = 0.5 * (exp_k2_dx - exp_neg_k2_dx)
        cosh_k2_dx = 0.5 * (exp_k2_dx + exp_neg_k2_dx)
        print(f"sinh(k2 * dx): {sinh_k2_dx}, cosh(k2 * dx): {cosh_k2_dx}")

        if np.isclose(k2, 0):
            sinh_k2_dx = dx
            cosh_k2_dx = 1
            print("k2 is close to 0; using approximations.")

        M = np.array([
            [cosh_k2_dx, sinh_k2_dx / k2],
            [k2 * sinh_k2_dx, cosh_k2_dx]
        ], dtype=complex)
    else:
        k2 = np.sqrt(2 * m * (V_segment - E)) / hbar * 1j
        print(f"k2 (imaginary): {k2}")

        cosh_k2_dx = np.cos(k2.imag * dx)
        sinh_k2_dx = 1j * np.sin(k2.imag * dx)
        print(f"cosh(k2.imag * dx): {cosh_k2_dx}, sinh(k2.imag * dx): {sinh_k2_dx}")

        if np.isclose(k2.imag, 0):
            sinh_k2_dx = dx
            cosh_k2_dx = 1
            print("k2.imag is close to 0; using approximations.")

        M = np.array([
            [cosh_k2_dx, sinh_k2_dx / k2.imag],
            [k2.imag * sinh_k2_dx, cosh_k2_dx]
        ], dtype=complex)

    print("Resulting transfer matrix M:")
    print(M)
    return M

# Print the particle energy E in Joules
print(f"Particle energy E in Joules: {E} J")


# Calculate the transfer matrix for each segment
for i in range(N - 1):
    V_segment = (V_smoothed[i] + V_smoothed[i + 1]) / 2
    M_segment = transfer_matrix_segment(V_segment, dx)
    M_total = np.dot(M_segment, M_total)

    # Check for large values in M_total and normalize if necessary
    if np.any(np.abs(M_total) > 1e3):
        print(f"Large values in M_total at step {i}, scaling down.")
        M_total /= np.max(np.abs(M_total))

# Print the total transfer matrix
print("Total Transfer Matrix (M_total):")
print(M_total)

