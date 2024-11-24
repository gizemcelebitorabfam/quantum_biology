import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


# Define the three transfer matrices

M1 = np.array([[1, 5e-14],
               [2.95e+9, 1]])

M2 = np.array([[1, 5e-14],
               [6.69e+9, 1]])

M3 = np.array([[1, 5e-14],
               [5.2e+9, 1]])

# Compute the total transfer matrix
M_total = np.dot(M3, np.dot(M2, M1))

print("Total Transfer Matrix:\n", M_total)
