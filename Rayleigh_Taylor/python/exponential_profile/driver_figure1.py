import numpy as np
import matplotlib.pyplot as plt
from get_eigenvalues import get_eigenvalues

# Translated from Matlab to Python with ChatGPT

# Clear screen equivalent (optional)
print("Plot of the top eigenvalue as a function of k")

# Parameters
p = {}
p["rho_m"] = 0.1
p["rho_p"] = 1.0  # bigger than rho_m
p["h"] = 1.5
p["L"] = 1.0

# Don't change
p["g"] = 9.8
p["a"] = np.log(p["rho_p"] / p["rho_m"]) / p["h"]

profile_ratio = lambda x: p["a"]

# Atwood number
p["Atwood_number"] = (p["rho_p"] - p["rho_m"]) / (p["rho_p"] + p["rho_m"])
print(f"\nAtwood number = {p['Atwood_number']:.4g}")

bound = np.sqrt(p["a"] * p["g"])

# Compute eigenvalues
num_k_indices = 200
num_eig_indices = 1

M = np.zeros((num_k_indices, num_eig_indices))
k_vals = np.zeros(num_k_indices)


for j in range(num_k_indices):
    for index in range(num_eig_indices):
        k = (j + 1) / p["L"]  # MATLAB indexing starts at 1
        k_vals[j] = k
        M[j, index] = get_eigenvalues(index + 1, k, p["a"], p["h"])

# Plot
plt.figure()

for j in range(M.shape[1]):
    plt.plot(k_vals, M[:, j], '.k', markersize=8)

plt.plot(
    [0, num_k_indices],
    [bound, bound],
    '--b',
    linewidth=2
)

mx = np.max(M)

plt.axis([
    -1,
    num_k_indices,
    -0.1 * mx,
    1.1 * mx
])

plt.xlabel('k', fontsize=18)
plt.ylabel(r'$\lambda$', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.show()