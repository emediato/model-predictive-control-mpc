import numpy as np

# Assuming X and x(k) are given matrices
X = np.array([[2, 1],
              [1, 3]])

x_k = np.array([[1],
                [2]])

# Construct the block matrix M
M = np.block([[X, x_k],
              [x_k.T, 1]])

# Extract submatrices
X = M[:2, :2]
x_k = M[:2, 2:3]

# Calculate the Schur complement
S = X - np.dot(x_k, x_k.T)

# Check if M is negative definite
if np.all(np.linalg.eigvals(M) < 0):
    print("M is negative definite.")
else:
    print("M is not negative definite.")
