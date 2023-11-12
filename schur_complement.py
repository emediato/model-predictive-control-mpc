import numpy as np

# Define matrices A, B, C, and D
A = np.array([[2, 1],
              [1, 3]])

B = np.array([[1, -1],
              [0, 2]])

C = np.array([[3, 2],
              [4, 1]])

D = np.array([[2, 0],
              [0, 1]])

# Calculate the Schur complement
Schur_complement = A - np.dot(np.dot(B, np.linalg.inv(D)), C)

# Display the results
print("Matrix A:")
print(A)

print("\nMatrix B:")
print(B)

print("\nMatrix C:")
print(C)

print("\nMatrix D:")
print(D)

print("\nSchur Complement:")
print(Schur_complement)
