#################################################################
#
#   main.py
#   main file for the demonstration of the ART algorithm
#   written by: Walter Simson
#               Chair for Computer Aided Medical Procedures
#               & Augmented Reality
#               Technical University of Munich
#               27.10.2017
#   based on the work of Maximilian Baust
#
#################################################################

import warnings

import timeit

import matplotlib.pyplot as plt
import numpy as np

from art import art,art_anim
from helper import load_data

# Clean up
plt.close('all')

# Load system of equations
# A is the system matrix, b is the right hand side,
# x is the true solution, i.e. Ax=b
A, b, x = load_data("system.mat")

# Set number of iterations
iterations = 50

# Initialize plots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(6, 6))
plt.tight_layout()
ax1.set_title('np.linalg.solve')
ax2.set_title('np.linalg.qr')
ax3.set_title('ART')
ax4.set_title('True Solution')

# Plot true solution
x = x.reshape((60, 60))
ax4.imshow(x, extent=[0, 1, 0, 1], cmap='bone')

# # Solve LSE with numpy solver
# start_time = timeit.default_timer()
# x_np = np.linalg.solve(A, b)
# elapsed = timeit.default_timer() - start_time
# print('numpy took %s seconds', elapsed)
# x_np = x_np.reshape((60, 60))
# # # Warn like MATLAB - Takes too long!
# # x_cond = np.linalg.cond(A)
# # warnings.warn("Warning: Matrix is close to singular or badly scaled." +
# #               " Results may be inaccurate. Condition = {0}.".format(x_cond))

# ax1.imshow(x_np, extent=[0, 1, 0, 1], cmap='bone')

# # Solve LSE with QR-Algorithm
# start_time = timeit.default_timer()
# q, r = np.linalg.qr(A)
# y = np.matmul(q.T, b)
# x_qr = np.linalg.solve(r, y)
# elapsed = timeit.default_timer() - start_time
# print('QR took %s seconds', elapsed)
# x_qr = x_qr.reshape((60, 60))
# ax2.imshow(x_qr, extent=[0, 1, 0, 1], cmap='bone')

# Solve LSE with ART
# start_time = timeit.default_timer()
# x_art = art(A, b, iterations,1)
# elapsed = timeit.default_timer() - start_time
# print('ART took %s seconds', elapsed)

plt.show()

# x_art = art_anim(A,b,50,0.1);

fig , plts = plt.subplots(nrows=9, ncols=11,figsize=(6, 6))
plt.tight_layout();

for i in range(25,250,25):
    for j in [x * .5 for x in range(10)]:
        solution = art(A,b,i,j)
        plts[i][j].set_title("I:",i,"  R:",j)
        plts[i][j].imshow(solution, extent=[0, 1, 0, 1], cmap='bone')
plt.show();