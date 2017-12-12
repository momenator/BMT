# -*- coding: utf-8 -*-
# randomWalks.m
#
# written by: Maximilian Baust & Ruediger Goebl
#               Chair for Computer Aided Medical Procedures
#               & Augmented Reality
#               Technische Universität München
#               11-19-2013
#
########################################

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from pprint import pprint


def weight(w,sigma):
    w = np.exp(-(w/sigma)**2)
    # w = 1/( 1 + (w/sigma)**4 )
    # w = 1 # Averaging.
    return w


def assemble(numbering, I, sigma):
    num_unknowns = int(round(np.max(numbering)))+1
    b = np.zeros(num_unknowns)
    A = np.zeros([num_unknowns, num_unknowns])

    for i in range(1, I.shape[0] - 1):
        for j in range(1, I.shape[1] - 1):
            if numbering[i, j] + 1 > 0:

                W = 0

                # generate entries for northern neighbor
                # compute weight
                w = weight(float(I[i, j]) - float(I[i - 1, j]), sigma)
                pprint(w)
                # update sum
                W = W + w # u*(SUM OF ALL WEIGHTS)

                # check, if there is a northern neighbor
                # with unknown potential value
                if numbering[i - 1, j] >= 0:
                    # generate entry
                    A[numbering[i, j], numbering[i - 1, j]] = -w

                # check, if there is a northern neighbor
                # with known potential value 1
                if numbering[i - 1, j] == -2:
                    # modify right hand side
                    b[numbering[i, j]] = b[numbering[i, j]] + w



                # check for southern neighbor
                # compute weight
                w = weight (float(I[i, j]) - float(I[i + 1, j]), sigma)

                # update sum
                W = W + w

                # with unknown potential value
                if numbering[i + 1, j] >= 0:
                    # generate entry
                    A[numbering[i, j], numbering[i + 1, j]] = -w

                # check, if there is a northern neighbor
                # with known potential value 1
                if numbering[i + 1, j] == -2:
                    # modify right hand side
                    b[numbering[i, j]] = b[numbering[i, j]] + w

                # check for western neighbor
                # compute weight
                w = weight (float(I[i, j]) - float(I[i, j - 1]), sigma)

                # update sum
                W = W + w

                # with unknown potential value
                if numbering[i, j - 1] >= 0:
                    # generate entry
                    A[numbering[i, j], numbering[i, j - 1]] = -w

                # check, if there is a northern neighbor
                # with known potential value 1
                if numbering[i, j - 1] == -2:
                    # modify right hand side
                    b[numbering[i, j]] = b[numbering[i, j]] + w

                # check for eastern neighbor
                # compute weight
                w = weight (float(I[i, j]) - float(I[i, j + 1]), sigma)

                # update sum
                W = W + w

                # with unknown potential value
                if numbering[i, j + 1] >= 0:
                    # generate entry
                    A[numbering[i, j], numbering[i, j + 1]] = -w

                # check, if there is a northern neighbor
                # with known potential value 1
                if numbering[i, j + 1] == -2:
                    # modify right hand side
                    b[numbering[i, j]] = b[numbering[i, j+1]] + w

                # finally, generate entry for the center pixel
                A[numbering[i, j], numbering[i, j]] = W
                '''
                That unknown number generator is cool! 
                '''
    return A, b

#ART Implementation from exercise 2.
def art(A, b, iterations):
    # For help with numpy (the numerical programming library for Python) check out this resource:
    # https://www.safaribooksonline.com/library/view/python-for-data/9781449323592/ch04.html
    x = np.zeros ((A.shape[1]))
    # Initialize variables squared_norm (see numpy.zeros)

    # A.shape = (3600, 3600)
    # b.shape = (3600, 1)

    # Iterate over rows and compute the squared norm row-wise (we will need this in a second)
    # Hint: look into ranges

    # Iterate over iterations
    for it in range (iterations):
        # Iterate over matrix rows
        for j, row in enumerate (A):
            sqrd_norm = (np.linalg.norm (row) ** 2)
            # x' = x + correction
            if sqrd_norm <= 0:
                continue
            else:
                x += ((b[j] - np.dot (row, x)) * np.conjugate (row) / sqrd_norm)

    return x



# choose sigma
# Higher sigma Results in a smoother segment boundary.
# sigma = 0.0005
# sigma = 0.05
sigma = 0.1
# sigma = 1.0
sigma = 10.0
# sigma = 1000.0
# sigma = -0.100
# sigma = -0.1
use_small_images = False

# load image
# Why div 255?
if use_small_images:
    I = np.array(Image.open("mri_small.png"))
else:
    I = np.array(Image.open("mri.png"))

# show image and select foreground
plt.figure()
plt.imshow(I, cmap="gray")
plt.title('selected foreground region')
# This could in the real world be determined by user
if use_small_images:
    R = np.array(Image.open("R_small.png"))
else:
    R = np.array(Image.open("R.png"))
R = R / np.max(R)
R = R.astype(int)
plt.contour(R, levels=[0.5], colors=['y'], linewidths=[2])

# generate mask where pixel values shall be computed
mask = np.zeros(I.shape, dtype=np.int)
mask[4:-5,4:-5] = 1
mask[R == 1] = 0

# generate numbering

# we use -1 to mark pixels with known potential 0 volts
numbering = - np.ones(I.shape, dtype=np.int)
num_unknown_pixels = np.sum(mask)
numbering[mask == 1] = np.asarray(range(0, num_unknown_pixels))

# pprint(numbering)

# we use -2 to mark pixels with known potential 1 volts
numbering[R == 1] = -2
plt.figure()
plt.imshow(numbering, cmap='gray')
plt.title('single pixel numbering')

# assemble random walker equation system
[A, b] = assemble(numbering, I, sigma)

# pprint(A)
# pprint(b)
# compute solution
# ART implementation from exercise 2.
x = art(A, b, 50)

# Cant use linalg solver for small values of sigma since that produces a singular matrix.
# x = np.linalg.solve(A, b)

print(x)

# generate segmentation
seg = np.zeros(I.shape)
seg[R == 1] = 1
seg[numbering >= 0] = x

# visualize segmentation
plt.figure()
plt.imshow(I, cmap="gray")
plt.contour(seg, levels=[0.5], colors=['y'], linewidths=[2])

# leftovers from Matlab!
# imagesc(I)
#colormap gray
#axis equal tight off
#hold on
#contour(seg,[0.5, 0.5],'y','LineWidth',2);
#hold off

plt.show()