import numpy as np


def art(A, b, iterations):

    # For help with numpy (the numerical programming library for Python) check out this resource:
    # https://www.safaribooksonline.com/library/view/python-for-data/9781449323592/ch04.html
    x = np.zeros((A.shape[1]))
    # Initialize variables squared_norm (see numpy.zeros)

    # A.shape = (3600, 3600)
    # b.shape = (3600, 1)

    # Iterate over rows and compute the squared norm row-wise (we will need this in a second)
    # Hint: look into ranges

    # Iterate over iterations
    for it in range(iterations):
        # Iterate over matrix rows
        for j, row in enumerate(A):
            sqrd_norm = (np.linalg.norm(row) ** 2)
            # x' = x + correction
            if sqrd_norm <= 0:
                continue
            else:
                x += (b[j][0] - np.dot(row, x)) * np.conjugate(row) / sqrd_norm
            
    return x
