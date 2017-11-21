import numpy as np
import timeit

def projection(u, x): return (np.dot(x, u) / np.dot(u, u)) * u

""" 
	Implementation of QR decomposition algorithm 
	Input: invertible matrix A (M, N)
	Output: 2 matrix - ortogonal matrix Q (M, N) and 
	upper right triangular matrix R (N, N)

	M and N should be the same
"""
def my_qr(A):
	N = A.shape[1]
	# probably better to declare with np array
	q = []
	R = []
	
	""" 
		iterate columns of A
			compute vector u
			compute orthonormal basis q for each column, this column q forms the matrix Q
			compute for nth row for N
	"""
	for i in range(N):
		# compute for ui
		ui = A[:, i]
		curr_r_row = []
		
		if i != 0:
			ui = A[:, i] - sum([ projection(qi, A[:, i]) for qi in q ])

		ui_norm = np.linalg.norm(ui)

		# compute for qi and append to array q
		qi = ui / ui_norm
		q.append(qi)
		
		# build row n of R
		for j in range(N):
			if i == j:
				curr_r_row.append(ui_norm)
			elif j < i:
				curr_r_row.append(0)
			else:
				curr_r_row.append(np.dot(A[:,j], qi))
			
		R.append(curr_r_row)

	# format u that forms Q
	Q = np.array(q).T
	R = np.array(R)

	return [Q, R]

"""
A = np.array([[15, 3, 3, 0], [-1, 4, 3, 56], [-9, 2, 3, 2], [1, 0, 2, 7]])
my_qr(A)

[Q1, R1] = my_qr(A)
print(Q1, R1)

[Q, R] = np.linalg.qr(A)
print('linalg q ', Q)
print('linalg r ', R)
"""
ep = 1e-9
A = np.array([[ep, 1], [0, ep], [1, 1]])
b = np.array([1, 2, 0])

start_time = timeit.default_timer()
x = np.linalg.solve(np.matmul(A.T,A),np.matmul(A.T,b))
elapsed = timeit.default_timer() - start_time
#print('using numpy solve : {}, took {} seconds'.format(x, elapsed))

start_time = timeit.default_timer()
q, r = np.linalg.qr(A)
print('NUMPY')
print(q, r)
y = np.matmul(q.T, b)
x_qr = np.linalg.solve(r, y)
elapsed = timeit.default_timer() - start_time
print('using numpy qr : {}, took {} seconds'.format(x_qr, elapsed))

start_time = timeit.default_timer()
q2, r2 = my_qr(A)
print('OWN')
print(q2, r2)
y2 = np.matmul(q2.T, b)
x_qr2 = np.linalg.solve(r2, y2)
elapsed = timeit.default_timer() - start_time
print('using my_qr : {}, took {} seconds'.format(x_qr2, elapsed))