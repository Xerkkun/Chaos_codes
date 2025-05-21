import numpy as np

A = np.eye( 5 )

A = np.delete(A, [3,4], axis = 0)
print(A)
