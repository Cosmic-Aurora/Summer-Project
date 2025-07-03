import numpy as np

arr1 = np.array([[1,2],[1,3],[2,3]])
arr2 = np.array([[2,3], [1,4]])

print(np.intersect1d(arr1, arr2, assume_unique = True))
print(arr1.includes([1,2]))