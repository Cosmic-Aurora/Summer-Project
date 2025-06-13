import numpy as np

arr = np.array([0, 0, 0, 0, 0, 0, 8, 83, 120, 111, 31, 37, 10, 0, 0, 0, 0, 0, 0, 0])

sak = np.matrix([[0,1,0,4],[3,0,2,1],[5,6,0,0]])
sak = np.pad(sak,2)
print(sak)
s = sak != 0
print(s)
print(np.argmax(s[::-1],axis = 0))
print(np.argmax(s[::],axis = 1))


m = arr!=0
print(m)
print(m.argmax())

arr = np.trim_zeros(arr)
print(arr)