import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

L = 1
E_I = 1

def SOR(A, b, Imax, eps, x0, omega):
    D = np.diag(np.diag(A))
    M = np.dot((1/omega), D) + np.tril(A, -1)
    N = M - A
    r = np.dot(A,x0)-b
    x = x0
    i = 0
    err = 1 + eps
    res = []
    while((i < Imax) and ((la.norm(r)) >= eps)):
        x = np.dot(np.dot((la.inv(M)),N),x) + np.dot((la.inv(M)),b)
        r = np.dot(A,x) - b
        err = la.norm(r,2)
        res.append(err)
        i = i+1
    return (x,i,res)

def CreateMatrixK():
	matrixK = np.array([
		[1,0,0,0,0,0,0,0,0,0,0,0],
		[-4,5,-4,1,0,0,0,0,0,0,0,0],
		[1,-4,6,-4,1,0,0,0,0,0,0,0],
		[0,1,-4,6,-4,1,0,0,0,0,0,0],
		[0,0,1,-4,6,-4,1,0,0,0,0,0],
		[0,0,0,1,-4,6,-4,1,0,0,0,0],
		[0,0,0,0,1,-4,6,-4,1,0,0,0],
		[0,0,0,0,0,1,-4,6,-4,1,0,0],
		[0,0,0,0,0,0,1,-4,6,-4,1,0],
		[0,0,0,0,0,0,0,1,-4,6,-4,1],
		[0,0,0,0,0,0,0,0,1,-4,5,-4],
		[0,0,0,0,0,0,0,0,0,0,0,1]
		])
	return matrixK

def CreateVectorF(n):
	vectorF = np.array(0)
	i = 1
	while (i < n):
		newElement = (((18 + (18**2) * ((L/n)-(L/n)**2))/E_I ) * ((L/n)**4))
		vectorF = np.append(vectorF, newElement)
		i = i + 1 
	vectorF = np.append(vectorF, 0)
	return vectorF

def CreateVectorU(n):
	vectorU = np.array(0)
	i = 1
	while (i <= n):
		newElement = 0
		vectorU = np.append(vectorU, newElement)
		i = i + 1
	return vectorU

n = 11
K = CreateMatrixK()
u = CreateVectorU(n)  
f = CreateVectorF(n)
print(K)
print(u)
print(f)

def Graph():
	i = 1.0
	while (i < 2):
	    [x,i,res] = SOR(K,f,100,10**(-3),u,i)
	    plt.plot(range(len(res)),res)
	    plt.yscale('log')
	    plt.show()
	    i = i + 0.05

Graph()

