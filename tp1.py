import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

L = 1
E_I = 1
nes = [5, 10, 100]

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
    print("n = " + str(len(b)) + " converge en " + str(i) + " pasos")
    return (x,i,res)

def CreateLineMatrixK(j, n):
	if (j == 1):
		newLine = [1]
		i = 1
		while (i < n):
			newLine.append(0)
			i = i + 1
		return newLine
	if (j == 2):
		newLine = [-4, 5, -4, 1]
		i = 4
		while (i < n):
			newLine.append(0)
			i = i + 1
		return newLine
	if (j >= 3 and j <= (n-2)):
		newLine = []
		i = 0
		while (i < (j-3)):
			newLine.append(0)
			i = i + 1
		newLine.append(1)
		newLine.append(-4)
		newLine.append(6)
		newLine.append(-4)
		newLine.append(1)
		i = i + 5
		while (i < n):
			newLine.append(0)
			i = i + 1
		return newLine
	if (j == (n-1)):
		newLine = []
		i = 0
		while (i < (n - 4)):
			newLine.append(0)
			i = i + 1
		newLine.append(1)
		newLine.append(-4)
		newLine.append(5)
		newLine.append(-4)
		i = i + 4
		return newLine	
	if (j == (n)):
		newLine = []
		i = 0
		while (i < (n - 1)):
			newLine.append(0)
			i = i + 1
		newLine.append(1)
		return newLine

def CreateMatrixK(n):
	matrixK = []
	i = 1
	while (i <= n):
		newLine = CreateLineMatrixK(i,n)
		matrixK.append(newLine)
		i = i + 1 
	return matrixK

def CreateVectorF(n):
	vectorF = []
	vectorF.append(0)
	i = 1
	while (i < n-1):
		newElement = (((18 + (18**2) * ((i*L/n)-(i*L/n)**2))/E_I ) * ((L/n)**4))
		vectorF.append(newElement)
		i = i + 1 
	vectorF.append(0)
	return vectorF

def CreateVectorU(n):
	vectorU = []
	vectorU.append(0)
	i = 1
	while (i < n):
		newElement = 1
		vectorU.append(newElement)
		i = i + 1
	return vectorU

x = 0
while (x < len(nes)):
	K = CreateMatrixK(nes[x])
	u = CreateVectorU(nes[x])
	f = CreateVectorF(nes[x])

	def Graph():
		i = 1.0
		while (i < 2):
		    [x,i,res] = SOR(K, f, 10000, 10**(-2), u, i)
		    plt.plot(range(len(res)),res)
		    plt.yscale('log')
		    plt.show()
		    i = i + 0.05

	Graph()
	x = x + 1

