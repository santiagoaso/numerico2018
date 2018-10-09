import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

L = 1
E_I = 1
nes = [5, 10, 100]
nro_grupo = 18

def SOR(K, f, Imax, eps, u, omega):
    D = np.diag(np.diag(K))
    M = np.dot((1 / omega), D) + np.tril(K, -1)
    N = M - K
    r = np.dot(K, u) - f
    x = u
    i = 0
    err = 1 + eps
    res = []
    cut = err
    while((i < Imax) and (cut >= eps)):
        x_old = x
        x = np.dot(np.dot((la.inv(M)), N), x) + np.dot((la.inv(M)), f)
        r = np.dot(K, x) - f
        err = la.norm(r, np.inf)
        res.append(err)
        cut = (la.norm(x - x_old, np.inf)) / (la.norm(x, np.inf))
        i = i + 1
    # print("n = " + str(len(f)) + " converge en " + str(i) + " pasos para omega: " + str(omega))
    return (x, i, res)

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
		newElement = (((nro_grupo + (nro_grupo**2) * ((i * L/n) - (i * L/n)**2)) / E_I ) * ((L/n)**4))
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
	current_n = nes[x]
	K = CreateMatrixK(current_n)
	u = CreateVectorU(current_n)
	f = CreateVectorF(current_n)

	def CalculateBestOmega():
		l = 1.0
		bestOmega = 0
		recordSteps = 999
		while (l < 2):
		    [x, i, res] = SOR(K, f, 10000, 10**(-2), u, l)
		    if (len(res) < recordSteps):
		    	bestOmega = l
		    	recordSteps = len(res)
		    l = l + 0.05
		return [bestOmega, recordSteps]

	[bestOmega, recordSteps] = CalculateBestOmega()

	def Graph(omega):
		[x, i, res] = SOR(K, f, 10000, 10**(-4), u, omega)
		plt.plot(range(len(res)),res)
		plt.title("n = " + str(current_n) + ". Omega = " + str(bestOmega) + ". Steps: " + str(len(res)))
		plt.yscale('log')
		plt.ylabel('Infinity Norm (log)')
		plt.xlabel('Step')
		plt.show()

	Graph(bestOmega)

	x = x + 1

