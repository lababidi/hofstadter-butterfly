import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt

q_max = 30

def  gcd(a,b): return  gcd(b, a%b) if b else a
def unique(p,q): return gcd(p,q)==1
def krond(a,b): return (a==b)*1

J1 = 2.0
J2 = 2.0
res = []
for q in range(2,q_max):
	for p in range(1,q):
		if unique(p,q):
			kx = 0 #should 0 to 2pi/q
			ky = 0 #should be 0 to pi/q
			qr = np.arange(0,q,1)
			i_ , j_  = np.meshgrid(qr,qr)
			phi = 1.0 * p / q
			H1 = krond(i_,j_)*2*np.cos(ky + 2* np.pi * (j_+1) * phi)
			H2 = krond(i_+1,j_)+krond(i_,j_+1)+krond(i_+q-1,j_)*np.exp(-1j*(q*kx))+krond(i_,j_+q-1)*np.exp(1j*(q*kx));
			ee = np.dot(linalg.expm(-1j*H1*J1) , linalg.expm(-1j*H2*J2))
			eigs = np.angle(linalg.eigvals(ee))/np.pi
			for eig in eigs:
				res.append((phi,eig))
#plt.scatter(*zip(*res),s=.2)
#plt.show()
