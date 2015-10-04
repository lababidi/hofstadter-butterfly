import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Calculate Hofstadter\'s butterfly.')
parser.add_argument("--file_output", help="Name of file to output", default="butterfly.data")
parser.add_argument("--q_max", help="The maximum value of q (p/q rational)", default=100, type=int)
parser.add_argument("--filetype", help="(csv,tsv)", default="tsv")
args = parser.parse_args()
separation_char = ',' if args.filetype == "csv" else '\t'

upper_bound = 7.0/19
lower_bound = 4.0/11
def gcd(a,b): return  gcd(b, a%b) if b else a
def unique(p,q): return gcd(p,q)==1
def krond(a,b): return (a==b)*1
def in_range(p,q): return p>=lower_bound *q and p <= upper_bound * q
J1 = 1.0
J2 = 1.0
res = set()
with open(args.file_output,'w') as f:
    if args.filetype=="csv":
	f.write('x,y')
for q in range(1,args.q_max):
	print(q)
	kr = np.arccos(np.arange(-1.0,1.0,.2))/q
	for p in range(2,q):
		if unique(p,q): # and in_range(p,q):
			dk = .1*np.pi/q
			phi = 1.0 * p / q
			for kx in kr:
				ki = len(kr)
				for ky in np.nditer(kr[0:ki]):
					qr = np.arange(0,q,1)
					i_ , j_  = np.meshgrid(qr,qr)
					H1 = krond(i_,j_)*2*np.cos(ky + 2* np.pi * (j_+1) * phi)
					H2 = krond(i_+1,j_)+krond(i_,j_+1)+krond(i_+q-1,j_)*np.exp(-1j*(q*kx))+krond(i_,j_+q-1)*np.exp(1j*(q*kx));
					#ee = np.dot(linalg.expm(-1j*H1*J1) , linalg.expm(-1j*H2*J2))
					eigs = linalg.eigvals(np.add(H1,H2))
					for eig in eigs:
						res.add((phi,eig.real))
	with open(args.file_output,'a') as f:
		f.writelines([str(r[0])+separation_char+str(r[1])+'\n' for r in res])
        res.clear()
