from math import gcd
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import argparse
from numpy import isclose as k_d

parser = argparse.ArgumentParser(description='Calculate Hofstadter\'s butterfly.')
parser.add_argument("--file_output", help="Name of file to output", default="data/")
parser.add_argument("--q_max", help="The maximum value of q (p/q rational)", default=40, type=int)
parser.add_argument("--upper_bound", help="The maximum value p/q", default=0.0, type=float)
parser.add_argument("--lower_bound", help="The minimum value of p/q", default=1.0, type=float)
parser.add_argument("--filetype", help="(csv,tsv)", default="csv")
parser.add_argument("--showplot", help="show plots while code runs", default=False, type=bool)
args = parser.parse_args()
separation_char = ',' if args.filetype == "csv" else '\t'

def in_range(p_: int, q_: int): return args.lower_bound * q_ <= p_ <= args.upper_bound * q_


J1, J2 = 1.0, 1.0
res, total_res = set(), []

for J1 in np.arange(1.90, 4.01, .1):
    J2=J1
    fname = args.file_output + str(J1) + "-" + str(J2) + ".data"
    with open(fname, 'w') as f:
        f.write('x' + separation_char + 'y\n')
    for q in range(2, args.q_max):
        print(q)
        kr = np.arange(-0.50, .51, .1) * np.pi / q
        kyr = np.arange(-0.0, 1.01, .1) * np.pi / q
        kxr = np.arange(0.0, 1.01, .1) * np.pi / q
        i, j = np.meshgrid(np.arange(0, q, 1), np.arange(0, q, 1))
        one = np.ones_like(i)
        q1 = q * one
        for p in range(1, q):
            if gcd(p, q) == 1 and p <= q / 2:  # and in_range(p,q):
                for kx in kxr:
                    for ky in kyr:
                        H1 = k_d(i, j) * 2 * np.cos(ky + 2 * np.pi * (p / q) * (j + 1))
                        H2 = k_d(i + one, j) + k_d(i + q1 - one, j) * np.exp(-1j * q * kx)
                        H2 += k_d(i, j + one) + k_d(i, j + q1 - one) * np.exp(1j * q * kx)
                        # eigs = linalg.eigvals(np.add(H1, H2))
                        eigs = linalg.eigvals(
                            np.dot(linalg.expm(-1j * H1 * J1 / np.pi), linalg.expm(-1j * H2 * J2 / np.pi)))
                        eigs = np.arctan2(eigs.imag, eigs.real)
                        for eig in eigs:
                            res.add(((1.0 * p / q), eig.real))
                            res.add((1 - (1.0 * p / q), eig.real))
                            res.add(((1.0 * p / q), -eig.real))
                            res.add((1 - (1.0 * p / q), -eig.real))
                            if args.showplot:
                                total_res.append(((1.0 * p / q), eig.real))
                                total_res.append((1 - (1.0 * p / q), eig.real))
        with open(fname, 'a') as f:
            f.writelines([str(r[0]) + separation_char + str(r[1]) + '\n' for r in res])
        res.clear()
        if args.showplot and q > 19 and q % 5 == 0:
            plt.scatter(*zip(*total_res), s=.2)
            plt.show()
