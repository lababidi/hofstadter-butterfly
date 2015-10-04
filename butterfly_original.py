import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt

_pi = np.pi
q_max = 40
upper_bound = 1
lower_bound = 0
# upper_bound = 7.0 / 19
# lower_bound = 4.0 / 11


def gcd(a, b):
    return gcd(b, a % b) if b else a


def unique(p, q):
    return gcd(p, q) == 1


def krond(a, b):
    return (a == b) * 1


def in_range(p, q):
    return lower_bound * q <= p <= upper_bound * q


def main():
    j1 = 1.0
    j2 = 1.0
    res = []
    with open('butterfly2.csv', 'w') as f:
        f.write('x,y')
    for q in range(1, q_max):
        print(q)
        kr = np.arccos(np.arange(-1.0, 1.0, .2)) / q
        for p in range(2, q):
            if unique(p, q) and in_range(p, q):
                # dk = .1 * np.pi / q
                phi = 1.0 * p / q
                for kx in kr:
                    ki = len(kr)
                    for ky in np.nditer(kr[0:ki]):
                        qr = np.arange(0, q, 1)
                        i_, j_ = np.meshgrid(qr, qr)
                        h1 = krond(i_, j_) * 2 * np.cos(ky + 2 * np.pi * (j_ + 1) * phi)
                        h2 = krond(i_ + 1, j_) + krond(i_, j_ + 1) + krond(i_ + q - 1, j_) * np.exp(-1j * (q * kx))
                        h2 += krond(i_, j_ + q - 1) * np.exp(1j * (q * kx))
                        # ee = np.dot(linalg.expm(-1j*h1*j1) , linalg.expm(-1j*h2*j2))
                        eigs = linalg.eigvals(np.add(j1 * h1, j2 * h2))
                        for eig in eigs:
                            res.append((phi, eig))
        with open('butterfly2.csv', 'a') as f:
            f.writelines([str(r[0]) + ',' + str(r[1]) + '\n' for r in res])
        # res = []
    print(res)
    plt.scatter(*zip(*res), s=.2)
    plt.show()

if __name__ == '__main__':
    main()
