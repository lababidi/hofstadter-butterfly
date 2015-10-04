import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt

q_max = 50


def gcd(a, b):
    return gcd(b, a % b) if b else a


def unique(p, q):
    return gcd(p, q) == 1


def krond(a, b):
    return (a == b) * 1


def main():
    j1 = 1.0
    j2 = 1.0
    res = []
    res_unique = set()
    for q in range(2, q_max):
        for p in range(1, q):
            if unique(p, q):
                kx = 0  # should 0 to 2pi/q
                ky = 0  # should be 0 to pi/q
                qr = np.arange(0, q, 1)
                i, j = np.meshgrid(qr, qr)
                phi = 1.0 * p / q
                h1 = krond(i, j) * 2 * np.cos(ky + 2 * np.pi * (j + 1) * phi)
                h2 = krond(i + 1, j) + krond(i + q - 1, j) * np.exp(-1j * (q * kx))
                h2 += krond(i, j + 1) + krond(i, j + q - 1) * np.exp(1j * (q * kx))
                ee = np.dot(linalg.expm(-1j * h1 * j1), linalg.expm(-1j * h2 * j2))
                eigs = np.angle(linalg.eigvals(ee)) / np.pi
                # phis = np.ones_like(eigs)
                for eig in eigs:
                    res.append((phi, eig))
                    res_unique.add((phi, eig))
    plt.scatter(*zip(*res), s=.2)
    plt.show()
    print(len(res_unique))


if __name__ == '__main__':
    main()
