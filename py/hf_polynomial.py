# Hofstadter Butterfly Fractal
# http://en.wikipedia.org/wiki/Hofstadter%27s_butterfly
# Wolfgang Kinzel/Georg Reents,"Physics by Computer" Springer Press (1998)
# FB36 - 20130922
import math
import matplotlib.pyplot as plt
import numpy as np

imgSize = 400

pixels = np.zeros((imgSize, imgSize))


def gcd(a, b):  # Greatest Common Divisor
    return a if b == 0 else gcd(b, a % b)

def poly_cal(poly, poly_old, n):
    for m in range(2, q // 2):
        poly_new = (2.0 * math.cos(sigma * m) - e) * poly - poly_old
        if poly * poly_new < 0.0:
            n += 1
        poly_old = poly
        poly = poly_new
    return poly, poly_old, n


pi2 = math.pi * 2.0
MAXX = imgSize + 1
MAXY = imgSize + 1
qmax = 100  # imgSize
for q in range(2, qmax):
    print(str(100 * q / qmax).zfill(2) + "%")
    for p in range(1, q):
        if gcd(p, q) == 1:
            sigma = pi2 * p / q
            n_old = 0
            ie = 0
            for ie in range(0, MAXY + 2):
                e = 8.0 * ie / MAXY - 4.0 - 4.0 / MAXY
                n = 0
                poly_old = 1.0
                poly = 2.0 * math.cos(sigma) - e
                if poly_old * poly < 0.0:
                    n += 1

                # for m in range(2, q // 2):
                #     poly_new = (2.0 * math.cos(sigma * m) - e) * poly - poly_old
                #     if poly * poly_new < 0.0:
                #         n += 1
                #     poly_old = poly
                #     poly = poly_new
                poly, poly_old, n = poly_cal(poly, poly_old, n)

                poly_old = 1.0
                poly = 2.0 - e
                if poly_old * poly < 0.0:
                    n += 1
                poly_new = (2.0 * math.cos(sigma) - e) * poly - 2.0 * poly_old

                if poly * poly_new < 0.0:
                    n += 1
                poly_old = poly
                poly = poly_new

                poly, poly_old, n = poly_cal(poly, poly_old, n)
                # for m in range(2, q // 2):
                #     poly_new = (2.0 * math.cos(sigma * m) - e) * poly - poly_old
                #     if poly * poly_new < 0.0:
                #         n += 1
                #     poly_old = poly
                #     poly = poly_new

                poly_new = (2.0 * math.cos(sigma * q / 2.0) - e) * poly - 2.0 * poly_old
                if poly * poly_new < 0.0:
                    n += 1

                poly_old = 1.0
                poly = 2.0 - e
                if poly_old * poly < 0.0:
                    n += 1
                poly_new = (2.0 * math.cos(sigma) - e) * poly - 2.0 * poly_old
                if poly * poly_new < 0.0:
                    n += 1
                poly_old = poly
                poly = poly_new

                poly, poly_old, n = poly_cal(poly, poly_old, n)
                # for m in range(2, q // 2):
                #     poly_new = (2.0 * math.cos(sigma * m) - e) * poly - poly_old
                #     if poly * poly_new < 0.0:
                #         n += 1
                #     poly_old = poly
                #     poly = poly_new

                poly_new = (2.0 * math.cos(sigma * q / 2.0) - e) * poly - 2.0 * poly_old
                if poly * poly_new < 0.0:
                    n += 1
                if n > n_old:
                    # pixels[int(MAXY - ie)][int(MAXX * p / q)] = -1  # (255, 255, 255)
                    pixels[int(MAXX * p / q)][int(MAXY - ie)] = -1  # (255, 255, 255)
                n_old = n

plt.matshow(pixels, cmap=plt.cm.gray)
plt.show()
