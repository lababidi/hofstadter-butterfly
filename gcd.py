import numpy as np
import matplotlib.pyplot as plt
def gcd(a,b): return  gcd(b, a%b) if b else a
gcd_v = np.vectorize(gcd)
n = 1000
a = [[p]*n for p in range(1,n+1)]
print(gcd_v(a,np.transpose(a)))
