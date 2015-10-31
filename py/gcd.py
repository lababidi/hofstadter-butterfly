"""GCD script"""

from __future__ import print_function

#pylint: disable=no-member

import numpy as np

def gcd(a, b):
    """Do the GCD calculation"""
    return gcd(b, a%b) if b else a

def main():
    """Run the script"""
    gcd_v = np.vectorize(gcd)
    count = 1000
    array = [[p] * count for p in range(1, count+1)]
    print(gcd_v(array, np.transpose(array)))

if __name__ == "__main__":
    main()
