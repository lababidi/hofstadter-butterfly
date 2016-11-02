""" A simple Mandlebrot """

def mandlebrot(b, z):
    """Computes the mandlebrot"""
    return b * b + z

Z = .075 + .043*1j
B = Z
DEC = 10
while DEC > 0:
    B = mandlebrot(B, Z)
    print(B, abs(B))
    DEC -= 1
