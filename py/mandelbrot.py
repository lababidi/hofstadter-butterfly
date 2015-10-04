import cmath

def m(b,z):
	return b*b+z

z = .075+ .043*1j
b = z
dec = 10
while dec>0:
	b = m(b,z)
	print b,abs(b)
	dec-=1

