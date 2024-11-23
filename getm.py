#!/usr/bin/env python3
"""
output refractive index given frequency/wavelength and temperature
inputs: -fghz FinGHz -wl WLum -tc Tc -tk TK -water -ice
when both wavelength and frequecy are both given, the one last given will be used.
when both Tc and TK are given, the one last given will be used.
for class MET5471 - G. Liu 2024.10.15
"""
import sys
import riwaterice as ri
#speed of light in vacuum
c=2.99792458e8
wl=30.e3
tc=10
what=0

narg = len(sys.argv)
i = 0
while i<narg:
	if sys.argv[i] == "-fghz":
		i += 1
		fghz = float(sys.argv[i])
		wl = c/fghz*1.e-3
	elif sys.argv[i] == "-wl":
		i += 1
		wl = float(sys.argv[i])
	elif sys.argv[i] == "-tc":
		i += 1
		tc = float(sys.argv[i])
	elif sys.argv[i] == "-tk":
		i += 1
		tc =float(sys.argv[i]) - 273.16
	elif sys.argv[i] == "-ice":
		what=1
	elif sys.argv[i] == "-water":
		what=0
	elif sys.argv[i] == "-h":
		print("Usage: getm.py [-fghz F | -wl lambda(um)] [-tc Tc | -tk TK] [-water]|-ice",file=sys.stderr)
		exit()
	i += 1

if what == 0:
	m = ri.mwater(wl,tc)
	print("Refractive index for water at t=%.2f C and wl=%.2f um:" % (tc, wl),file=sys.stderr)  
else:
	m = ri.mice(wl,tc)
	print("Refractive index for ice at t=%.2f C and wl=%.2f um:" % (tc, wl),file=sys.stderr)

#print(m)
if abs(m.imag)>0.1:
	print("%.4f,%.4f" % (m.real, m.imag))
else:
	print("%.4f,%.4e" % (m.real, m.imag))
