#!/usr/bin/env python3
"""
modules to compute refractive indices of water and ice
1.Water for microwave and far infrared (0 - 25 THz):
		Ellison, W. J., J. Phys. Chem. Ref. Data, 2007: Permittivity of pure water, at standard atmospheric
			 pressure, over the frequency range 0-25 THz and the temperature range 0 - 100C,36.
			 doi:10.1063/1.2360986
2.Water for visible and infrared (0.2 to 200 um):
		Hale, G. M. and M. R. Querry, 1973: Optical constants of water in the 200-nm to 200-um wavelength region.
		App. Opt., 12(3), 555-563. doi: 10.1364/ao.12.000555
3. Ice for visible and infrared: 
		Warren, S. G., and R. E. Brandt (2008), Optical constants of ice from the ultraviolet to the microwave: 
		A revised compilation. J. Geophys. Res., 113, D14220, doi:10.1029/2007JD009744
4. Ice for microwave
		Ma ̈tzler,C.(2006),Microwavedielectricpropertiesofice,inThermal Microwave Radiation—Applications for 
		Remote Sensing, Electromagn, Waves Ser., vol. 52, edited by C. Ma ̈tzler et al., chap. 5, 455 – 462,
		Inst. Eng. Technol., Stevenage, U. K.
"""
import numpy
from numpy import exp, log, loadtxt
import cmath

# for water. note: in vis-ir, temperature t is not used. only for 25C
# for ice. note: in vis_ir, temperature t is not used, only for -7C
# wavelength wl in um, temperature t in C
def mwater(wl,t):
	if wl<100:
		return mwater_visir(wl)
	else:
		return mwater_mw(wl,t)
def mice(wl,t):
	if wl<100.:
		return mice_visir(wl)
	else:
		return mice_mw(wl,t)
# water fo vis-ir
def mwater_visir(wl):
	wdata = loadtxt("WOP_1973_HaleQuerry_umkn.dat",skiprows=1)
	wls = wdata[:,0]
	ns = wdata[:,2]
	ks = wdata[:,1]
	nl = len(wls)
	m = complex(0.0,0.0)
	i = -1
	for j in range(nl-1):
		if wl>= wls[j] and wl< wls[j+1]:
			i = j
	if i<0:
		return m
	n = ns[i] + (ns[i+1]-ns[i])/(wls[i+1]-wls[i])*(wl-wls[i])
	k = 0.1 #ks[i] + (ks[i+1]-ks[i])/(wls[i+1]-wls[i])*(wl-wls[i]) # !!!!!!!!!!!!! k modified here. -cansu
	m = complex(n, -k)
	return m

# ice for vis-ir
def mice_visir(wl):
	idata = loadtxt("IOP_2008_WarrenBrandt_umnk.dat",skiprows=1)
	wls = idata[:,0]
	ns = idata[:,1]
	ks = idata[:,2]
	nl = len(wls)
	m = complex(0.0,0.0)
	i = -1
	for j in range(nl-1):
		if wl>= wls[j] and wl< wls[j+1]:
			i = j
	if i<0:
		return m
#	n = ns[i] + (ns[i+1]-ns[i])/(wls[i+1]-wls[i])*(wl-wls[i])
	n = ns[i] + (ns[i+1]-ns[i])/(log(wls[i+1])-log(wls[i]))*(log(wl)-log(wls[i]))
#	k = ks[i] + (ks[i+1]-ks[i])/(wls[i+1]-wls[i])*(wl-wls[i])
	logk = log(ks[i]) + (log(ks[i+1])-log(ks[i]))/(log(wls[i+1])-log(wls[i]))*(log(wl)-log(wls[i]))
	k = exp(logk)
	m = complex(n, -k)
	return m

# ice for microwave
def mice_mw(wl,t):
	c=2.99792458e8
# um to GHz
	fghz = c/wl*1.e-3
# C to K
	tk = t + 273.16
	if tk>243.:
		er = 3.1884+9.1e-4*(tk-273)
	else:
		er = 3.1611+4.3e-4*(tk-243)
	theta = 300./tk - 1
	alpha = (0.00504+0.0062*theta)*exp(-22.1*theta)
	dbeta = exp(-9.963+0.0372*(tk-273.16))
	B1=0.0207
	B2=1.16E-11
	b=335.
	betam=B1/tk*exp(b/tk)/(exp(b/tk)-1)**2+B2*fghz**2
	beta=betam+dbeta
	ei=alpha/fghz+beta*fghz
	e=complex(er,-ei)
	m=cmath.sqrt(e)
	return m

# water for microwave
def mwater_mw(wl,t):
	c=2.99792458e8
# frequency in Hz
	v = c/wl*1.e6
#constants in Table 2
	a1 = 79.23882
	a2 = 3.815866
	a3 = 1.634967
	tc = 133.1383
	b1 = 0.004300598
	b2 = 0.01117295
	b3 = 0.006841548
	c1 = 1.382264e-13
	c2 = 3.510354e-16
	c3 = 6.30035e-15
	d1 = 652.7648
	d2 = 1249.533
	d3 = 405.5169
	p0 = 0.8379692
	p1 = -0.006118594
	p2 = -0.000012936798
	p3 = 4235901000000.
	p4 = -14260880000.
	p5 = 273815700.
	p6 = -1246943.
	p7 = 9.618642e-14
	p8 = 1.795786e-16
	p9 = -9.310017e-18
	p10 = 1.655473e-19
	p11 = 0.6165532
	p12 = 0.007238532
	p13 = -0.00009523366
	p14 = 15983170000000.
	p15 = -74413570000.
	p16 = 497448000.
	p17 = 2.882476e-14
	p18 = -3.142118e-16
	p19 = 3.528051e-18

	es = 87.9144-0.404399*t+9.58726e-4*t*t-1.32802e-6*t*t*t
	dt1 = a1*exp(-b1*t)
	dt2 = a2*exp(-b2*t)
	dt3 = a3*exp(-b3*t)
	tt1 = c1*exp(d1/(t+tc))
	tt2 = c2*exp(d2/(t+tc))
	tt3 = c3*exp(d3/(t+tc))

	dt4 = p0+p1*t+p2*t*t
	ft0 = p3+p4*t+p5*t*t+p6*t*t*t
	tt4 = p7+p8*t+p9*t*t+p10*t*t*t
	dt5 = p11+p12*t+p13*t*t
	ft1 = p14+p15*t+p16*t*t
	tt5 = p17+p18*t+p19*t*t

	pi = numpy.pi
	er = es - (2*pi*v)**2*(tt1**2*dt1/(1+(2*pi*v*tt1)**2) + tt2**2*dt2/(1+(2*pi*v*tt2)**2) + tt3**2*dt3/(1+(2*pi*v*tt3)**2))
	- (2*pi*tt4)**2*dt4/2*(v*(ft0+v)/(1+(2*pi*tt4*(ft0+v))**2) - v*(ft0-v)/(1+(2*pi*tt4*(ft0-v))**2))
	- (2*pi*tt5)**2*dt5/2*(v*(ft1+v)/(1+(2*pi*tt5*(ft1+v))**2) - v*(ft1-v)/(1+(2*pi*tt5*(ft1-v))**2))
	ei = 2*pi*v*(tt1*dt1/(1+(2*pi*v*tt1)**2) + tt2*dt2/(1+(2*pi*v*tt2)**2) + tt3*dt3/(1+(2*pi*v*tt3)**2))
	+ pi*v*tt4*dt4*(1/(1+(2*pi*tt4*(ft0+v))**2) + 1/(1+(2*pi*tt4*(ft0-v))**2))
	+ pi*v*tt5*dt5*(1/(1+(2*pi*tt5*(ft1+v))**2) + 1/(1+(2*pi*tt5*(ft1-v))**2))
	e =complex(er,-ei)
	m = cmath.sqrt(e)

	return m
