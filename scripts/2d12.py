import numpy as np
import scipy.constants as c

a = c.physical_constants["Bohr radius"][0]

def cartesian_to_spherical(x, y, z):
	r = round((x**2+y**2+z**2)**0.5,5)
	theta = np.arctan2(y,x)
	z = absolute(z)
	if r==0:
		phi = np.arccos(0)
	else:
		phi = np.arccos(z/r)
	return (r, phi, theta)

def return_cartesian(x,y,z):
	return x, y, z

## Q3 absolute / mag -----------------------

def absolute(cnumber):
	 a = np.real(cnumber)
	 b = np.imag(cnumber)
	 mag = (a**2 + b**2)**.5
	 return mag

## Q4 --------------------------------------

def angular_wave_func(m,l,x,y,z):
	r = (x**2+y**2+z**2)**0.5
	try:
		if l==0:
			y_a = (1/(4*np.pi))**0.5
		elif l==1:
			if m==-1:
				y_a = (y/r)*(3/(4*np.pi))**0.5
			elif m==0:
				y_a = (z/r)*(3/(4*np.pi))**0.5
			elif m==1:
				y_a = (x/r)*(3/(4*np.pi))**0.5
		elif l==2:
			if m==-2:
				y_a = (x*y/r**2)*(15/(4*np.pi))**0.5
			elif m==-1:
				y_a = (y*z/r**2)*(15/(4*np.pi))**0.5
			elif m==0:
				y_a = ((-x**2-y**2+2*z**2)/r**2)*(5/(16*np.pi))**0.5
			elif m==1:
				y_a = (z*x/r**2)*(15/(4*np.pi))**0.5
			elif m==2:
				y_a = ((x**2-y**2)/r**2)*(15/(16*np.pi))**0.5
		elif l==3:
			if m==-3:
				y_a = (((3*x**2-y**2)*y)/r**3)*(35/(32*np.pi))**0.5
			elif m==-2:
				y_a = ((x*y*z)/r**3)*(105/(4*np.pi))**0.5
			elif m==-1:
				y_a = ((y*(4*z**2-x**2-y**2))/r**3)*(21/(32*np.pi))**0.5
			elif m==0:
				y_a = ((z*(2*z**2-3*x**2-3*y**2))/r**3)*(7/(16*np.pi))**0.5
			elif m==1:
				y_a = ((x*(4*z**2-x**2-y**2))/r**3)*(21/(32*np.pi))**0.5
			elif m==2:
				y_a = (((x**2-y**2)*z)/r**3)*(105/(16*np.pi))**0.5
			elif m==3:
				y_a = (((x**2-3*y**2)*x)/r**3)*(35/(32*np.pi))**0.5
		else:
			return None      
		return y_a
	except ZeroDivisionError:
		return None

def radial_wave_func(n,l,r):
	if n==1:
		r_w = (2/(a**1.5))*np.exp(-r/a)
	elif n==2:
		if l==0:
			r_w = (1/(2**0.5))*(a**-1.5)*(1-r/(2*a))*np.exp(-r/(2*a))
		elif l==1:
			r_w = (1/(24**0.5))*(a**-1.5)*(r/a)*np.exp(-r/(2*a))
	elif n==3:
		if l==0:
			r_w = (2/(81*3**0.5))*(a**-1.5)*(27-18*(r/a)+2*(r/a)**2)*np.exp(-r/(3*a))
		elif l==1:
			r_w = (8/(27*6**0.5))*(a**-1.5)*(1-r/(6*a))*(r/a)*np.exp(-r/(3*a))
		elif l==2:
			r_w = (4/(81*30**0.5))*(a**-1.5)*(r/a)**2*np.exp(-r/(3*a))
	elif n==4:
		if l==0:
			r_w = 0.25*(a**-1.5)*(1-0.75*(r/a)+(1/8)*(r/a)**2-(1/192)*(r/a)**3)*np.exp(-r/(4*a))
		elif l==1:
			r_w = ((5**0.5)/(16*3**0.5))*(a**-1.5)*(r/a)*(1-0.25*(r/a)+(1/80)*(r/a)**2)*np.exp(-r/(4*a))
		elif l==2:
			r_w = (1/(64*5**0.5))*(a**-1.5)*(r/a)**2*(1-(1/12)*(r/a))*np.exp(-r/(4*a))
		elif l==3:
			r_w = (1/(768*35**0.5))*(a**-1.5)*(r/a)**3*np.exp(-r/(4*a))
	else:
		return None
	r_w=r_w/(a**-1.5)
	return r_w

## Q6 --------------------------------------
def linspace(start, stop, num=50):
	ls = [round(start + x*(stop-start)/(num-1),5) for x in range(num)]
	return ls

## Q7 --------------------------------------

def meshgrid(x, y, z):
	lx = [[[y[i] for k in range(len(z))] for j in range(len(x))] for i in range(len(y))]
	ly = [[[x[j] for k in range(len(z))] for j in range(len(x))] for i in range(len(y))]
	lz = [[[z[k] for k in range(len(z))] for j in range(len(x))] for i in range(len(y))]

	return lx, ly, lz
	
## Q8

def hydrogen_wave_func(n,l,m,roa,nx,ny,nz):

	ls_mesh = np.mgrid[-1*roa:roa:(nx*1j),-1*roa:roa:(ny*1j),-1*roa:roa:(nz*1j)]
	
	x, y, z = np.vectorize(return_cartesian)(*ls_mesh)
	r,phi,theta = np.vectorize(cartesian_to_spherical)(*ls_mesh)

	ang_mag = np.vectorize(angular_wave_func)(m,l,x, y, z)

	r *= a
	rad_mag = np.vectorize(radial_wave_func)(n,l,r)

	density = np.vectorize(np.absolute)(ang_mag * rad_mag)

	density = np.vectorize(round)(density**2,5)

	return ls_mesh[0],ls_mesh[1],ls_mesh[2],density

x, y, z, mag = hydrogen_wave_func(4, 0, 0, 30, 100, 100, 100)

print('Writing to x......')
x.dump('xdata432.dat')
print('Writing to y......')
y.dump('ydata432.dat')
print('Writing to z......')
z.dump('zdata432.dat')
print('Writing to mag......')
mag.dump('magdata432.dat')

# print(x)
# print(y)
# print(z)
# print(mag)
