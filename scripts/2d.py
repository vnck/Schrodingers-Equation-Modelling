## Q1 --------------------------------------
import scipy.constants as c

def deg_to_rad(deg):
 rad = round(deg*c.pi/180, 5)
 return(rad)

def rad_to_deg(rad):
 deg = round(rad*180/c.pi, 5)
 return(deg)

## Q2 --------------------------------------
import numpy as np

def spherical_to_cartesian(r,theta,phi):
	x = r * np.cos(phi) * np.sin(theta)
	y = r * np.sin(phi) * np.sin(theta)
	z = r * np.cos(theta)
	return x, y, z

def cartesian_to_spherical(x, y, z):
	r = np.sqrt(x**2 + y**2 + z**2)
	if x!=0:
		phi=np.arctan(y/x)
	else:
		eps=1e-10
		phi=np.arctan(y/eps)
	if z == 0 and r == 0:
		theta = None
	else:
		theta = np.arccos(z/r)
	return r, theta, phi

## Q3 absolute / mag -----------------------
import numpy as np

def absolute(cnumber):
	 a = np.real(cnumber)
	 b = np.imag(cnumber)
	 mag = (a**2 + b**2)**.5
	 return mag

## Q4 --------------------------------------
import numpy as np
import scipy.constants as c

def angular_wave_func(m, l, theta, phi):
	if l == 0:
		if m == 0:
			y = np.sqrt(1/(4 * c.pi))
	elif l == 1:
		if m == 1:
			y = -1 * np.sqrt(3 / (8 * c.pi)) * np.sin(theta) * np.exp(phi * 1j)
		elif m == 0:
			y = np.sqrt(3 / (4 * c.pi)) * np.cos(theta)
		elif m == -1:
			y = np.sqrt(3 / (8 * c.pi)) * np.sin(theta) * np.exp(phi * -1j)
	elif l == 2:
		if m == 2:
			y1 = np.sqrt(15 / (32 * c.pi))
			y2 = np.sin(theta)**2
			y3 = np.exp(phi * 2j)
			y = y1 * y2 * y3
		elif m == 1:
			y1 = -1 * np.sqrt(15 / (18 * c.pi))
			y2 = np.sin(theta) * np.cos(theta)
			y3 = np.exp(phi * 1j)
			y = y1 * y2 * y3
		elif m == 0:
			y1 = np.sqrt(5 / (16 * c.pi))
			y2 = 3 * np.cos(theta)**2 - 1
			y3 = 1
			y = y1 * y2 * y3
		elif m == -1:
			y1 = np.sqrt(15 / (8 * c.pi))
			y2 = np.sin(theta) * np.cos(theta)
			y3 = np.exp(phi * -1j)
			y = y1 * y2 * y3
		elif m == -2:
			y1 = np.sqrt(15 / (32 * c.pi))
			y2 = np.sin(theta)**2
			y3 = np.exp(phi * -2j)
			y = y1 * y2 * y3
	elif l == 3:
		if m == 3:
			y1 = np.sqrt(35 / (64 * c.pi))
			y2 = np.sin(theta)**3
			y3 = np.exp(phi * 3j)
			y = -1 * y1 * y2 * y3
		elif m == 2:
			y1 = np.sqrt(105 / (32 * c.pi))
			y2 = np.sin(theta)**2 * np.cos(theta)
			y3 = np.exp(phi * 2j)
			y =  y1 * y2 * y3
		elif m == 1:
			y1 = np.sqrt(21 / (64 * c.pi))
			y2 = np.sin(theta) * (5 * np.cos(theta)**2 - 1)
			y3 = np.exp(phi * 1j)
			y =  -1 * y1 * y2 * y3
		elif m == 0:
			y1 = np.sqrt(7 / (16 * c.pi))
			y2 = 5 * np.cos(theta)**3 - 3 * np.cos(theta)
			y3 = 1
			y =  y1 * y2 * y3
		elif m == -1:
			y1 = np.sqrt(21 / (64 * c.pi))
			y2 = np.sin(theta) * (5 * np.cos(theta)**2 - 1)
			y3 = np.exp(phi * -1j)
			y =  y1 * y2 * y3
		elif m == -2:
			y1 = np.sqrt(105 / (32 * c.pi))
			y2 = np.cos(theta) * np.sin(theta)**2
			y3 = np.exp(phi * -2j)
			y =  y1 * y2 * y3
		elif m == -3:
			y1 = np.sqrt(35 / (64 * c.pi))
			y2 = np.sin(theta)**3
			y3 = np.exp(phi * -3j)
			y =  y1 * y2 * y3
	return np.round(y, 5)

## Q5 --------------------------------------
import numpy as np
import scipy.constants as c

a = c.physical_constants['Bohr radius'][0]

def radial_wave_func(n, l, r):
	if n == 1 and l == 0:
		y = 2 / np.sqrt(a**3) * np.exp(-r / a)
	elif n == 2:
		if l == 0:
			y = 1/2**.5 * a**-1.5 * (1-r/(2*a)) * np.exp(-r / (2*a))
		elif l == 1:
			y = 1/24**.5 * a**-1.5 * (r/a) * np.exp(-r / (2*a))
	elif n == 3:
		if l == 0:
			y = 2/(81*3**.5) * a**-1.5 * (27-18*(r/a)+2(r/a)**2) * np.exp(-r/(3*a))
		elif l == 1:
			y = 8/(27*6**.5) * a**-1.5 * (1 - (r/(6*a))) * (r/a) * np.exp(-r / (3*a))
		elif l == 2:
			y = 4/(81*30**.5) * a**-1.5 * (r/a)**2 * np.exp(-r/(3*a))
	elif n == 4:
		if l == 0:
			y = 1/4 * a**-1.5 * (1 - (3/4)*(r/a) + (1/8)*(r/a)**2 - (1/192)*(r/a)**3) * np.exp(-r/(4*a))
		elif l == 1:
			y = 5**.5/(16*3**.5) * a**-1.5 * (r/a) * (1 - (1/4)*(r/a) + (1/80)*(r/a)**2) * np.exp(-r/(4*a))
		elif l == 2:
			y = 1/(64*5**.5) * a**-1.5 * (r/a)**2 * (1 - (1/12)*(r/a)) * np.exp(-r/(4*a))
		elif l == 3:
			y = 1/(768*35**.5) * a**-1.5 * (r/a)**3 * np.exp(-r/(4*a))
	return np.round(y/a**-1.5, 5)

## Q6 --------------------------------------
def linspace(start, stop, num=50):
	ls = [round(start + x*(stop-start)/(num-1),5) for x in range(num)]
	return ls

## Q7 --------------------------------------
y = [4,5,6,7]
x = [1,2,3]
z = [8,9]

def meshgrid(x, y, z):
	lx = [[[y[i] for k in range(len(z))] for j in range(len(x))] for i in range(len(y))]
	ly = [[[x[j] for k in range(len(z))] for j in range(len(x))] for i in range(len(y))]
	lz = [[[z[k] for k in range(len(z))] for j in range(len(x))] for i in range(len(y))]

	return lx, ly, lz
	
## Q8

def norm_angular_wave_func(m,l,theta,phi):
	if m < 0:
		ylm = angular_wave_func(m, l, theta, phi)
		ylnm = angular_wave_func(-m, l, theta, phi)
		y_lm = j/2**.5 * (ylm - (-1)**m * ylnm)
	elif m == 0:
		y_lm = angular_wave_func(m, l, theta, phi)
	elif m > 0:
		ylm = angular_wave_func(m, l, theta, phi)
		ylnm = angular_wave_func(-m, l, theta, phi)
		y_lm = 1/2**.5 * (ylnm + (-1)**m * ylm)
	return y_lm


def hydrogen_wave_func(n,l,m,roa,nx,ny,nz):

	ls_x = linspace(-roa,roa,nx)
	ls_y = linspace(-roa,roa,ny)
	ls_z = linspace(-roa,roa,nz)

	ls_mesh = np.meshgrid(ls_x, ls_y,ls_z)
	
	r, theta, phi = np.vectorize(cartesian_to_spherical)(*ls_mesh)
	r *= a

	ang_mag = np.vectorize(norm_angular_wave_func)(m,l,theta,phi)
	rad_mag = np.vectorize(radial_wave_func)(n,l,r)

	density = absolute(ang_mag * rad_mag)**2

	density = np.round(np.double(density,5))

	return ls_mesh[0],ls_mesh[1],ls_mesh[2],density


# r,theta,phi = cartesian_to_spherical(x,y,z)

# if m < 0:
# 	ylm = angular_wave_func(m, l, theta, phi)
# 	ylnm = angular_wave_func(-m, l, theta, phi)
# 	y_lm = j/2**.5 * (ylm - (-1)**m * ylnm)
# elif m == 0:
# 	y_lm = angular_wave_func(m, l, theta, phi)
# elif m > 0:
# 	ylm = angular_wave_func(m, l, theta, phi)
# 	ylnm = angular_wave_func(-m, l, theta, phi)
# 	y_lm = 1/2**.5 * (ylnm + (-1)**m * ylm)

# rnl = radial_wave_func(n, l, r)

# psy = rnl * y_lm
# mag = absolute(psy) ** 2
# mag_list[b][c][coord] = mag

# x = mg[0][n][n]
# y = mg[1][n][n]
# x = mg[2][n][n]

# ([
#   [
#    [-8.0, -8.0, -8.0],
#    [0.0, 0.0, 0.0],
#    [8.0, 8.0, 8.0]
#   ],
#   [
#    [-8.0, -8.0, -8.0],
#    [0.0, 0.0, 0.0],
#    [8.0, 8.0, 8.0]
#   ],
#   [
#    [-8.0, -8.0, -8.0],
#    [0.0, 0.0, 0.0],
#    [8.0, 8.0, 8.0]
#   ]
#  ],
#  [
#   [
#    [-8.0, -8.0, -8.0],
#    [-8.0, -8.0, -8.0],
#    [-8.0, -8.0, -8.0]
#   ],
#   [
#    [0.0, 0.0, 0.0],
#    [0.0, 0.0, 0.0],
#    [0.0, 0.0, 0.0]
#   ],
#   [
#    [8.0, 8.0, 8.0],
#    [8.0, 8.0, 8.0],
#    [8.0, 8.0, 8.0]
#   ]
#  ],
#  [
#   [
#    [-8.0, 0.0, 8.0],
#    [-8.0, 0.0, 8.0],
#    [-8.0, 0.0, 8.0]
#   ],
#   [
#    [-8.0, 0.0, 8.0],
#    [-8.0, 0.0, 8.0],
#    [-8.0, 0.0, 8.0]
#   ],
#   [
#    [-8.0, 0.0, 8.0],
#    [-8.0, 0.0, 8.0],
#    [-8.0, 0.0, 8.0]
#   ]
#  ])

# [[[-8.0, -8.0, -8.0], [0.0, 0.0, 0.0], [8.0, 8.0, 8.0]], [[-8.0, -8.0, -8.0], [0.0, 0.0, 0.0], [8.0, 8.0, 8.0]], [[-8.0, -8.0, -8.0], [0.0
# , 0.0, 0.0], [8.0, 8.0, 8.0]]] [[[-8.0, -8.0, -8.0], [-8.0, -8.0, -8.0], [-8.0, -8.0, -8.0]], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0
# , 0.0]], [[8.0, 8.0, 8.0], [8.0, 8.0, 8.0], [8.0, 8.0, 8.0]]] [[[-8.0, 0.0, 8.0], [-8.0, 0.0, 8.0], [-8.0, 0.0, 8.0]], [[-8.0, 0.0, 8.0], 
# [-8.0, 0.0, 8.0], [-8.0, 0.0, 8.0]], [[-8.0, 0.0, 8.0], [-8.0, 0.0, 8.0], [-8.0, 0.0, 8.0]]] 

print(hydrogen_wave_func(2,1,1,8,3,3,3))
print('\n ######## \n')
print(hydrogen_wave_func(2,1,1,5,3,4,2))

# [012][012][0123][01]

# [
# 	[
# 		[-5.0, -5.0],
# 		[-5.0, -5.0],
# 		[-5.0, -5.0],
# 		[-5.0, -5.0]
# 	],
# 	[
# 		[0.0, 0.0],
# 		[0.0, 0.0],
# 		[0.0, 0.0],
# 		[0.0, 0.0]
# 	],
# 	[
# 		[5.0, 5.0],
# 		[5.0, 5.0],
# 		[5.0, 5.0],
# 		[5.0, 5.0]
# 	]
# ],
# [
# 	[
# 		[-5.0, -5.0],
# 		[-1.66667, -1.66667],
# 		[1.66667, 1.66667],
# 		[5.0, 5.0]
# 	], 
# 	[
# 		[-5.0, -5.0],
# 		[-1.66667, -1.66667],
# 		[1.66667, 1.66667],
# 		[5.0, 5.0]
# 	],
# 	[
# 		[-5.0, -5.0],
# 		[-1.66667, -1.66667],
# 		[1.66667, 1.66667],
# 		[5.0,5.0]
# 	]
# ],
# [
# 	[
# 		[-5.0, 5.0],
# 		[-5.0, 5.0],
# 		[-5.0, 5.0],
# 		[-5.0, 5.0]
# 	],
# 	[
# 		[-5.0, 5.0],
# 		[-5.0, 5.0], 
# 		[-5.0, 5.0],
# 		[-5.0, 5.0]
# 	],
# 	[
# 		[-5.0, 5.0],
# 		[-5.0, 5.0],
# 		[-5.0, 5.0],
# 		[-5.0, 5.0]
# 	]
# ]
