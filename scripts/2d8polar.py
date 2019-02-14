from sympy import diff, exp, factorial, symbols, cos, simplify
from math import e, pi
from numpy import arctan2, arccos, sqrt, vectorize, array, meshgrid, absolute, real, imag, warnings
import scipy.constants as c
from functools import wraps

a = c.physical_constants["Bohr radius"][0]


def cartesian_to_spherical(x, y, z):
    r = round((x**2+y**2+z**2)**0.5,5)
    theta = arctan2(y,x)
    z = absolute(z)
    if r==0:
        phi = arccos(0)
    else:
        phi = arccos(z/r)
    return (r, phi, theta)

def return_cartesian(x,y,z):
    return x, y, z

## Q3 absolute / mag -----------------------

def absolute(cnumber):
     a = real(cnumber)
     b = imag(cnumber)
     mag = (a**2 + b**2)**.5
     return mag

## Q4 --------------------------------------

def assoc_legendre(mv,lv):
    m, l, x = symbols('m l x')

    def P_l():
        t1 = 1 / (2**l * factorial(l))
        t2 = diff( (x**2 - 1)**l, x, int(lv) )

        return simplify(t1 * t2)

    def legendre(theta):
        xv = cos(theta)

        t1 = (1 - x**2)**(abs(m)/2)
        t2 = diff( P_l(), x, int(abs(mv)))
        expression = simplify(t1*t2)

        ans = expression.subs({m:mv, l:lv, x:xv})
        return ans

    return legendre

def assoc_laguerre(pv,qmpv):

    qv = qmpv + pv
    p, q, x = symbols('p q x')

    def L_q():
        t1 = exp(x)
        t2exp = exp(-x) * x**qv
        t2 = diff( t2exp, x, int(qv))

        return simplify(t1*t2)

    def laguerre(xv):
        t1 = (-1)**p
        t2 = diff( L_q(), x, int(pv))

        expression = simplify(t1*t2)

        ans = expression.subs({p:pv, q:qv, x:xv})
        return ans

    return laguerre

def angular_wave_func(m,l,theta,phi):
    pml = assoc_legendre(m,l)
    
    epsilon = (-1)**m if m > 0 else 1

    t1 = (2*l + 1) / (4*pi)
    t2 = factorial(l - abs(m)) / factorial(l + abs(m))
    t12 = (t1 * t2)**0.5

    t3 = e**(1j*m*phi) * pml(theta)

    ans = complex(epsilon * t12 * t3, 0j)

    return ans

def radial_wave_func(n, l, r):
    p = 2*l + 1
    qmp = n - l - 1
    Lpqmp = assoc_laguerre(p, qmp)

    t1 = (2 / (n*a))**3
    t2 = factorial(qmp) / (2*n*(factorial(n+l)**3))

    t12 = (t1 * t2)**0.5
    t3 = e**(-r / (n*a))
    t4 = (2*r / (n*a))**l
    t5 = Lpqmp(2*r / (n*a))

    ans = (t12*t3*t4*t5)/a**-1.5

    return ans

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

    ls_x = linspace(-roa,roa,nx)
    ls_y = linspace(-roa,roa,ny)
    ls_z = linspace(-roa,roa,nz)

    ls_mesh = meshgrid(ls_y,ls_x,ls_z)
    
    r,theta,phi = vectorize(cartesian_to_spherical)(*ls_mesh)
    r *= a

    ang_mag = vectorize(norm_angular_wave_func)(m,l,theta,phi)
    rad_mag = vectorize(radial_wave_func)(n,l,r)

    density = vectorize(absolute)(ang_mag * rad_mag)

    density = vectorize(round)(density**2,5)

    return array(ls_mesh[0]),array(ls_mesh[1]),array(ls_mesh[2]),density

file = open("printout2.txt","w")

file.write(str(hydrogen_wave_func(2,1,1,8,3, 3, 3)))                                                                                      
file.write(str(hydrogen_wave_func(2,1,1,5,3,4,2)))                                                                                        
file.write(str(hydrogen_wave_func(2,0,0,3,5,4,3)))                                                                                        
file.write(str(hydrogen_wave_func(3,1,0,10,6,6,6)))           

file.close()
