import numpy as np
from mayavi import mlab

x = np.load('xdata432.dat')
y = np.load('ydata432.dat')
z = np.load('zdata432.dat')
density = np.load('magdata432.dat')

figure = mlab.figure('DensityPlot')

# you should modify the parameters
pts = mlab.contour3d(density ,contours =20, opacity =0.5)
mlab.axes()
mlab.show()
