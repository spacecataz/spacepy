#!/usr/bin/env python
'''
Some functions for the generation of a dipole field.

Copyright 2010 Los Alamos National Security, LLC.
'''

import numpy as np
import math
import pylab as plt

def b_mag(x,y,z):
    '''
    For a position *x*, *y*, *z* in units of planetary radius, return the
    strength of the dipole magnetic field in nanoTesla.
    
    Input and output values are in SM coordinates.
    '''
    r = np.sqrt(x**2 + y**2 + z**2)
    cos = z/r

    return 30400 * np.sqrt(1+3*cos**2)/r**3
    
def b_hat(x, y, z):
    ''' 
    For locations at X, Y, and Z in SM coordinates, calculate the direction
    of the magnetic field in X, Y, and Z about a grid for all X-Y-Z
    combinations.  Three arrays of size `x.size,y.size,z.size` are returned.
    These arrays are useful for tracing field lines, plotting quiver plots,
    etc.

    Note that the returned matrices are generated using Numpy's *meshgrid*
    function using 'ij' indexing.  For quiver plotting or other 3D handling,
    location meshes should be constructed using the same 'ij' indexing.

    If X, Y, and Z are not the same size, any single-entry dimensions
    will be pruned.

    Examples
    ========
    # Create 2D fields in the Y=0 plane:
    >>> import numpy as np
    >>> from spacepy.pybats.dipole import b_hat

    >>> x = np.arange(-100.0, 101.0, 5.0)
    >>> y = 0
    >>> z = np.arange(-100.0, 101.0, 5.0)

    >>> xgrid, ygrid, zgrid = np.meshgrid(x,y,z)
    >>> bx, by, bz = b_hat(x,y,z)
    >>> bx.shape

    (41, 41)

    # Create 3D fields, plot quivers on 3D axes:
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from spacepy.pybats.dipole import b_hat

    >>> x = np.arange(-10, 11, 2.0)
    >>> y = np.arange(  1, 11, 2.0)
    >>> z = np.arange(-10, 11, 2.0)

    >>> xgrid, ygrid, zgrid = np.meshgrid(x,y,z, indexing='ij')
    >>> bx, by, bz = b_hat(x,y,z)

    >>> ax = plt.subplot(111, projection='3d')
    >>> ax.quiver(xgrid,ygrid,zgrid, bx,by,bz, length=4, alpha=.2)

    '''

    # Create meshgrids of inputs:
    xgrid, ygrid, zgrid = np.meshgrid(x,y,z, indexing='ij')

    # Calculate geometrical factors:
    r = np.sqrt(xgrid**2 + ygrid**2 + zgrid**2) # radial distance
    r_xy = np.sqrt(xgrid**2 + ygrid**2)         # radial distance from Z
    phi  = np.arctan2(ygrid, xgrid)             # Azimuth
    sin, cos = zgrid/r, r_xy/r                  # Cosine(colat)

    # Normalization factor: b_hat = B/|B|:
    denom = 1/np.sqrt(1.0 + 3.0*cos**2)

    b_x = -denom*cos*sin*np.cos(phi)
    b_y = -denom*cos*sin*np.sin(phi)
    b_z = denom*(cos**2-1/3)
        
    return(b_x.squeeze(), b_y.squeeze(), b_z.squeeze())

def b_line(x, y, z, npoints=30):
    '''
    For a starting X, Y, Z point in SM coordinates, return x, y, and z
    vectors that trace the dipole field line that passes through the given 
    point.
    ''' 
    npoints = float(npoints)

    # Radii of interest:
    r = np.sqrt(x**2+y**2+z**2) # Radial distance to starting point.
    r_xy = np.sqrt(x**2 + y**2) # Radial distance in equatorial plane.
    theta = np.arctan2(r_xy,z)  # Colatitude.
    phi = np.arctan2(y,x)       # Local time angle from Noon.
    R = r/(np.sin(theta)**2)    # Maximum R along field at equator.

    #if x<0:
    #    theta = np.arange(np.pi, 2.0*np.pi, np.pi/npoints)
    #else:
    theta = np.arange(0, np.pi, np.pi/npoints)
    r_vec = R * np.sin(theta)**2

    xy = r_vec * np.sin(theta)
    x_out = xy * np.cos(phi)
    y_out = xy * np.sin(phi)
    z_out = r_vec * np.cos(theta)

    return (x_out, y_out, z_out)

def test2d():
    '''
    A quick test of the dipole field functions.
    '''

    # 2D fields:
    x = np.arange(-100.0, 101.0, 5.0)
    y = 0
    z = np.arange(-100.0, 101.0, 5.0)

    xgrid, ygrid, zgrid = np.meshgrid(x,y,z)#, indexing='ij')
    x_vec, y_vec, z_vec = b_hat(x,y,z)

    fig = plt.figure(figsize=(10,8))
    ax1 = plt.subplot(111)

    #ax1.quiver(xgrid[:,0,:], zgrid[:,0,:], x_vec[:,0,:], z_vec[:,0,:])
    ax1.quiver(xgrid, zgrid, x_vec, z_vec)

    for i in range(-120, 121, 10):
        (x,y,z) = b_line(float(i), 0.0, 0.0, 100)
        ax1.plot(x, z, 'b')
    for theta in np.arange(np.pi/2.0, 3.0*np.pi/2.0, np.pi/100.0):
        x = np.sin(theta)
        z = np.cos(theta)
        x,y,z = b_line(x, 0, z, 100)
        ax1.plot(x, z,'r')
    ax1.set_xlim( [-100, 100] )
    ax1.set_ylim( [-100, 100] )
    plt.title('Unit vectors for an arbitrary dipole field')
    
    fig.show()

def test3d():
    from spacepy.plot import set_aspect3d

    ax = plt.subplot(111, projection='3d')
    for theta in np.arange(0, 2*np.pi, np.pi/8.0):
        x = np.cos(theta)
        y = np.sin(theta)
        for r, c in zip( [3,6,9], ['b','r','g']):
            bx,by,bz = b_line(r*x, r*y, 0, 100)
            ax.plot(bx, by, bz, c, alpha=.5)

    size = 10
    x = np.arange(-size, size+1, 2.0)
    y = np.arange(-size, size+1, 2.0)
    z = np.arange(-size, size+1, 2.0)
    xgrid, ygrid, zgrid = np.meshgrid(x,y,z, indexing='ij')
    bx, by, bz = b_hat(x,y,z)
    ax.quiver(xgrid,ygrid,zgrid, bx,by,bz, length=4, alpha=.2)
    set_aspect3d(ax)
    ax.set_xlabel("X", size=20)
    ax.set_ylabel("Y", size=20)
    ax.set_zlabel("Z", size=20)
    
if __name__ == '__main__':
    test2d()
    test3d()
