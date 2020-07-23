"""In order to use the folowing functions you will need the functions VectorInt2D and VectorInt3D respectively from my repository titled vector interpolation"""

import numpy as np
import matplotlib.pyplot as plt


def VectorFieldTrace2D(fgrid, x0, y0, h, t0, tf):
    
    """Traces a vector field on a 2D grid where the value of the function is known at the vertices on the grid. Starts
    from a position (x0,y0) and traces the vector field for a given time period tf-t0, for a given step size h. 
    fgrid must be a 2D numpy array"""
    
    N = int((tf - t0) / h)

    t = np.linspace(t0, tf, N + 1)

    x = x0

    y = y0

    xn = [x0]

    yn = [y0]


    for i in range(N):
        k1 = VectorInt2D(fgrid, x, y)
        xn1 = x + 0.5 * h * k1
        yn1 = y + 0.5 * h * k1
        k2 = VectorInt2D(fgrid, xn1, yn1)
        xn2 = x + 0.5 * h * k2
        yn2 = y + 0.5 * h * k2
        k3 = VectorInt2D(fgrid, xn2, yn2)
        xn3 = x + h * k3
        yn3 = y + h * k3
        k4 = VectorInt2D(fgrid, xn3, yn3)
        x = x + (1 / 6) * h * (k1 + 2 * k2 + 2 * k3 + k4)
        y = y + (1 / 6) * h * (k1 + 2 * k2 + 2 * k3 + k4)
        xn.append(x)
        yn.append(y)
        u = xn[i] - xn[i-1]
        v = yn[i] - yn[i - 1]
    fig, ax = plt.subplots()
    q = ax.quiver(xn, yn, u, v, color='blue')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    return xn, yn

def VectorFieldTrace3D(fgrid, x0, y0, z0, h, t0, tf):
    """Traces a vector field on a 3D lattice where the value of the function is known at the vertices on the lattice. 
    Starts from a position (x0,y0,z0) and traces the vector field for a given time period tf-t0, for a given step size 
    h. fgrid must be a 3D numpy array"""
    
    N = int((tf - t0) / h)

    t = np.linspace(t0, tf, N + 1)

    x = x0

    y = y0

    z = z0

    xn = [x0]

    yn = [y0]

    zn = [z0]

    for i in range(N):
        k1 = VectorInt3D(fgrid, x, y, z)
        xn1 = x + 0.5 * h * k1
        yn1 = y + 0.5 * h * k1
        zn1 = z + 0.5 * h * k1
        k2 = VectorInt3D(fgrid, xn1, yn1, zn1)
        xn2 = x + 0.5 * h * k2
        yn2 = y + 0.5 * h * k2
        zn2 = z + 0.5 * h * k2
        k3 = VectorInt3D(fgrid, xn2, yn2, zn2)
        xn3 = x + h * k3
        yn3 = y + h * k3
        zn3 = z + h * k3
        k4 = VectorInt3D(fgrid, xn3, yn3, zn3)
        x = x + (1 / 6) * h * (k1 + 2 * k2 + 2 * k3 + k4)
        y = y + (1 / 6) * h * (k1 + 2 * k2 + 2 * k3 + k4)
        z = z + (1 / 6) * h * (k1 + 2 * k2 + 2 * k3 + k4)
        xn.append(x)
        yn.append(y)
        zn.append(y)
        u = xn[i] - xn[i - 1]
        v = yn[i] - yn[i - 1]
        k = zn[i] - zn[i - 1]
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.quiver(xn, yn, zn, u, v, k)
    plt.xlabel('x')
    plt.ylabel('y')
    ax.set_zlabel('z')
    return xn, yn, zn
