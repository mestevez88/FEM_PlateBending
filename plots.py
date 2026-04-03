import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def Plot_PlateNodes(mesh_x, mesh_y, EX, EY, LX, LY, Z, caption):
    X = np.reshape(mesh_x, (EX + 1, EY + 1))
    Y = np.reshape(mesh_y, (EX + 1, EY + 1))
    scale = .25
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    # Plot the surface with deflection coordinates negative
    dplot = np.reshape(Z, ((EX + 1)*(EY + 1), 3))
    dplot = dplot[:,0]
    dplot = np.reshape(dplot, (EX + 1, EY + 1))
    surf = ax.plot_surface(X, Y, -dplot, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    # Customize the z axis.
    ax.set_box_aspect([1, LY/LX, scale])
    # ax.zaxis.set_major_locator(LinearLocator(10))
    # ax.zaxis.set_major_formatter('{x:.02f}')
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.25, aspect=5)
    plt.title(caption)
    plt.show()
    
def Plot_PlateCenter(mesh_x, mesh_y, EX, EY, LX, LY, Z, caption):
    'Ploy values on the center of each plate, normally comming from averaging values on Gauss Points'
    X = np.reshape(mesh_x, (EX + 1, EY + 1))
    Y = np.reshape(mesh_y, (EX + 1, EY + 1))
    X = X[0:EX, 0:EY]
    Y = Y[0:EX, 0:EY]
    "Translate coordinates to plate center"
    X = X + .5*(LX/EX)
    Y = Y + .5*(LY/EY)
    Z = np.reshape(Z, (EX, EY))
    scale = .25
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    # Plot the surface with deflection coordinates negative
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    # Customize the z axis.
    ax.set_box_aspect([1, LY/LX, scale])
    # ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter('{x:.02f}')
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.25, aspect=5)
    plt.title(caption)
    plt.show()