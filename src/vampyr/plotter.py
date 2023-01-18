import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def plot_surface_xy(x=0.0, y=0.0, z=0.0, length=1.0):
    rx = [x, x + length]
    ry = [y, y + length]
    X, Y = np.meshgrid(rx, ry)
    r_1 = [z, z]
    r_2 = [z + length, z + length]
    p1, _ = np.meshgrid(r_1, r_1)
    p2, _ = np.meshgrid(r_2, r_2)
    return (X, Y, p1), (X, Y, p2)


def plot_surface_yz(x=0.0, y=0.0, z=0.0, length=1.0):
    ry = [y, y + length]
    rz = [z, z + length]
    Y, Z = np.meshgrid(ry, rz)
    r_1 = [x, x]
    r_2 = [x + length, x + length]
    p1, _ = np.meshgrid(r_1, r_1)
    p2, _ = np.meshgrid(r_2, r_2)
    return (p1, Y, Z), (p2, Y, Z)


def plot_surface_xz(x=0.0, y=0.0, z=0.0, length=1.0):
    rx = [x, x + length]
    rz = [z, z + length]
    X, Z = np.meshgrid(rx, rz)
    r_1 = [y, y]
    r_2 = [y + length, y + length]
    p1, _ = np.meshgrid(r_1, r_1)
    p2, _ = np.meshgrid(r_2, r_2)

    return (X, p1, Z), (X, p2, Z)


def plot_cube(corner, length):
    x, y, z = corner[0], corner[1], corner[2]
    a, b = plot_surface_xy(x, y, z, length)
    c, d = plot_surface_yz(x, y, z, length)
    e, f = plot_surface_xz(x, y, z, length)
    return (a, b, c, d, e, f)


def grid_plotter(tree=None, dpi=150, lw=0.08, color=(1, 0, 0, 0.01)):
    assert len(tree.MRA().world().upperBounds()) == 3, "basis plotter only works for 3D FunctionTrees"
    fig, ax = plt.subplots(figsize=(6, 6), dpi=150, subplot_kw={'projection': '3d'})
    ax.grid(False)
    ax.axis("off")

    corners = [tree.fetchEndNode(i).lowerBounds() for i in range(tree.nEndNodes())]
    lengths = [
        tree.fetchEndNode(i).upperBounds()[0] - tree.fetchEndNode(i).lowerBounds()[0]
        for i in range(tree.nEndNodes())
    ]

    for i in range(tree.nEndNodes()):
        data = plot_cube(corners[i], lengths[i])
        for d in data:
            ax.plot_surface(d[0], d[1], d[2], color=color, edgecolor='black', lw=lw)

    return fig, ax


def representation_vs_basis(tree, type="scaling"):
    assert len(tree.MRA().world().upperBounds()) == 1, "basis plotter only works for 1D FunctionTrees"
    mra = tree.MRA()
    k = mra.basis().scalingOrder()
    n = tree.depth() - 1

    upper_bound = mra.world().upperBound(0)
    lower_bound = mra.world().lowerBound(0)
    x = np.arange(lower_bound, upper_bound, 0.001)
    y = [tree([_]) for _ in x] # Plot f1 to f4 to see how the function representation improves.

    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    ax[1].title.set_text(f"Basis")
    ax[0].plot(x, y)
    ax[0].set_xlim(lower_bound, upper_bound)
    ax[1].set_xlim(lower_bound, upper_bound)
    basis_polys = None
    if type.lower() == "scaling":
        ax[0].title.set_text(f"Function tree: $f^{n}(x)$")
        basis_polys = mra.basis().scaling
    elif type.lower() == "wavelet":
        ax[0].title.set_text(f"Function tree: $df^{n}(x)$")
        basis_polys = mra.basis().wavelet
    else:
        raise Exception(f"{type} : type not recogniced")
    for _ in range(tree.nEndNodes()):
        idx = tree.fetchEndNode(_).index()
        n = idx.scale()
        l = idx.translation()[0]

        for i in range(k+1):
            y = [basis_polys(i=i, l=l, n=n)([x]) if basis_polys(i=i, l=l, n=n)([x]) != 0.0 else np.nan for x in x]
            ax[1].plot(x, y)

    return fig, ax


