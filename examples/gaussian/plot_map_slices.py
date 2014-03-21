import numpy as np
import matplotlib.pyplot as plt

nx = 7
ny = 7
nz = 7
nsk = 16

m = np.fromfile("map.bin").reshape((nx, ny, nz))
sc = np.fromfile("skewer_coords.bin").reshape((2, nsk))

e = (0, 1, 0, 1)
vmin = m.min()
vmax = m.max()

for i in xrange(nz):
    print i
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)

    im = ax.imshow(m[:, :, i], vmin=vmin, vmax=vmax, extent=e, aspect="equal")
    p = ax.plot(sc[0], sc[1], ls="", marker=".", ms=5, c="black", alpha=0.5)
    cb = plt.colorbar(im)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")

    im.set_rasterized(True)
    cb.solids.set_rasterized(True)

    plt.tight_layout()
    plt.savefig("zslice_%i.pdf" % i, dpi=100)
    plt.close()
