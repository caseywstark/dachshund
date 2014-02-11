
import numpy as np
import matplotlib.pyplot as plt

n = 17
m = np.fromfile("map.bin").reshape((n, n, n))

e = (0, 1, 0, 1)
vmin = m.min()
vmax = m.max()

for i in xrange(n):
    print i
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)

    im = ax.imshow(m[:, :, i], vmin=vmin, vmax=vmax, extent=e, aspect="equal")
    cb = plt.colorbar(im)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")

    im.set_rasterized(True)
    cb.solids.set_rasterized(True)

    plt.tight_layout()
    plt.savefig("zslice_%i.pdf" % i, dpi=100)
    plt.close()
