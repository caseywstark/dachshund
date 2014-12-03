from matplotlib import rcParams

rcParams["image.interpolation"] = "nearest"
rcParams["image.aspect"] = "equal"

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

l = 1.0
n = 11
gs = (n, n, n)

signal = np.fromfile("signal.bin").reshape(gs)
m = np.fromfile("map.bin").reshape(gs)

s_min = signal.min()
s_max = signal.max()

for ix in xrange(0, n):
    print "plotting slice ix = %i" % ix
    fig = plt.figure(figsize=(8, 5))
    gs = GridSpec(2, 2, height_ratios=(1, 30))

    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1])
    cax1 = fig.add_subplot(gs[0, 0])
    cax2 = fig.add_subplot(gs[0, 1])

    # signal
    im1 = ax1.imshow(signal[ix], vmin=s_min, vmax=s_max, extent=(0, l, 0, l))
    ax1.set_xlabel(r"$z$")

    cb1 = plt.colorbar(im1, cax=cax1, orientation="horizontal")
    cb1.solids.set_rasterized(True)
    cax1.tick_params(bottom=False, labelbottom=False, top=True, labeltop=True)
    cb1.set_label(r"$s$", labelpad=-40)

    # reconstruction
    im2 = ax2.imshow(m[ix], vmin=s_min, vmax=s_max, extent=(0, l, 0, l))
    ax2.set_xlabel(r"$z$")
    ax2.tick_params(labelleft=False)

    cb2 = plt.colorbar(im2, cax=cax2, orientation="horizontal")
    cb2.solids.set_rasterized(True)
    cax2.tick_params(bottom=False, labelbottom=False, top=True, labeltop=True)
    cb2.set_label(r"$\mathrm{reconstruction} \, s$", labelpad=-40)

    ax1.set_xlim(0, l)
    ax1.set_ylim(0, l)
    ax2.set_xlim(0, l)
    ax2.set_ylim(0, l)

    plt.tight_layout()
    plt.savefig("map_x%02d.pdf" % ix, dpi=100)
    plt.close()
