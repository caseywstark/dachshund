import numpy as np
import matplotlib.pyplot as plt

lx = 1.0
ly = 1.0
lz = 1.0
nx = 7
ny = 7
nz = 7

sx = np.fromfile("skewer_x.bin")
sy = np.fromfile("skewer_y.bin")
ideal_signal = np.fromfile("ideal_signal.bin").reshape(nx, ny, nz)
m = np.fromfile("map.bin").reshape(nx, ny, nz)
mcov_nd = np.fromfile("map_covar_nd.bin").reshape((nx, ny, nz))

s_min = ideal_signal.min()
s_max = ideal_signal.max()
mc_min = mcov_nd.min()
mc_max = mcov_nd.max()

# Figure params
fig_w = 8.0
fig_h = 3.0
aspect = fig_w / fig_h

l = 0.08
b = 0.14
im_w = 0.22
im_h = aspect * im_w
im_wpad = 0.0
im_hpad = 0.0
cb_w = im_w
cb_h = 0.02
cb_wpad = (im_w - cb_w) / 2
lp = -37

for ix in xrange(0, nx):
    print "plotting slice ix = %i" % ix
    fig = plt.figure(figsize=(fig_w, fig_h))

    ax1 = fig.add_axes((l, b, im_w, im_h))
    ax2 = fig.add_axes((l + im_w + im_wpad, b, im_w, im_h))
    ax3 = fig.add_axes((l + 2*(im_w + im_wpad), b, im_w, im_h))
    ax4 = fig.add_axes((l + 3*(im_w + im_wpad), b, im_w, im_h))

    cax1 = fig.add_axes((l + cb_wpad, b + im_h + im_hpad, cb_w, cb_h))
    cax2 = fig.add_axes((l + im_w + im_wpad + cb_wpad, b + im_h + im_hpad, cb_w, cb_h))
    cax3 = fig.add_axes((l + 2*(im_w + im_wpad) + cb_wpad, b + im_h + im_hpad, cb_w, cb_h))
    cax4 = fig.add_axes((l + 3*(im_w + im_wpad) + cb_wpad, b + im_h + im_hpad, cb_w, cb_h))

    # skewer position plot
    xx = lx / nx * (ix + 0.5)
    ax1.axvline(xx, ls="-", lw=1, c="k", alpha=0.5)
    s = ax1.scatter(sx, sy, c=np.ones(sx.size), marker="s", s=5, linewidths=0.0, alpha=0.7)
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$y$")

    cb1 = plt.colorbar(s, cax=cax1, orientation="horizontal")
    cb1.solids.set_rasterized(True)
    cax1.tick_params(bottom=False, labelbottom=False, top=True, labeltop=True)
    cb1.set_label(r"$\rm skewer \, pos/weights$", labelpad=lp)

    # simulation delta f
    im2 = ax2.imshow(ideal_signal[ix], vmin=s_min, vmax=s_max,
        extent=(0, lz, 0, ly), aspect="equal")
    im2.set_rasterized(True)
    ax2.set_xlabel(r"$z$")
    ax2.tick_params(labelleft=False)

    cb2 = plt.colorbar(im2, cax=cax2, orientation="horizontal")
    cb2.solids.set_rasterized(True)
    cax2.tick_params(bottom=False, labelbottom=False, top=True, labeltop=True)
    cb2.set_label(r"$\mathrm{ideal} \, s$", labelpad=lp)

    # reconstruction
    im3 = ax3.imshow(m[ix], vmin=s_min, vmax=s_max, extent=(0, lz, 0, ly),
        aspect="equal")
    im3.set_rasterized(True)
    ax3.set_xlabel(r"$z$")
    ax3.tick_params(labelleft=False)

    cb3 = plt.colorbar(im3, cax=cax3, orientation="horizontal")
    cb3.solids.set_rasterized(True)
    cax3.tick_params(bottom=False, labelbottom=False, top=True, labeltop=True)
    cb3.set_label(r"$\mathrm{reconstruction} \, s$", labelpad=lp)

    # covar plot
    im4 = ax4.imshow(mcov_nd[ix], vmin=mc_min, vmax=mc_max,
        extent=(0, lz, 0, ly), aspect="equal")
    im4.set_rasterized(True)
    ax4.set_xlabel(r"$z$")
    ax4.tick_params(labelleft=False)

    cb4 = plt.colorbar(im4, cax=cax4, orientation="horizontal")
    cb4.solids.set_rasterized(True)
    cb4.set_label(r"$\lim_{N \gg S} M_{ii}$", labelpad=lp)
    cax4.tick_params(bottom=False, labelbottom=False, top=True, labeltop=True)

    ax1.set_xlim(0, lx)
    ax1.set_ylim(0, ly)
    ax2.set_xlim(0, lz)
    ax2.set_ylim(0, ly)
    ax3.set_xlim(0, lz)
    ax3.set_ylim(0, ly)
    ax4.set_xlim(0, lz)
    ax4.set_ylim(0, ly)

    #ax2.set_xticks([20, 40, 60, 80, 100])
    #ax3.set_xticks([20, 40, 60, 80, 100])
    #ax4.set_xticks([20, 40, 60, 80, 100])
    #cb1.set_ticks([0, 250, 500])
    #cb2.set_ticks([-0.8, -0.4, 0])
    #cb3.set_ticks([-0.8, -0.4, 0])
    #cb4.set_ticks([0, 1, 2])

    plt.savefig("map_x%02d.pdf" % ix, dpi=100)
    plt.close()
