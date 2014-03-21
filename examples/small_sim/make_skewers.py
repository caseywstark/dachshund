import h5py
import numpy as np

#snap_path = "l10_n128/z25.h5"
snap_path = "/project/projectdirs/m1796/Analysis/convergence/runs/l40_n2048/z25.h5"
f = h5py.File(snap_path, "r")

# sim params...
L = f["domain"].attrs["size"][0]
N = f["domain"].attrs["shape"][0]
dx = L / N
dy = dx

tau_ds = f["derived_fields/tau_red"]

# skewer params
num_skewers = 128
num_pixels = 64
skewer_coords = np.random.randint(N, size=(num_skewers, 2))

dz = L / N
dz_pix = L / num_pixels
sigma_z = dz_pix / 2
z_cell = np.linspace(dz/2, L-dz/2, num=N)

print "Making skewers."
skewer_x = []
skewer_y = []
data = []
for ix, iy in skewer_coords:
    print ix, iy

    skewer_x.append(dx * (ix + 0.5))
    skewer_y.append(dy * (iy + 0.5))

    tau = tau_ds[ix, iy].astype(np.float64)
    flux = np.exp(-tau)

    for ipix in xrange(num_pixels):
        z = dz_pix * (ipix + 0.5)
        z_sep = z - z_cell
        w = np.exp(-(z_sep/sigma_z)**2)
        # w.sum should be 1, but why not
        fw = (w * flux).sum() / w.sum()
        data.append(fw)

# array cast's
skewer_x = np.array(skewer_x)
skewer_y = np.array(skewer_y)
data = np.array(data)

# random weights
weights = 0.2 + 0.8 * np.random.random(data.shape)

print "Writing input files."
skewer_x.tofile("skewer_x.bin")
skewer_y.tofile("skewer_y.bin")
data.tofile("pixel_data.bin")
weights.tofile("pixel_weights.bin")

