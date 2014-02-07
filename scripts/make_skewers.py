import subprocess
import h5py
import numpy as np
from scipy.signal import decimate

#snap_path = "l10_n128/z25.h5"
snap_path = "/Users/caseywstark/remote/alcc/Analysis/convergence/runs/l10_n128/z25.h5"
f = h5py.File(snap_path, "r")

# sim params...
L = f["domain"].attrs["size"][0]
N = f["domain"].attrs["shape"][0]
dx = L / N
dy = dx

# skewer params
num_skewer_1d = 4
skewer_step = N / num_skewer_1d
skewer_indexes = np.arange(skewer_step / 2, N, skewer_step)
skewer_downsample = 4

print "Reading tau dataset."

tau_ds = f["derived_fields/tau_red"]
tau = tau_ds.value.astype(np.float64)
f = np.exp(-tau)
mean_flux = f.mean()

print "Making skewers."
skewer_x = []
skewer_y = []
delta_f = []
for ix in skewer_indexes:
    for iy in skewer_indexes:
        tau = tau_ds[ix, iy].astype(np.float64)
        tau = decimate(tau, skewer_downsample)

        f = np.exp(-tau)
        df = f / mean_flux - 1.0

        skewer_x.append(dx * (ix + 0.5))
        skewer_y.append(dy * (iy + 0.5))
        delta_f.append(df)

# array cast's
skewer_x = np.array(skewer_x)
skewer_y = np.array(skewer_y)
delta_f = np.array(delta_f)
noise = 0.05 * np.ones(delta_f.shape)

print "Writing skewer file."

# write out individual files because python's binary writing is terrible.
skewer_x.tofile("skewer_x")
skewer_y.tofile("skewer_y")
delta_f.tofile("delta_f")
noise.tofile("noise")

# cat files together
# grumble grumble subprocess
c = "cat skewer_x skewer_y delta_f noise > skewers.bin"
process = subprocess.call(c, shell=True)

# clean up after one's self...
c = "rm skewer_x skewer_y delta_f noise"
process = subprocess.call(c, shell=True)
