import numpy as np

l = 1.0
n = 4

num_skewers = n*n
num_pixels = n
pix_n = num_skewers * num_pixels

dx = l / n

# skewer pos.
sx = np.zeros([num_skewers])
sy = np.zeros([num_skewers])

i = 0
for ix in xrange(n):
    for iy in xrange(n):
        sx[i] = dx * (ix + 0.5)
        sy[i] = dx * (iy + 0.5)
        i += 1

# gaussian center
c = 0.5

# set data and weights.
ii = np.mgrid[0:n, 0:n, 0:n]
xx = dx * (ii + 0.5)
rr = xx - c
r2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2]
g = np.exp(-r2/0.2**2)

w = np.ones(g.shape)

# write out
sx.tofile("skewer_x.bin")
sy.tofile("skewer_y.bin")
g.tofile("pixel_data.bin")
w.tofile("pixel_weights.bin")

