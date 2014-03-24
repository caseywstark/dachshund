import numpy as np

# domain size
l = 1.0
# gaussian center and size
c = 0.5
sigma = 0.3

# full map
map_n = 7
dx = l / map_n
ii = np.mgrid[0:map_n, 0:map_n, 0:map_n]
xx = dx * (ii + 0.5)
rr = xx - c
r2 = rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2]
g = np.exp(-r2/sigma**2)
g_mean = g.mean()
dg = g / g_mean - 1

dg.tofile("ideal_signal.bin")

# skewer setup

n = 4
num_skewers = n*n
num_pixels = n
pix_n = num_skewers * num_pixels

dx = l / n
dy = l / n
dz = l / num_pixels

# skewer pos.
sx = np.zeros(num_skewers)
sy = np.zeros(num_skewers)
d = np.zeros((num_skewers, num_pixels))

isk = 0
for ix in xrange(n):
    for iy in xrange(n):
        x = dx * (ix + 0.5)
        y = dx * (iy + 0.5)
        for iz in xrange(num_pixels):
            z = dz * (iz + 0.5)
            rx = x - c
            ry = y - c
            rz = z - c
            r2 = rx*rx + ry*ry + rz*rz
            d[isk, iz] = np.exp(-r2 / sigma**2)
        sx[isk] = x
        sy[isk] = y
        isk += 1

dd = d / g_mean - 1
w = np.ones(d.shape)

# write out
sx.tofile("skewer_x.bin")
sy.tofile("skewer_y.bin")
dd.tofile("pixel_data.bin")
w.tofile("pixel_weights.bin")

var_s = dd.std()**2
print "var_s = %g" % var_s

x_perp_2_list = []
x_para_2_list = []
s_list = []
ss_list = []

x2_list = []
cs_list = []

i = 0
for ixi in xrange(n):
    xi = dx * (ixi + 0.5)
    for iyi in xrange(n):
        yi = dy * (iyi + 0.5)
        for izi in xrange(num_pixels):
            zi = dz * (izi + 0.5)
            i = ixi * n + iyi
            di = dd[i, izi]

            for ixj in xrange(n):
                xj = dx * (ixj + 0.5)
                dx_ij = xj - xi

                for iyj in xrange(n):
                    yj = dy * (iyj + 0.5)
                    dy_ij = yj - yi
                    x_perp_2 = (dx_ij*dx_ij + dy_ij*dy_ij)

                    for izj in xrange(num_pixels):
                        zj = dz * (izj + 0.5)
                        dz_ij = zj - zi
                        x_para_2 = dz_ij*dz_ij

                        j = ixj * n + iyj
                        dj = dd[j, izj]
                        sij = di * dj

                        ssij = var_s * np.exp(-x_para_2/sigma) * np.exp(-x_perp_2/sigma)

                        x_perp_2_list.append(x_perp_2)
                        x_para_2_list.append(x_para_2)
                        s_list.append(sij)
                        ss_list.append(ssij)

                        x2 = dx_ij*dx_ij + dy_ij*dy_ij + dz_ij*dz_ij
                        x2_list.append(x2)

x_perp_2 = np.array(x_perp_2_list)
x_para_2 = np.array(x_para_2_list)
x2 = np.array(x2_list)
s = np.array(s_list)
ss = np.array(ss_list)
