"""
Generate input files for a problem where the signal is

s(x) = a exp[-(x - x0)^2 / b^2]

"""

import numpy as np

# gaussian signal
def g_signal(x, x0, a, b):
    return a * np.exp(-( (x - x0)**2 ).sum(axis=1) / b**2)

# domain size
l = 1.0
map_n = 11

# pixels
num_pixels = 500
pix_x = l * np.random.random((num_pixels, 3))

# gaussian params
a = 1.0
x0 = np.array([0.5, 0.5, 0.5])
b = 0.15

# signal
pix_g = g_signal(pix_x, x0, a, b)
mean_g = 0.5
pix_s = pix_g / mean_g - 1.0
# noise
pix_nstd = 0.1 * a
pix_n = pix_nstd * np.ones(num_pixels)
# just take d = s
pix_d = pix_s

# pixel data altogether
p = np.vstack([pix_x.T[0], pix_x.T[1], pix_x.T[2], pix_n, pix_d]).T
p.tofile("pixel_data.bin")

# write config file
cf = open("run.cfg", "w")
cf.write("lx = %f\n" % l)
cf.write("ly = %f\n" % l)
cf.write("lz = %f\n" % l)
cf.write("num_pixels = %i\n" % num_pixels)
cf.write("map_nx = %i\n" % map_n)
cf.write("map_ny = %i\n" % map_n)
cf.write("map_nz = %i\n" % map_n)
cf.write("corr_var_s = 1.0\n")
cf.write("corr_l_perp = %f\n" % b)
cf.write("corr_l_para = %f\n" % b)
cf.write("pcg_tol = 1.0e-5\n")
cf.write("pcg_max_iter = 1000\n")

# full map
dx = l / map_n
ii = np.mgrid[0:map_n, 0:map_n, 0:map_n]
xx = dx * (ii + 0.5)
xx = xx.reshape(3, map_n**3).T
g = g_signal(xx, x0, a, b)
s = g / mean_g - 1.0
s.tofile("signal.bin")
