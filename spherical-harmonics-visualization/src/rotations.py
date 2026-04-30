import numpy as np
from numpy.f2py.auxfuncs import throw_error
from scipy.linalg import block_diag
from math import *
from scipy.linalg import expm

# wigner matrices

rot_x_pi_over_2_band_2 = np.array([
    [0, 0, 0, -1, 0],
    [0, -1, 0, 0, 0],
    [0, 0, -1/2, 0, -sqrt(3)/2],
    [1, 0, 0, 0, 0],
    [0, 0, -sqrt(3)/2, 0, 1/2]
])

rotate_z_5d = lambda gamma: np.array([
    [cos(2 * gamma), 0, 0, 0, sin(2 * gamma)],
    [0, cos(gamma), 0, sin(gamma), 0],
    [0, 0, 1, 0, 0],
    [0, -sin(gamma), 0, cos(gamma), 0],
    [-sin(2 * gamma), 0, 0, 0, cos(2 * gamma)]
])

rot_x_pi_over_two_band_4 = np.array([
    [0, 0, 0, 0, 0, sqrt(14)/4, 0, - sqrt(2)/4, 0],
    [0, -3/4, 0, sqrt(7)/4, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, sqrt(2)/4, 0, sqrt(14)/4, 0],
    [0, sqrt(7)/4, 0, 3/4, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 3/8, 0, sqrt(5)/4, 0, sqrt(35)/8],
    [-sqrt(14)/4, 0, -sqrt(2)/4, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, sqrt(5)/4, 0, 1/2, 0, -sqrt(7)/4],
    [sqrt(2)/4, 0, -sqrt(14)/4, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, sqrt(35)/8, 0, -sqrt(7)/4, 0, 1/8]
])

rotate_z_9d = lambda gamma: np.array([
    [cos(4 * gamma), 0, 0, 0, 0, 0, 0, 0, sin(4 * gamma)],
    [0, cos(3 * gamma), 0, 0, 0, 0, 0 , sin(3 * gamma), 0],
    [0, 0, cos(2 * gamma), 0, 0, 0, sin(2 * gamma), 0, 0],
    [0, 0, 0, cos(gamma), 0, sin(gamma), 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, -sin(gamma), 0, cos(gamma), 0, 0, 0],
    [0, 0, -sin(2 * gamma), 0, 0, 0, cos(2 * gamma), 0, 0],
    [0, -sin(3 * gamma), 0, 0, 0, 0, 0, cos(3 * gamma), 0],
    [-sin(4 * gamma), 0, 0, 0, 0, 0, 0, 0, cos(4 * gamma)]
])

# lie angular momentum operators

L_x_5d = np.array([
    [0, 0, 0, -1, 0],
    [0, 0, -sqrt(3), 0, -1],
    [0, sqrt(3), 0, 0, 0],
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0]
])

L_y_5d = -np.array([
    [0, 1, 0, 0, 0],
    [-1, 0, 0, 0, 0],
    [0, 0, 0, -sqrt(3), 0],
    [0, 0, sqrt(3), 0, -1],
    [0, 0, 0, 1, 0]
])

L_z_5d = np.array([
    [0, 0, 0, 0, 2],
    [0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0],
    [0, -1, 0, 0, 0],
    [-2, 0, 0, 0, 0]
])

L_x_9d = np.array([
    [0, 0, 0, 0, 0, 0, 0, -sqrt(2), 0],
    [0, 0, 0, 0, 0, 0, -sqrt(7/2), 0, -sqrt(2)],
    [0, 0, 0, 0, 0, -3/sqrt(2), 0, -sqrt(7/2), 0],
    [0, 0, 0, 0, -sqrt(10), 0, -3/sqrt(2), 0, 0],
    [0, 0, 0, sqrt(10), 0, 0, 0, 0, 0],
    [0, 0, 3/sqrt(2), 0, 0, 0, 0, 0, 0],
    [0, sqrt(7/2), 0, 3/sqrt(2), 0, 0, 0, 0, 0],
    [sqrt(2), 0, sqrt(7/2), 0, 0, 0, 0, 0, 0],
    [0, sqrt(2), 0, 0, 0, 0, 0, 0, 0]
])

L_y_9d = -np.array([
    [0, sqrt(2), 0, 0, 0, 0, 0, 0, 0],
    [-sqrt(2), 0, sqrt(7/2), 0, 0, 0, 0, 0, 0],
    [0, -sqrt(7/2), 0, 3/sqrt(2), 0, 0, 0, 0, 0],
    [0, 0, -3/sqrt(2), 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, -sqrt(10), 0, 0, 0],
    [0, 0, 0, 0, sqrt(10), 0, -3/sqrt(2), 0, 0],
    [0, 0, 0, 0, 0, 3/sqrt(2), 0, -sqrt(7/2), 0],
    [0, 0, 0, 0, 0, 0, sqrt(7/2), 0, -sqrt(2)],
    [0, 0, 0, 0, 0, 0, 0, sqrt(2), 0]
])

L_z_9d = np.array([
    [0, 0, 0, 0, 0, 0, 0, 0, 4],
    [0, 0, 0, 0, 0, 0, 0, 3, 0],
    [0, 0, 0, 0, 0, 0, 2, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1, 0, 0, 0, 0, 0],
    [0, 0, -2, 0, 0, 0, 0, 0, 0],
    [0, -3, 0, 0, 0, 0, 0, 0, 0],
    [-4, 0, 0, 0, 0, 0, 0, 0, 0]
])

rotate_y_9d = lambda beta: rot_x_pi_over_two_band_4 @ rotate_z_9d(beta) @ rot_x_pi_over_two_band_4.T
rotate_x_9d = lambda alpha: rotate_y_9d(pi / 2).T @ rotate_z_9d(alpha) @ rotate_y_9d(pi / 2)

rot_x_pi_over_two_bands_combined = block_diag(1, rot_x_pi_over_2_band_2, rot_x_pi_over_two_band_4)

rotate_z_15d = lambda gamma: block_diag(1, rotate_z_5d(gamma), rotate_z_9d(gamma))
rotate_y_15d = lambda beta: rot_x_pi_over_two_bands_combined @ rotate_z_15d(beta) @ rot_x_pi_over_two_bands_combined.T
rotate_x_15d = lambda alpha: rotate_y_15d(pi / 2).T @ rotate_z_15d(alpha) @ rotate_y_15d(pi / 2)

def rotate_octahedral(a_0, alpha, beta, gamma):
    return rotate_x_9d(alpha) @ rotate_y_9d(beta) @ rotate_z_9d(gamma) @ a_0

lie_rotate_x_15d = lambda alpha: expm(alpha * block_diag(0, L_x_5d, L_x_9d))
lie_rotate_y_15d = lambda beta: expm(beta * block_diag(0, L_y_5d, L_y_9d))
lie_rotate_z_15d = lambda gamma: expm(gamma * block_diag(0, L_z_5d, L_z_9d))

def rotate_odeco(a_0, alpha, beta, gamma):
    R = rotate_x_15d(alpha) @ rotate_y_15d(beta) @ rotate_z_15d(gamma)
    Rl = lie_rotate_x_15d(alpha) @ lie_rotate_y_15d(beta) @ lie_rotate_z_15d(gamma)

    if not np.allclose(R, Rl, rtol=1e-7):
        raise Exception("Rotation matrices via lie algebra and wigner do not match up.")

    return lie_rotate_x_15d(alpha) @ lie_rotate_y_15d(beta) @ lie_rotate_z_15d(gamma) @ a_0
