import numpy as np
from scipy.linalg import block_diag
from math import *


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

rotate_y_9d = lambda beta: rot_x_pi_over_two_band_4 @ rotate_z_9d(beta) @ rot_x_pi_over_two_band_4.T
rotate_x_9d = lambda alpha: rotate_y_9d(pi / 2).T @ rotate_z_9d(alpha) @ rotate_y_9d(pi / 2)

rot_x_pi_over_two_bands_combined = block_diag(1, rot_x_pi_over_2_band_2, rot_x_pi_over_two_band_4)

rotate_z_15d = lambda gamma: block_diag(1, rotate_z_5d(gamma), rotate_z_9d(gamma))
rotate_y_15d = lambda beta: rot_x_pi_over_two_bands_combined @ rotate_z_15d(beta) @ rot_x_pi_over_two_bands_combined.T
rotate_x_15d = lambda alpha: rotate_y_15d(pi / 2).T @ rotate_z_15d(alpha) @ rotate_y_15d(pi / 2)

def rotate_octahedral(a_0, alpha, beta, gamma):
    return rotate_x_9d(alpha) @ rotate_y_9d(beta) @ rotate_z_9d(gamma) @ a_0


def rotate_odeco(a_0, alpha, beta, gamma):
    return rotate_x_15d(alpha) @ rotate_y_15d(beta) @ rotate_z_15d(gamma) @ a_0