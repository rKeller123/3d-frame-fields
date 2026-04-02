
from math import cos, sin, sqrt, pi
import numpy as np

rotate_sh_z = lambda gamma: np.array([
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

R_x_B_pi_over_two = np.array([
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

rotate_sh_y = lambda beta: R_x_B_pi_over_two @ rotate_sh_z(beta) @ R_x_B_pi_over_two.T
rotate_sh_x = lambda alpha: rotate_sh_y(pi / 2).T @ rotate_sh_z(alpha) @ rotate_sh_y(pi / 2)

rotate_sh = lambda alpha, beta, gamma : rotate_sh_x(alpha) @ rotate_sh_y(beta) @ rotate_sh_z(gamma)