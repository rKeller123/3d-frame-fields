import sympy as sp
from sympy import sqrt, Integer
from scipy.linalg import expm
import numpy as np

c, s, c2, c3, c4, s2, s3, s4 = sp.symbols("c,s,c2,c3,c4,s2,s3,s4")

L_x_5d = sp.Matrix([
    [0, 0, 0, -1, 0],
    [0, 0, -sqrt(3), 0, -1],
    [0, sqrt(3), 0, 0, 0],
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0]
])

L_y_5d = sp.Matrix([
    [0, 1, 0, 0, 0],
    [-1, 0, 0, 0, 0],
    [0, 0, 0, -sqrt(3), 0],
    [0, 0, sqrt(3), 0, -1],
    [0, 0, 0, 1, 0]
])

L_z_5d = sp.Matrix([
    [0, 0, 0, 0, 2],
    [0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0],
    [0, -1, 0, 0, 0],
    [-2, 0, 0, 0, 0]
])

L_x_9d = sp.Matrix([
    [0, 0, 0, 0, 0, 0, 0, -sqrt(2), 0],
    [0, 0, 0, 0, 0, 0, -sqrt(Integer(7)/2), 0, -sqrt(2)],
    [0, 0, 0, 0, 0, -Integer(3)/sqrt(2), 0, -sqrt(Integer(7)/2), 0],
    [0, 0, 0, 0, -sqrt(10), 0, -Integer(3)/sqrt(2), 0, 0],
    [0, 0, 0, sqrt(10), 0, 0, 0, 0, 0],
    [0, 0, 3/sqrt(2), 0, 0, 0, 0, 0, 0],
    [0, sqrt(Integer(7)/2), 0, Integer(3)/sqrt(2), 0, 0, 0, 0, 0],
    [sqrt(2), 0, sqrt(Integer(7)/2), 0, 0, 0, 0, 0, 0],
    [0, sqrt(2), 0, 0, 0, 0, 0, 0, 0]
])

L_y_9d = sp.Matrix([
    [0, sqrt(2), 0, 0, 0, 0, 0, 0, 0],
    [-sqrt(2), 0, sqrt(Integer(7)/2), 0, 0, 0, 0, 0, 0],
    [0, -sqrt(Integer(7)/2), 0, Integer(3)/sqrt(2), 0, 0, 0, 0, 0],
    [0, 0, -Integer(3)/sqrt(2), 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, -sqrt(10), 0, 0, 0],
    [0, 0, 0, 0, sqrt(10), 0, -3/sqrt(2), 0, 0],
    [0, 0, 0, 0, 0, 3/sqrt(2), 0, -sqrt(Integer(7)/2), 0],
    [0, 0, 0, 0, 0, 0, sqrt(Integer(7)/2), 0, -sqrt(2)],
    [0, 0, 0, 0, 0, 0, 0, sqrt(2), 0]
])

L_z_9d = sp.Matrix([
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

L_x = sp.diag(
    Integer(0),
    L_x_5d,
    L_x_9d
)

L_y = sp.diag(
    Integer(0),
    L_y_5d,
    L_y_9d
)

L_z = sp.diag(
    Integer(0),
    L_z_5d,
    L_z_9d
)

rotate_z_5d = sp.Matrix([
    [c2, 0, 0, 0, s2],
    [0, c, 0, s, 0],
    [0, 0, 1, 0, 0],
    [0, -s, 0, c, 0],
    [-s2, 0, 0, 0, c2]
])

rotate_z_9d = sp.Matrix([
    [c4, 0, 0, 0, 0, 0, 0, 0, s4],
    [0, c3, 0, 0, 0, 0, 0 , s3, 0],
    [0, 0, c2, 0, 0, 0, s2, 0, 0],
    [0, 0, 0, c, 0, s, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, -s, 0, c, 0, 0, 0],
    [0, 0, -s2, 0, 0, 0, c2, 0, 0],
    [0, -s3, 0, 0, 0, 0, 0, c3, 0],
    [-s4, 0, 0, 0, 0, 0, 0, 0, c4]
])

rotate_z = sp.diag(
    sp.Integer(1),
    rotate_z_5d,
    rotate_z_9d
)

rot_x_pi_over_two_band_2 = sp.Matrix([
    [0, 0, 0, -1, 0],
    [0, -1, 0, 0, 0],
    [0, 0, -Integer(1)/2, 0, -sqrt(3)/2],
    [1, 0, 0, 0, 0],
    [0, 0, -sqrt(3)/2, 0, Integer(1)/2]
])

rot_x_pi_over_two_band_4 = sp.Matrix([
    [0, 0, 0, 0, 0, sqrt(14)/4, 0, - sqrt(2)/4, 0],
    [0, -Integer(3)/4, 0, sqrt(7)/4, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, sqrt(2)/4, 0, sqrt(14)/4, 0],
    [0, sqrt(7)/4, 0, Integer(3)/4, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, Integer(3)/8, 0, sqrt(5)/4, 0, sqrt(35)/8],
    [-sqrt(14)/4, 0, -sqrt(2)/4, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, sqrt(5)/4, 0, 1/2, 0, -sqrt(7)/4],
    [sqrt(2)/4, 0, -sqrt(14)/4, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, sqrt(35)/8, 0, -sqrt(7)/4, 0, Integer(1)/8]
])

rot_x_pi_2 = sp.diag(
    Integer(1),
    rot_x_pi_over_two_band_2,
    rot_x_pi_over_two_band_4
)


rotate_y = rot_x_pi_2 * rotate_z.T * rot_x_pi_2.T

subs = {
    s: 1,
    c: 0,
    s2: 0,
    c2: -1,
    s3: -1,
    c3: 0,
    s4: 0,
    c4: 1,
}

rot_y_pi_over_two = rotate_y.subs(subs)

rotate_x = rot_y_pi_over_two.T * rotate_z.T * rot_y_pi_over_two


def get_trig_subs(theta):
    return {
        c: np.cos(theta),  s: np.sin(theta),
        c2: np.cos(2*theta), s2: np.sin(2*theta),
        c3: np.cos(3*theta), s3: np.sin(3*theta),
        c4: np.cos(4*theta), s4: np.sin(4*theta),
    }

get_rot_z = lambda theta: rotate_z.subs(get_trig_subs(theta)).evalf()
get_rot_y = lambda theta: rotate_y.subs(get_trig_subs(theta)).evalf()
get_rot_x = lambda theta: rotate_x.subs(get_trig_subs(theta)).evalf()

def test_lie_vs_wig(angle, lie_matrix, wig_matrix_func):
    wig_eval = wig_matrix_func(angle)
    wig_num = np.array(wig_eval, dtype=np.float64)
    
    lie_num = np.array(lie_matrix.evalf(), dtype=np.float64)
    rot_from_lie = expm(angle * lie_num)
    
    # 3. Calculate Frobenius norm of difference
    diff_norm = np.linalg.norm(rot_from_lie - wig_num, ord='fro')
    
    return diff_norm

# --- Example Usage ---

n = 24
for i in range(n):
    theta = i * np.pi / n

    err_z = test_lie_vs_wig(theta, L_z, get_rot_z)
    err_y = test_lie_vs_wig(theta, L_y, get_rot_y)
    err_x = test_lie_vs_wig(theta, L_x, get_rot_x)

    print(f"Error for Z-rotation: {err_z:.2e}")
    print(f"Error for Y-rotation: {err_y:.2e}")
    print(f"Error for X-rotation: {err_x:.2e}")