import math
import cmath

def legendre(l: int, m: int, x: float):
    assert l >= 0, f"l is negative: {l}"
    assert -l <= m <= l, f"m is not in the range required: -{l} <= {m} <= {l}"
    m_abs = abs(m)
    sum_val = 0
    for k in range((l - m_abs) // 2 + 1):
        sum_val += (
        (math.pow(-1, k) / (math.factorial(k) * math.factorial(l - k))) *
        (math.factorial(2 * (l - k)) / math.factorial(l - 2*k - m)) *
        math.pow(x, l - 2 * k - m)
        )
    leg_abs = (math.pow(-1, m_abs) / (math.pow(2, l))) * math.pow(1 - x * x, m_abs / 2) * sum_val

    if m < 0:
        return math.pow(-1, m) * (math.factorial(l - m_abs) / math.factorial(l + m_abs)) * leg_abs
    return leg_abs

def sph_harm(l: int, m: int, theta: float, phi: float):
    assert l >= 0, f"l is negative: {l}"
    assert -l <= m <= l, f"m is not in the range required: -{l} <= {m} <= {l}"
    return math.sqrt(
        (2 * l + 1) / (4 * math.pi) *
        math.factorial(l - m) / math.factorial(l + m)
    ) * legendre(l, m, math.cos(theta)) * cmath.exp(1j * m * phi)

def sph_harm_real(l: int, m: int, theta: float, phi: float):
    assert l >= 0, f"l is negative: {l}"
    assert -l <= m <= l, f"m is not in the range required: -{l} <= {m} <= {l}"
    if m < 0:
        return math.sqrt(2) * math.pow(-1, m) * sph_harm(l, abs(m), theta, phi).imag
    elif m == 0:
        return sph_harm(l, 0, theta, phi).real
    else:
        return math.sqrt(2) * math.pow(-1, m) * sph_harm(l, m, theta, phi).real
