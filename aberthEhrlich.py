import time
import sys

from numpy import linspace
import math
sys.setrecursionlimit(10 ** 6)


def round_complex(complex_num):
    return [round(num.real, 2) + round(num.imag, 2) * 1j for num in complex_num]


def read_coefficients(filename):
    """
    Reads a coefficients array from a text file
    :param filename: the text file
    :return: an array of all the coefficients in the file
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
        coefficients = [float(line.strip()) for line in lines]
    return coefficients


def is_close_to_zero(value, epsilon):
    # print(abs(value.real))
    # print(abs(value.imag))
    # print(abs(value.real) < epsilon and abs(value.imag) < epsilon)
    return abs(value.real) < epsilon and abs(value.imag) < epsilon


def normalize_degree(angle):
    """
    Normalize an angle to be within the range [0, 360) degrees.

    Parameters:
    angle (float): The angle in degrees.

    Returns:
    float: The normalized angle.
    """
    while angle >= 360:
        angle -= 360
    while angle < 0:
        angle += 360
    return angle


def cos(x, precision=50):
    term = 1
    cos_x = term
    for n in range(1, precision):
        term *= -x ** 2 / ((2 * n - 1) * (2 * n))
        cos_x += term
    return cos_x


def sin(x, precision=50):
    term = x
    sin_x = term
    for n in range(1, precision):
        term *= -x ** 2 / ((2 * n) * (2 * n + 1))
        sin_x += term
    return sin_x


def pi():
    return 3.141592653589793


# function to solve
def f(x):
    return 4 * x ** 2 - 10 * x + 4


# derivative of function
def df(x):
    return 8 * x - 10
    #return 8 * x ** 3 - 9 * x ** 2 + 4 * x - 7


# def poly_val(coefficients, x, n=None):
#     if n is None:
#         n = len(coefficients)
#
#     if n == 1:
#         return coefficients[0]
#
#     return coefficients[0] + x * poly_val(coefficients[1:], x, n - 1)


# def poly_val(p, x, n):  # n = deg(p)
#     if n == 0:
#         return p[0]
#     else:
#         return p[n] + poly_val(p, x, n - 1) * x

def poly_val(p, x):
    s = 0
    for i in range(len(p)):
        s = s * x + p[i]
    return s


def frac_val(p, q, x):
    #deg_p = len(p) - 1
    #deg_q = len(q) - 1
    if abs(x) <= 1:
        return poly_val(p, x) / poly_val(q, x)
    else:
        return x * (poly_val(p[::-1], 1 / x) / poly_val(q[::-1], 1 / x))


def f_x_and_df_x(x, coefficients):
    """
    Evaluate a polynomial and its derivative at a given value of x,
    using Horner's method

    :param x: The value at which to evaluate f() and df()
    :param coefficients: The coefficients of the polynomial,
    ordered from the highest degree term to the constant term.

    :return: f(x) and df(x)
    """
    n = len(coefficients)
    if n == 0:
        return 0, 0

    result = complex(coefficients[0])
    derivative = 0

    for i in range(1, n):
        derivative = derivative * x + result
        result = result * x + coefficients[i]

    return result, derivative


def polynomial_derivative_coefficients(coefficients):
    # The length of the coefficients array
    n = len(coefficients)

    # If the polynomial is a constant (degree 0), its derivative is 0
    if n == 1:
        return [0]

    # Compute the coefficients of the derivative
    derivative_coeffs = [(n - 1 - i) * coefficients[i] for i in range(n - 1)]

    return derivative_coeffs


def get_approximations(coefficients):
    """
    Calculates the initial root approximations for the given polynomial coefficients
    :param coefficients: the coefficients array
    :return: the array of approximations in complex number form
    """
    # Polynomial's degree
    deg = len(coefficients) - 1

    step = 2 * pi() / deg

    # Constant coefficient (lower bound)
    p_0 = coefficients[deg]

    # Leading coefficient (upper bound)
    p_n = coefficients[0]

    # Complex plane's radius
    r = abs((p_0 / p_n) ** (1 / deg))

    approximations = []
    angles = [0 + x * (2 * pi()) / (deg-1) for x in range(deg)]
    return [complex(r * cos(angle), r * sin(angle)) for angle in angles]
    # # The approximations array
    #
    # for i in range(deg):
    #     # The angle that is sent to the sine and cosine functions
    #     angle = i * step
    #
    #     cos_trig = cos(angle)
    #     sin_trig = sin(angle)
    #
    #     # the approximation of the current coefficient
    #     [0 + angle * (2 * pi() ) / (length - 1) for x in range(deg)]
    #
    #     approximations.append(complex(r * cos_trig, r * sin_trig))
    #
    # return approximations


def aberth_method(coefficients, epsilon=0.0001):
    approximations = get_approximations(coefficients)

    # results array
    roots = []
    n = len(approximations)
    valid_roots = [False] * n
    # valid_root_cnt = 0
    for i in range(50):
        while True:
            valid = 0
            for k in range(n):
                # if valid_roots[k]:
                #     continue
                # root approximation in index k
                z_k = approximations[k]

                # f(z_k) / df(z_k)
                frac = frac_val(coefficients, derivative_coefficients, z_k)

                offset = frac / (1 - (frac * sum(1 / (z_k - z_j) for j, z_j in enumerate(approximations) if j != k)))

                # if round(offset.real, 14) == 0 and round(offset.imag, 14) == 0:
                if is_close_to_zero(offset, epsilon):
                    # valid_roots[k] = True
                    valid += 1
                # print(approximations[k] - offset)
                approximations[k] -= offset
            # print('valid cnt: ', sum(valid_roots))
            # if all(valid_roots):
            if valid == len(approximations):
                break
        roots = approximations
    return roots


# coefficients = [4, -10, 4]

coefficients = read_coefficients('poly_coeff(997).txt')  # Polynomial's coefficients
derivative_coefficients = polynomial_derivative_coefficients(coefficients)  # Polynomial derivative's coefficients

start = time.time()
result = (aberth_method(coefficients))
end = time.time()
print(end - start)

print(result)
print(round_complex(result))
