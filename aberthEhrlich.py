import time
from mpmath import mp, mpc

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
     return 4*x**2 - 10*x + 4
    # 2x^4 - 3x^3 + 2x^2 - 7x + 8
    #return 2 * x ** 4 - 3 * x ** 3 + 2 * x ** 2 - 7 * x + 8


# derivative of function
def df(x):
    return 8 * x - 10
    #return 8 * x ** 3 - 9 * x ** 2 + 4 * x - 7


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


def get_approximations(coefficients, epsilon=0.000001):
    """
    Calculates the initial root approximations for the given polynomial coefficients
    :param coefficients: the coefficients array
    :return: the array of approximations in complex number form
    """
    # Polynomial's degree
    deg = len(coefficients) - 1

    # Constant coefficient (lower bound)
    p_0 = coefficients[0]

    # Leading coefficient (upper bound)
    p_n = coefficients[deg]

    # Complex plane's radius
    r = abs((p_0 / p_n) ** (1 / deg))

    # The approximations array
    approximations = []

    for i in range(deg):
        # The angle that is sent to the sine and cosine functions
        angle = normalize_degree((2 * pi() * i))

        cos_trig = cos(angle)
        # cos_trig = 0 if is_close_to_zero(cos_trig, epsilon) else cos_trig

        sin_trig = sin(angle)
        # sin_trig = 0 if is_close_to_zero(sin_trig, epsilon) else cos_trig

        # the approximation of the current coefficient
        approximations.append(complex(r * cos_trig, r * sin_trig))

    return approximations


def compute_new_roots(approximations):
    # updated roots
    roots = []
    # number of roots
    n = len(approximations)
    while True:
        valid = 0
        for k in range(n):
            # root approximation in index k
            z_k = approximations[k]

            # results of calling function and its derivative with the approximation as parameter
            # f_i, df_i = f_x_and_df_x(z_k, coefficients)

            f_i = f(z_k)
            df_i = df(z_k)

            ratio = f_i / df_i

            offset = ratio / (1 - (ratio * sum(1/(z_k - z_j)
                                               for j, z_j in enumerate(approximations) if j != k)))

            # if round(offset.real, 14) == 0 and round(offset.imag, 14) == 0:
            if is_close_to_zero(offset, 0.000001):
                valid += 1
            approximations[k] -= offset
        if valid == len(approximations):
            break

    return approximations


def aberth_method(coefficients, epsilon=0.000001):
    approximations = get_approximations(coefficients)

    new_roots = []

    for i in range(1000):
        new_roots = compute_new_roots(approximations)

    return new_roots


# the functions coefficients sorted by ascending degree size
# coefficients = [8, -7, 2, -3, 2]
coefficients = [4, -10, 4]
# coefficients = [3, 2, -1]
# coefficients = read_coefficients('poly_coeff(997).txt')
print(aberth_method(coefficients))

# start = time.time()
# print(cos(30))
# end = time.time()
# print(end - start)


