import time
import sys

sys.setrecursionlimit(2000)


def round_complex(complex_num):
    return [round(num.real, 3) + round(num.imag, 3) * 1j for num in complex_num]


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


def cos(x, precision=50):
    # ( (-1)^n ) * x^2n / (2n)!
    term = 1
    cos_x = term

    # Precompute x^2
    x_squared = x * x

    # Alternating sign for the series
    sign = -1

    for n in range(1, precision):
        term *= x_squared / ((2 * n - 1) * (2 * n))
        cos_x += sign * term
        sign = -sign
    return cos_x


def sin(x, precision=50):
    # ( (-1)^n ) * x^(2n+1) / (2n+1)!
    term = x
    sin_x = term

    # Precompute x^2
    x_squared = x * x

    # Alternating sign for the series
    sign = -1

    for n in range(1, precision):
        term *= x_squared / ((2 * n) * (2 * n + 1))
        sin_x += sign * term
        sign = -sign
    return sin_x


def pi():
    return 3.141592653589793


def poly_val_rec(p, x, i=0):  # n = deg(p)
    if i == len(p):
        return 0
    else:
        return p[i] + x * poly_val_rec(p, x, i + 1)


def poly_val(p, x):
    s = 0
    for i in range(len(p)):
        s = s * x + p[i]
    return s


def frac_val(p, q, x):
    """
    Calculates the value of p(z_k)/p'(z_k) for the Aberth-Ehrlich algorithm
    :param p: The polynomial's coefficients
    :param q: The polynomial derivative's coefficients
    :param x: The value of z_k to insert as the function parameter
    :return: The result fraction
    """
    if abs(x) <= 1:
        return poly_val(p, x) / poly_val(q, x)
    else:
        # Handles x values larger than 1, where attempting to solve the function using the regular method
        # would result in a number overflow for large polynomials
        return x * (poly_val(p[::-1], 1 / x)) / poly_val(q[::-1], 1 / x)


def calculate_sigma(approximations, k):
    zk = approximations[k]
    sigma = 0.0
    for j, zj in enumerate(approximations):
        if j != k and zk != zj:
            sigma += 1 / (zk - zj)
    return sigma


def polynomial_derivative_coefficients(coefficients):
    """
    Calculates the coefficients of a polynomial's derivative based on the given coefficients
    :param coefficients: the polynomial's coefficients array
    :return: the derivative's coefficients array
    """
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

    # Precompute pi
    pi_val = pi()

    # Constant coefficient (lower bound)
    p_0 = coefficients[deg]

    # Leading coefficient (upper bound)
    p_n = coefficients[0]

    # Complex plane's radius
    r = abs((p_0 / p_n) ** (1 / deg))

    # Distributes the numbers across the complex plane
    angles = [x * (2 * pi_val) / (deg-1) for x in range(deg)]
    return [complex(r * cos(angle), r * sin(angle)) for angle in angles]


def aberth_method(coefficients, epsilon=0.0001):
    """
    Finds the roots of the given polynomial using the Aberth-Ehrlich algorithm
    :param coefficients: An array of the polynomial's coefficients in descending order (First element has
                         the highest degree, last element is constant)
    :param epsilon: The precision level of the algorithm
    :return: Array of roots
    """
    # Initial root approximations
    approximations = get_approximations(coefficients)

    # The process is repeated a number of times for better results
    for n in range(500):
        # Array of offsets (w_k)
        offsets = []

        # Finds the offset for every approximation
        for k, zk in enumerate(approximations):

            frac = frac_val(coefficients, derivative_coefficients, zk)

            sigma = calculate_sigma(approximations, k)

            denominator = 1 - frac * sigma

            offsets.append(frac / denominator)

        for i in range(len(approximations)):
            approximations[i] -= offsets[i]

        # Checks if all values are close to zero
        roots_converge = True
        for val in approximations:
            if not is_close_to_zero(poly_val(coefficients, val), epsilon):
                roots_converge = False
                break

        # Once all roots have converged, the process is done
        if roots_converge:
            print('roots converged at iteration ', n)
            break

    return approximations


# coefficients = [4, 6, 8, -10, 4]

coefficients = read_coefficients('poly_coeff(997).txt')  # Polynomial's coefficients
derivative_coefficients = polynomial_derivative_coefficients(coefficients)  # Polynomial derivative's coefficients

start = time.time()
result = (aberth_method(coefficients))
end = time.time()
print(end - start)

print(round_complex(result))
