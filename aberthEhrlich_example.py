# Implementation of Aberth's method in Python, written by Jonathan.
import time

from numpy import linspace
from math import sin, cos, pi


def round_complex(complex_num):
    return [round(num.real, 2) + round(num.imag, 2) * 1j for num in complex_num]


def read_coefficients(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        coefficients = [float(line.strip()) for line in lines]
    return coefficients


def polynomial_derivative_coefficients(coefficients):
    # The length of the coefficients array
    n = len(coefficients)

    # If the polynomial is a constant (degree 0), its derivative is 0
    if n == 1:
        return [0]

    # Compute the coefficients of the derivative
    derivative_coeffs = [(n - 1 - i) * coefficients[i] for i in range(n - 1)]

    return derivative_coeffs


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


# הפעולה מחזירה רשימה של קירובים התחלתיים, בהתאם לרשימת המקדמים אשר היא מקבלת
def random_approximations(coefficients):
    n = len(coefficients) - 1
    r = abs(coefficients[-1] / coefficients[0]) ** (1 / n)

    approximations = []

    for angle in linspace(0, 2*pi, n):
        cos_trig = cos(angle)
        sin_trig = sin(angle)
        approximations.append(complex(r * cos_trig, r * sin_trig))
    return approximations
    # return [complex(radius * cos(angle), radius * sin(angle)) for angle in linspace(0, 2 * pi, n)]


# הפעולה מחזירה האם מספר מרוכב כה קטן כך שהוא זניח וניתן להחשיבו כ-0
def negligible_complex(expression: complex, epsilon: float) -> bool:
    return abs(expression.real) < epsilon and abs(expression.imag) < epsilon


# הפעולה מקבלת את הפונקציה, את נגזרתה, ואת המקדמים של הפונקציה. בנוסף מקבלת אפסילון אשר קובע את ההפרש שלפיו מחשיבים
# מספר כזניח ואת מספר האיטרציות המקסימלי
def aberth_method(f_0, f_1, coefficients, epsilon=0.000001, nmax=100):
    try:
        random_guesses = round_complex(random_approximations(coefficients))
        #print(random_guesses)
        for n in range(nmax):
            offsets = []
            for k, zk in enumerate(random_guesses):

                # f, df = f_x_and_df_x(zk, coefficients)
                # m = f_0(zk) / f_1(zk)  # save it as m, so it won't be calculated many times
                # m = f / df
                m = frac_val(coefficients, derivative_coefficients, zk)

                sigma = sum(1 / (zk - zj) for j, zj in enumerate(random_guesses) if k != j and zk != zj)
                denominator = 1 - m * sigma
                offsets.append(m / denominator)
            random_guesses = [approximation - offset for approximation, offset in zip(random_guesses, offsets)]
            if all(negligible_complex(poly_val(coefficients, guess), epsilon) for guess in random_guesses):
                break
        return random_guesses
    except ValueError:
        return set()


# function to solve
# def f_0(x):
#     return 4*x**2 - 10*x + 4
#     # return 2 * x ** 2 - 7 * x + 8
#     #return 2 * x ** 4 - 3 * x ** 3 + 2 * x ** 2 - 7 * x + 8

def f_0(x, coefficients):
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

    for i in range(1, n):
        result = result * x + coefficients[i]

    return result


# derivative of function
def f_1(x):
    return 8*x - 10
    # return 4 * x - 7
    # return 8 * x ** 3 - 9 * x ** 2 + 4 * x - 7


# def function_from_file(x):
#

# the functions coefficients sorted by ascending degree size
# coefficients = [8, -7, 2, -3, 2]

coefficients = [4, 6, 8, -10, 4]
# coefficients = read_coefficients('poly_coeff(997).txt')
derivative_coefficients = polynomial_derivative_coefficients(coefficients)  # Polynomial derivative's coefficients

# print(random_approximations(coefficients))

# coefficients_from_file = read_coefficients('poly_coeff(997).txt')
# print(coefficients_from_file)

start = time.time()
# print(random_approximations(coefficients))
# result = aberth_method(f_0, f_1, coefficients)
result = (aberth_method(f_0, f_1, coefficients))
print(round_complex(result))
end = time.time()




# print(round_complex(result))
print(end - start)
