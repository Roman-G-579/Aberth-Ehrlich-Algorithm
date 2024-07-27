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
        random_guesses = random_approximations(coefficients)
        #print(random_guesses)
        for n in range(nmax):
            offsets = []
            for k, zk in enumerate(random_guesses):

                f, df = f_x_and_df_x(zk, coefficients)
                # m = f_0(zk) / f_1(zk)  # save it as m, so it won't be calculated many times
                m = f / df
                sigma = sum(1 / (zk - zj) for j, zj in enumerate(random_guesses) if k != j and zk != zj)
                denominator = 1 - m * sigma
                offsets.append(m / denominator)
            random_guesses = [approximation - offset for approximation, offset in zip(random_guesses, offsets)]
            if all(negligible_complex(f_0(guess, coefficients), epsilon) for guess in random_guesses):
                break
        return set(random_guesses)
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
# coefficients = [4, -10, 4]
coefficients = read_coefficients('poly_coeff(997).txt')
# print(random_approximations(coefficients))

# coefficients_from_file = read_coefficients('poly_coeff(997).txt')
# print(coefficients_from_file)

start = time.time()
# print(random_approximations(coefficients))
# result = aberth_method(f_0, f_1, coefficients)
print(aberth_method(f_0, f_1, coefficients))
end = time.time()




# print(round_complex(result))
print(end - start)
