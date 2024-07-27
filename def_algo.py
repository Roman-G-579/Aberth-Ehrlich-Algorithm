# Implementation of Aberth's method in Python, written by Jonathan.

from numpy import linspace
from math import sin,cos,pi
from typing import Callable

def random_approximations(coefficients):
    """ הפעולה מחזירה רשימה של קירובים התחלתיים, בהתאם לרשימת המקדמים אשר היא מקבלת"""
    n = len(coefficients) - 1
    radius = abs(coefficients[-1] / coefficients[0]) ** (1 / n)
    return [complex(radius * cos(angle), radius * sin(angle)) for angle in linspace(0, 2 * pi, n)]

def negligible_complex(expression: complex, epsilon:float)->bool:
    """ הפעולה מחזירה האם מספר מרוכב כה קטן כך שהוא זניח וניתן להחשיבו כ-0"""
    return abs(expression.real) < epsilon and abs(expression.imag) < epsilon

def aberth_method(f_0: Callable, f_1: Callable, coefficients, epsilon=0.000001, nmax=100000):
    """ הפעולה מקבלת את הפונקציה, את נגזרתה, ואת המקדמים של הפונקציה. בנוסף מקבלת אפסילון אשר קובע את ההפרש שלפיו מחשיבים מספר כזניח ואת מספר האיטרציות המקסימלי"""
    try:
        random_guesses = random_approximations(coefficients)
        for n in range(nmax):
            offsets = []
            for k, zk in enumerate(random_guesses):
                m = f_0(zk) / f_1(zk)  # save it as m, so it won't be calculated many times
                sigma = sum(1 / (zk - zj) for j, zj in enumerate(random_guesses) if k != j and zk != zj)
                denominator = 1 - m * sigma
                offsets.append(m / denominator)
            random_guesses = [approximation - offset for approximation, offset in zip(random_guesses, offsets)]
            if all(negligible_complex(f_0(guess), epsilon) for guess in random_guesses):
                break
        return set(random_guesses)
    except ValueError:
        return set()
    
def polynomial_function(coefficients):
    """Return a callable function to evaluate the polynomial with the given coefficients."""
    def poly(x):
        result = coefficients[0]
        for coefficient in coefficients[1:]:
            result = result * x + coefficient
        return result
    return poly

def polynomial_derivative_function(coefficients):
    """Return a callable function to evaluate the derivative of the polynomial with the given coefficients."""
    derivative_coeffs = [
        coefficients[i] * (len(coefficients) - i - 1)
        for i in range(len(coefficients) - 1)
    ]
    return polynomial_function(derivative_coeffs)

def read_coefficients_from_file(filename):
    """Read polynomial coefficients from a text file."""
    with open(filename, 'r') as file:
        coefficients = [float(line.strip()) for line in file.readlines()]
    return coefficients

coefficients = read_coefficients_from_file('test.txt')
f_0 = polynomial_function(coefficients)
f_1 = polynomial_derivative_function(coefficients)
roots = aberth_method(f_0, f_1, coefficients)
print(roots)
    
