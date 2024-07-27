from time import time

def linspace(start, stop, num):
    """Generate `num` evenly spaced values from `start` to `stop`."""
    if num == 1:
        return [start]
    step = (stop - start) / (num - 1)
    return [start + step * i for i in range(num)]

def cos(x):
    """Approximate the cosine function using the Taylor series expansion."""
    x = x % (2 * pi())
    result = 1
    term = 1
    x_squared = x * x
    for i in range(1, 50):
        term *= -x_squared / ((2 * i) * (2 * i - 1))
        result += term
    return result

def sin(x):
    """Approximate the sine function using the Taylor series expansion."""
    x = x % (2 * pi())
    result = x
    term = x
    x_squared = x * x
    for i in range(1, 50):
        term *= -x_squared / ((2 * i + 1) * (2 * i))
        result += term
    return result

def pi():
    """Calculate pi using the Leibniz formula for pi."""
    pi_approx = 0
    for k in range(100000):
        pi_approx += (-1) ** k / (2 * k + 1)
    return 4 * pi_approx

def random_approximations(coefficients):
    """Return a list of initial approximations based on the given coefficients."""
    n = len(coefficients) - 1
    radius = abs(coefficients[-1] / coefficients[0]) ** (1 / n)
    pi_val = pi()
    return [complex(radius * cos(angle), radius * sin(angle)) for angle in linspace(0, 2 * pi_val, n)]

def negligible_complex(expression, epsilon):
    """Check if a complex number is negligible (close to zero)."""
    return abs(expression.real) < epsilon and abs(expression.imag) < epsilon

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

def aberth_method(f_0, f_1, coefficients, epsilon=0.000001, nmax=100000):
    """Find the roots of a polynomial using Aberth's method."""
    try:
        random_guesses = random_approximations(coefficients)
        for n in range(nmax):
            offsets = []
            for k, zk in enumerate(random_guesses):
                f_val = f_0(zk)
                f_prime_val = f_1(zk)
                m = f_val / f_prime_val  # save it as m, so it won't be calculated many times
                sigma = sum(1 / (zk - zj) for j, zj in enumerate(random_guesses) if k != j and zk != zj)
                denominator = 1 - m * sigma
                offsets.append(m / denominator)
            random_guesses = [approximation - offset for approximation, offset in zip(random_guesses, offsets)]
            if all(negligible_complex(f_0(guess), epsilon) for guess in random_guesses):
                break
        return set(random_guesses)
    except ValueError:
        return set()

def read_coefficients_from_file(filename):
    """Read polynomial coefficients from a text file."""
    with open(filename, 'r') as file:
        coefficients = [float(line.strip()) for line in file.readlines()]
    return coefficients

# Example usage
filename = 'test.txt'
coefficients = read_coefficients_from_file(filename)
f_0 = polynomial_function(coefficients)
f_1 = polynomial_derivative_function(coefficients)
start_time = time()
roots = aberth_method(f_0, f_1, coefficients)
print("Time:", time() - start_time)
print("Roots:", roots)