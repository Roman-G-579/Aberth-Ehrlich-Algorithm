import subprocess

def call_scheme_function(func, args, output_file="output.txt"):
    """Call a Scheme function with arguments and write the result to a file."""
    result = subprocess.run(['racket', 'scheme_functions.rkt', func] + args + [output_file], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Scheme error: {result.stderr}")
    
    with open(output_file, 'r') as f:
        result = f.read().strip()
    
    return parse_result(result)

def parse_result(result):
    """Parse the result from Scheme."""
    if " " in result:
        elements = result.split(" ")
        if "," in elements[0]:
            return [complex(float(r), float(i)) for r, i in (elem.split(",") for elem in elements)]
        else:
            return [float(elem) for elem in elements]
    else:
        return float(result)

def polynomial_function(coefficients):
    """Return a callable function to evaluate the polynomial with the given coefficients."""
    def poly(x):
        return call_scheme_function('polynomial-eval', [",".join(map(str, coefficients)), str(x)])
    return poly

def polynomial_derivative_function(coefficients):
    """Return a callable function to evaluate the derivative of the polynomial with the given coefficients."""
    derivative_coeffs = call_scheme_function('polynomial-derivative', [",".join(map(str, coefficients))])
    return polynomial_function(derivative_coeffs)

def random_approximations(coefficients):
    """Generate initial approximations using Scheme."""
    return call_scheme_function('random-approximations', [",".join(map(str, coefficients))])

def negligible_complex(expression, epsilon):
    """Check if a complex number is negligible (close to zero)."""
    return abs(expression.real) < epsilon and abs(expression.imag) < epsilon

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
                sigma = sum(1 / (zk - zj) for j, zj in enumerate(random_guesses) if k != j)
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
roots = aberth_method(f_0, f_1, coefficients)
print("Roots:", roots)
