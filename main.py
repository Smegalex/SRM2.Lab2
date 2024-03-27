from sympy import diff, pprint, Interval, Abs, solve
from sympy.plotting import plot
from sympy.calculus.util import maximum
from sympy.abc import x, y
from math import factorial


def basic_checks(X: list, Y: list = None, funct=None):
    if Y == None and funct == None:
        raise (ValueError, "Недостатньо аргументів")
    elif Y == None and funct != None:
        Y = []
        for i in X:
            Y.append(funct.subs(x, i))

    if len(Y) != len(X):
        raise (
            ValueError, "Кількість аргументів X не відповідає кількості наданих значень Y")

    return X, Y


def multiply_polynomials(polynomial1: dict, polynomial2: dict):
    resulting_polynomial = {}
    for key1, value1 in polynomial1.items():
        key1 = int(key1[2:])  # removing x^ part
        for key2, value2 in polynomial2.items():
            key2 = int(key2[2:])  # removing x^ part
            new_key = f"x^{str(key1 + key2)}"
            new_value = value1 * value2
            if not new_key in resulting_polynomial:
                resulting_polynomial[new_key] = new_value
            else:
                resulting_polynomial[new_key] += new_value
    return resulting_polynomial


def add_polynomials(polynomial1: dict, polynomial2: dict):
    resulting_polynomial = {}
    for polynomial in [polynomial1, polynomial2]:
        for key, value in polynomial.items():
            if not key in resulting_polynomial:
                resulting_polynomial[key] = value
                continue
            resulting_polynomial[key] += value

    return resulting_polynomial


def from_dict_to_funct(polynomial: dict):
    returning_funct = 0
    for key, value in polynomial.items():
        value = round(value, 3)
        power = int(key[2:])
        if power == 0:
            returning_funct += value
            continue
        returning_funct += value*(x**power)
    return returning_funct


def Lagrange(X: list, Y: list = None, funct=None):
    X, Y = basic_checks(X, Y, funct)
    resulting_polynomial = {"x^0": 0}
    general_form = 0
    for i in range(len(X)):
        current_polynomial = None
        simplified_current_polynomial = {}
        current_divider = 1

        for j in X:
            if j == X[i]:
                continue
            if current_polynomial == None:
                current_polynomial = (x-j)
                simplified_current_polynomial = {"x^1": 1, "x^0": -j}
                continue
            current_polynomial = current_polynomial * (x-j)
            simplified_current_polynomial = multiply_polynomials(
                simplified_current_polynomial, {"x^1": 1, "x^0": -j})
        for k in X:
            if k == X[i]:
                continue
            current_divider = current_divider * (X[i]-k)

        general_form += round(
            Y[i]/current_divider, 3) * current_polynomial

        simplified_current_polynomial = multiply_polynomials({"x^0": round(
            Y[i]/current_divider, 3)}, simplified_current_polynomial)
        resulting_polynomial = add_polynomials(
            resulting_polynomial, simplified_current_polynomial)

    return general_form, from_dict_to_funct(resulting_polynomial)


def find_divided_difference(difference, X: list, Y: list, known_differences):
    if difference in known_differences:
        return known_differences[difference], known_differences

    buffer = difference
    difference = list(map(lambda a: int(a[1:]), difference[2:-1].split(",")))
    if (len(difference) == 1):
        return Y[difference[0]], known_differences
    difference_denominator = X[difference[-1]] - X[difference[0]]

    difference_numerator1 = "f("
    for i in range(difference[1], difference[-1]+1):
        difference_numerator1 += f"x{i}"
        if i != (difference[-1]):
            difference_numerator1 += ","
    difference_numerator1 += ")"

    if (difference_numerator1 in known_differences):
        difference_numerator1 = known_differences[difference_numerator1]
    else:
        difference_numerator1, known_differences = find_divided_difference(
            difference_numerator1, X, Y, known_differences)

    difference_numerator2 = "f("
    for i in range(difference[0], difference[-1]):
        difference_numerator2 += f"x{i}"
        if i != (difference[-1]-1):
            difference_numerator2 += ","
    difference_numerator2 += ")"

    if (difference_numerator2 in known_differences):
        difference_numerator2 = known_differences[difference_numerator2]
    else:
        difference_numerator2, known_differences = find_divided_difference(
            difference_numerator2, X, Y, known_differences)

    difference = round(
        (difference_numerator1 - difference_numerator2) / difference_denominator, 3)
    known_differences[buffer] = difference

    return difference, known_differences


def Newton(X: list, Y: list = None, funct=None):
    X, Y = basic_checks(X, Y, funct)
    found_divided_differences = {}
    returning_polynomial = {"x^0": 0}
    general_form = 0

    for i in range(len(X)):
        current_polynomial = None
        simplified_current_polynomial = None
        current_divided_difference = "f("
        for j in range(0, i+1):
            current_divided_difference += "x" + str(j)
            if j != i:
                current_divided_difference += ","
        current_divided_difference += ")"
        current_divided_difference, found_divided_differences = find_divided_difference(
            current_divided_difference, X, Y, found_divided_differences)
        for k in range(0, i+1):
            index = k
            k = X[k]
            if current_polynomial == None and i != 0:
                current_polynomial = (x-k)
                simplified_current_polynomial = {"x^1": 1, "x^0": -k}
                continue
            if current_polynomial == None and i == 0:
                current_polynomial = 1
                simplified_current_polynomial = {"x^0": 1}
                break
            if index == i:
                continue

            current_polynomial = current_polynomial*(x-k)
            simplified_current_polynomial = multiply_polynomials(
                simplified_current_polynomial, {"x^1": 1, "x^0": -k})

        general_form += current_divided_difference * current_polynomial
        simplified_current_polynomial = multiply_polynomials(
            simplified_current_polynomial, {"x^0": current_divided_difference})
        returning_polynomial = add_polynomials(
            returning_polynomial, simplified_current_polynomial)

    return general_form, from_dict_to_funct(returning_polynomial)


def find_interpolation_error(interpolation_polynomial, X: list, Xstar: float, funct) -> tuple[float, float]:
    derivative = funct
    for i in range(len(X)):
        derivative = diff(derivative, x)

    # maximumDer = maximum(Abs(derivative), x, Interval(X[0], X[-1]))

    critical_points = solve(derivative)
    values = []
    for point in critical_points:
        if point <= X[-1] and point >= X[0]:
            values.append(Abs(derivative.subs(x, point)))
    values += [Abs(derivative.subs(x, X[0])), Abs(derivative.subs(x, X[-1]))]

    maximumDer = max(values)

    denominator = factorial(len(X))
    polynomial = None
    for i in X:
        if polynomial == None:
            polynomial = (Xstar-i)
            continue
        polynomial = polynomial*(Xstar-i)

    actual_error = abs(funct.subs(x, Xstar) -
                       interpolation_polynomial.subs(x, Xstar))
    return round(((maximumDer/denominator)*polynomial), 1), round(actual_error, 3)


if __name__ == "__main__":
    funct = 1/x + x
    X1 = [0.1, 0.5, 0.9, 1.3]
    X2 = [0.1, 0.5, 1.1, 1.3]
    Y1 = [10.1, 2.5, 18.1/9, 2.69/1.3]
    Xstar = 0.8

    print(
        f"\nІнтерполяційні поліноми, розраховані для:\nX:{X1}\nY:{Y1}\n{''.ljust(50, '_')}")

    general_lagrange1, L1 = Lagrange(X1, Y1)
    print(
        f"Загальна форма многочлена Лагранжа:\n{general_lagrange1}\nІнтерполяційний поліном Лагранжа:\n")
    pprint(L1)

    general_newton1, N1 = Newton(X1, Y1)
    print(
        f"Загальна форма многочлена Ньютона:\n{general_newton1}\nІнтерполяційний поліном Ньютона:\n")
    pprint(N1)

    max_error, actual_error = find_interpolation_error(L1, X1, Xstar, funct)

    print(
        f"Інтерполяційна похибка:\n|f({Xstar})-L({Xstar})|<={max_error}\n|f({Xstar})-L({Xstar})|={actual_error}")

    print(
        f"\n\n\nІнтерполяційні поліноми, розраховані для:\nX:{X1}\nf(x):")
    pprint(funct)
    print(f"\n{''.ljust(50, '_')}")

    general_lagrange2, L2 = Lagrange(X2, funct=funct)
    print(
        f"Загальна форма многочлена Лагранжа:\n{general_lagrange2}\nІнтерполяційний поліном Лагранжа:\n")
    pprint(L2)

    general_newton2, N2 = Newton(X2, funct=funct)
    print(
        f"Загальна форма многочлена Ньютона:\n{general_newton2}\nІнтерполяційний поліном Ньютона:\n")
    pprint(N2)

    max_error, actual_error = find_interpolation_error(N2, X2, Xstar, funct)

    print(
        f"Інтерполяційна похибка:\n|f({Xstar})-L({Xstar})|<={max_error}\n|f({Xstar})-L({Xstar})|={actual_error}\n\n")

    plot(funct, L1, N2, xlim=[0, 1.5], ylim=[0, 12], legend=True)
