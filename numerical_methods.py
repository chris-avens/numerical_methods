def trapezoidal_rule(f, a, b, n=100000):
    """
    Approximates the definite integral of f(x) over the interval [a, b] using the trapezoidal rule.

    :param f: the function to integrate
    :param a: the lower limit of integration
    :param b: the upper limit of integration
    :param n: the number of sub intervals to use
    :return: The approximate value of the definite integral of f(x) over the interval [a, b].
    """
    if b < a:
        (a, b) = (b, a)
    root = secant_method(f=f, x0=a, x1=b)
    if a < root < b:
        response = input(f"there's a root (zero) between the range of integration at x={root}. Continue? (y/n)")
        if response == "n":
            print("procedure disrupted. please enter an interval without root in between")
            return 0

    h = (b - a) / n
    s = (f(a) + f(b)) / 2
    for i in range(1, n):
        s += f(a + i*h)
    return h * s


def simpsons_rule(f, a, b, n=100):
    """
    Approximates the definite integral of f(x) over the interval [a, b] using Simpson's rule.

    :param f: the function to integrate
    :param a: the lower limit of integration
    :param b: the upper limit of integration
    :param n: the number of sub intervals to use (must be even)
    :return: The approximate value of the definite integral of f(x) over the interval [a, b].
    """
    if b < a:
        (a, b) = (b, a)
    root = secant_method(f=f, x0=a, x1=b)
    if a < root < b:
        response = input(f"there's a root (zero) between the range of integration at x={root}. Continue? (y/n)")
        if response == "n":
            print("procedure disrupted. please enter an interval without root in between")
            return 0

    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n, 2):
        s += 4 * f(a + i*h)
    for i in range(2, n-1, 2):
        s += 2 * f(a + i*h)
    return h * s / 3


def bisection(f, x0, x1, tol=0.0000000001):
    """
    Finds a root of the function f using the Bisection method.
    Newton-Raphson is the common go-to for finding root, but it fails if:
    - function is not differentiable in all values of x
    - function has vertical asymptote
    Bisection method can overcome this limitations

    :param f: function to be evaluated the roots
    :param x0: x value at the left of assumed root
    :param x1: x value at the right of assumed root
    :param tol: tolerance or accuracy
    :return: root of the function in between x0 and x1
    """
    # iter_count = 0
    while abs(x1 - x0) > tol:
        # iter_count += 1
        # proceeding in this line means that the interval is not close enough together based on tolerance
        x2 = (x0 + x1) / 2
        # print(f'number of iterations: {iter_count}')
        # print(f'x0: {x0} | y0: {f(x0)}')
        # print(f'x1: {x1} | y1: {f(x1)}')
        # print(f'x2: {x2} | y2: {f(x2)}')
        if f(x2) == 0:
            return x2
        elif f(x0) * f(x2) < 0:
            # if x0 and x2 have the same sign, the estimate won't approach to zero
            x1 = x2
        elif f(x1) * f(x2) < 0:
            # if x1 and x2 have the same sign, the estimate won't approach to zero
            x0 = x2
        elif f(x0) * f((x0 + x2) / 2) < 0:
            # second iteration already
            # just in rare case that both of condition above did not satisfy
            x1 = x2
        else:
            x0 = x2

    if abs(f(x2)) < tol:
        return x2
    else:
        # seems like it needs to try a different range
        # if the function value does not oscillate around zero
        print(f"""an error causes the calculation not to converge to zero:
    - bisection method may be inappropriate to use in the function and/or the interval given
    - try using a different interval
    - the function may not have any zeroes in the given interval")
    - the tolerance is very low (tol={tol})
    - check the value of the function in last iteration: f({x2})={f(x2)}""")


def secant_method(f, x0, x1, tol=0.0000000001, max_iter=1000):
    """
    Finds a root of the function f using the secant method.
    It converges quicker than a linear rate, making it more convergent than the bisection method
    It does not necessitate the usage of the function’s derivative, which is not available in a number of applications.
    Unlike Newton’s technique, which requires two function evaluations in every iteration, it only requires one.

    :param f: the function to find a root of
    :param x0: the first initial guess for the root
    :param x1: the second initial guess for the root
    :param tol: the tolerance for the root
    :param max_iter: the maximum number of iterations allowed
    :return: the estimated root of the function f
    """
    i = 0
    while i < max_iter:
        fx0 = f(x0)
        fx1 = f(x1)
        if abs(fx1 - fx0) < tol:
            return x1
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0)  # estimates the root denoted by x2
        x0 = x1
        x1 = x2
        i += 1
    raise ValueError("Failed to converge after %d iterations" % max_iter)


def newton_raphson(f, df, x0, tol=0.0000000001, max_iter=1000):
    """
    Finds a root of the function f using the Newton-Raphson method.

    :param f: the function to find a root of
    :param df: the derivative of f
    :param x0: the initial guess for the root
    :param tol: the tolerance for the root
    :param max_iter: the maximum number of iterations allowed
    :return: the estimated root of the function f
    """

    x = x0
    for iter_count in range(max_iter):
        fx = f(x)
        if abs(fx) < tol:
            return x
        """calculates the value of the derivative df at x, and check if it's zero.
        If it is, we can't continue with the Newton-Raphson method and raise a ValueError."""
        dfx = df(x)
        if dfx == 0:
            raise ValueError(f"Derivative is zero at x = %f" % x)
        x = x - fx / dfx
    raise ValueError("Failed to converge after %d iterations" % max_iter)


def runge_kutta(f, x0=0.0, y0=0.0, h=0.01, n=100):
    """
    Solves the differential equation y' = f(x, y) using the fourth-order Runge-Kutta method.

    :param f: the function defining the differential equation y' = f(x, y)
    :param x0: the initial value of x
    :param y0: the initial value of y at x0
    :param h: the step size
    :param n: the number of steps to take
    :return: an array of x-values and corresponding y-values
    """
    x = [x0]
    y = [y0]
    # print(f'x0: {x[0]} | y0: {y[0]}')
    for i in range(n):
        # print(f'step {i+1}:')
        k1 = f(x[i], y[i])
        k2 = f(x[i]+h/2, y[i]+h*k1/2)
        k3 = f(x[i]+h/2, y[i]+h*k2/2)
        k4 = f(x[i]+h, y[i]+h*k3)
        y_next = y[i] + h*(k1+2*k2+2*k3+k4)/6
        x.append(round(x[i] + h, 15))
        # rounding to 15 decimal places because sometimes there are tendencies to change from 0.1 to 0.1000000000000001
        y.append(y_next)
        """
        Just delete the triple quotes and this line if want to print values
        print(f'k1: {k1}')
        print(f'k2: {k2}')
        print(f'k3: {k3}')
        print(f'k4: {k4}')
        print(f'x{i+1}: {"{:.2f}".format(x[i+1])} | y{i+1}: {y[i+1]}')
        if (i+1) % 10 == 0:  # 10 divides 100 (default n) to 10 in order to print 10 values only
            print(y[i+1])
        """
    return x, y
