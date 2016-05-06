# coding: utf-8
# Copyright (C) DATADVANCE, 2010-2015
#
"""1D test functions"""
import numpy as np
from numpy import array, append, ones, newaxis, absolute
from numpy import sin, cos, exp, sqrt, pi, arctan, tanh

def classic(x):
    """`classic`"""
    return (6 * x - 2)**2 * sin(12 * x - 4)

def sin2pix(x):
    """`sin2pix`"""
    return sin(2 * pi * x)

def sin10pix(x):
    """`sin10pix`"""
    return sin(10 * pi * x)

def airfoil(x):
    """`airfoil`"""
    return sqrt(x) * (1 - x) * (1.2 - x)

def heaviside(x):
    """`heaviside`"""
    return (x > 0.5).astype(float)

def kink(x):
    """`kink`"""
    return exp(-4 * absolute(x - 0.8))

def pressure1(x):
    """`pressure1`"""
    return (x**(0.3) * (1 - x) - 0.05 * arctan(30 * (x - 0.05)) +
            0.15 * exp(-50 * (x - 0.45)**2) - 0.20 * exp(-70 * (x - 0.9)**2))

def pressure2(x):
    """`pressure2`"""
    x = x * ((3 - x) / 2.0 + 0.01 * sin(pi * x))
    return (x**(0.55) * (1 - 0.99 * x)**0.48 -
            0.11 / (1 + 100 * (x - 0.3)**2) - 0.02 * exp(-200 * (x - 0.5)**2))

def f1(x):
    """`f1`"""
    tau = 0.36
    return (-4 * sin(3 * x)**2 * (x < tau) +
            sin((1.5 - 4 * x)) * (x > tau) + (x > 0.7))

def f2(x):
    """`f2`"""
    return (x * (x <= 0.3) + (0.6 - x) * (x > 0.3) * (x <= 0.7) +
            1 * (x - 0.8) * (x > 0.7))

def f3(x):
    """`f3`"""
    return (exp(-30 * (x - 0.5)**2) - 2 * exp(-100 * (x - 0.7)**2) -
            10 * exp(-50 * (x - 0.3)**2))

def f4(x):
    """`f4`"""
    return tanh(cos(0.3 * x + 30 * x**2))

def f5(x):
    """`f5`"""
    return 1 / (0.1 + absolute(x - 0.5)**2)

def f6(x):
    """`f6`"""
    return (1 - 4 * x**2) * (x < 0.5) + 2 * (x - 0.5) * (x > 0.5)

def sigmoid(x):
    """`sigmoid`"""
    # weights = array([15, -7])
    # points = append(x[:, newaxis], ones((x.shape[0], 1)), axis=1)
    # return 1.0 / (1.0 + exp(-points.dot(weights)))
    return 1.0 / (1.0 + exp(7 - 15 * x))

def triangle(x, width=0.25):
    """`triangle`"""
    return np.maximum(1 + np.minimum(0.5 - x, x - 0.5) / width, 0)

def get_functions():
    """Returns a dictionary"""
    functions_ = [classic, sin2pix, sin10pix, airfoil, heaviside, triangle,
                  kink, pressure1, pressure2, f1, f2, f3, f4, f5, f6, sigmoid]
    return {fn.__name__: fn for fn in functions_}
