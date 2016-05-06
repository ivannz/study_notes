#
# coding: utf-8
# Copyright (C) DATADVANCE, 2010-2015
#

import numpy as np
from numpy import zeros, ones, pi, tanh, exp, sin, cos, sqrt, arctan, absolute

# fucnlist = ['f1', 'f2', 'f3', 'f4', 'f5', 'branin', 'mystery', 'michalewicz', 'sqmichalewicz', 'rosenbrock', 'ellipsoidal',
#             'rastrigin', 'schwefel', 'himmelblau', 'six_hump_camel_back', 'GoldsteinPrice', 'Kink', 'CurvedKink', 'StraightDiscont',
#             'RoundDiscont', 'TwoHumps', 'aero1', 'aero2', 'aero3']

def heaviside(X):
    return (X >= 0).astype(np.float)

def f1(X):
    t = np.sum(X**2, axis=1) - 0.3
    Y = t**2 / (t**2 + 0.01)
    return Y

def f2(X):
    center1 = 0.5 + zeros((len(X), 2))
    center2 = zeros((len(X), 2))
    Y = tanh(0.3 * (5 * exp(- 4 * np.sum((X - 0.5 * center1)**2, axis=1)) - 7 * exp(- 40 * np.sum((X - center2)**2, axis=1))
                    - 3 * exp(- 10 * np.sum((X[:, 0] + 0.7)**2)) + np.sum(X, axis=1)))
    return Y

def f3(X):
    t = np.sum((X + 0.6)**2, axis=1) - 0.3
    Y = sin(t)**2 / tanh(t**2 + 0.4)
    return Y

def f4(X):
    Y = np.sum(X, axis=1) / (1 + 4 * np.sum(X**2, axis=1))
    return Y

def f5(X):
    center1 = 0.5 + zeros((len(X), 2))
    gamma = ones((len(X), 2))
    Y = (np.sum(gamma * X, axis=1) + 1 * (np.sum((X - 0 * center1)**2, axis=1) <= 0.5**2) -
         2 * (np.sum((X - 0.7)**2, axis=1) <= 1**2))
    return Y

def branin(X):
    X = (X + 1) / 2.0
    # rescale to  - 5<x1<10, 0<x2<15
    x1 = 15 * X[:, 0] - 5
    x2 = 15 * X[:, 1]
    Y = (
        x2 - 5.1 / 4 / pi**2 * x1**2 + 5 / pi * x1 - 6)**2 + 10 * (
            1 - 1 / 8 / pi) * cos(
                x1) + 10
    return Y

def mystery(X):
    X = (X + 1) / 2.0
    # rescale to 0<x1<5, 0<x2<5
    x1 = 5 * X[:, 0]
    x2 = 5 * X[:, 1]
    Y = 2 + 0.01 * (
        x2 - x1**2)**2 + (
            1 - x1)**2 + 2 * (
                2 - x2)**2 + 7 * sin(
                    0.5 * x1) * sin(
                        0.7 * x1 * x2)
    return Y

def michalewicz(X):
    X = (X + 1)/2
    # rescale to 0<x1<pi, 0<x2<pi
    x1 = pi * X[:, 0]
    x2 = pi * X[:, 1]
    Y = sin(x1) * sin(x1**2 / pi) + sin(x2) * sin(2 * x2**2 / pi)
    return Y

def sqmichalewicz(X):
    X = (X + 1) / 2.0
    # rescale to 0<x1<pi, 0<x2<pi
    x1 = pi * X[:, 0]
    x2 = pi * X[:, 1]
    Y = (sin(x1) * sin(x1**2 / pi) + sin(x2) * sin(2 * x2**2 / pi))**2
    return Y

def rosenbrock(X):
    X = (X + 1) / 2.0
    # rescale to  - 2.048<x<2.048
    c = 2.048
    x1 = c * (2 * X[:, 0] - 1)
    x2 = c * (2 * X[:, 1] - 1)
    Y = 100 * (x2 - x1**2)**2 + (1 - x1)**2
    return Y

def ellipsoidal(X):
    X = (X + 1) / 2.0
    # rescale to  - 1<x<1
    x1 = 2 * X[:, 0] - 1
    x2 = 2 * X[:, 1] - 1
    Y = x1**2 + 2 * x2**2
    return Y

def rastrigin(X):
    X = (X + 1) / 2.0
    # rescale to  - 5.12<x<5.12
    c = 5.12
    x1 = c * (2 * X[:, 0] - 1)
    x2 = c * (2 * X[:, 1] - 1)
    Y = 20 + (x1 - 10 * cos(2 * pi * x1)) + (x2 - 10 * cos(2 * pi * x2))
    return Y

def schwefel(X):
    X = (X + 1) / 2.0
    # rescale to  - 500<x<500
    c = 500
    x1 = c * (2 * X[:, 0] - 1)
    x2 = c * (2 * X[:, 1] - 1)
    Y = x1 * sin(sqrt(np.absolute(x1))) - x2 * sin(sqrt(np.absolute(x2)))
    return Y

def himmelblau(X):
    X = (X + 1) / 2.0
    # rescale to  - 6<x<6
    c = 6
    x1 = c * (2 * X[:, 0] - 1)
    x2 = c * (2 * X[:, 1] - 1)
    Y = (x1**2 + x2 - 11)**2 + (x1 + x2**2 - 7)**2
    return Y

def six_hump_camel_back(X):
    X = (X + 1) / 2.0
    # rescale to  - 3<x1<3,  - 2<x2<2
    x1 = 3 * (2 * X[:, 0] - 1)
    x2 = 2 * (2 * X[:, 1] - 1)
    a = 4.0
    b = 2.1
    c = 3.0
    Y = (a - b * (x1**2) + (x1**4) / c) * (
        x1**2) + x1 * x2 + a * ((x2**2) - 1) * (x2**2)
    return Y

def GoldsteinPrice(X):
    X = (X + 1) / 2.0
    # rescale to  - 2<x<2
    c = 2
    x1 = c * (2 * X[:, 0] - 1)
    x2 = c * (2 * X[:, 1] - 1)
    Y = ((1 + (x1 + x2 + 1)**2 * (19 - 14 * x1 + 3 * x1**2 - 14 * x2 + 6 * x1 * x2 + 3 * x2**2)) *
        (30 + (2 * x1 - 3 * x2)**2 * (18 - 32 * x1 + 12 * x1**2 + 48 * x2 - 36 * x1 * x2 + 27 * x2**2)))
    return Y

def Kink(X):
    X = (X + 1) / 2.0
    Y = np.min(
        np.vstack((X[:,
                     0]**2 + X[:,
                               1]**2,
                   2 * X[:,
                         0]**2)).T,
        axis=1)
    return Y

def CurvedKink(X):
    X = (X + 1) / 2.0
    Y = np.min(
        np.vstack((X[:,
                     0]**2 + X[:,
                               1]**2,
                   0.3 * ones(X.shape[0]))).T,
        axis=1)
    return Y

def StraightDiscont(X):
    X = (X + 1) / 2.0
    Y = heaviside(X[:, 0] - X[:, 1] - 0.0000001)
    return Y

def RoundDiscont(X):
    X = (X + 1) / 2.0
    Y = heaviside(1 / 7.0 - (X[:, 0] - 0.5)**2 - (X[:, 1] - 0.5)**2)
    return Y

def TwoHumps(X):
    X = (X + 1) / 2.0
    Y = exp(
        - 10 * ((X[:,
                   0] - 1 / 4)**2 + (X[:,
                                       1] - 1 / 4)**2)) + 2 * exp(-20 * ((X[:,
                                                                            0] - 3/4)**2 + (X[:,
                                                                                              1] - 1 / 2)**2))
    return Y

def aero1(X):
    X = (X + 1) / 2.0
    x1 = X[:, 0]
    x2 = X[:, 1]
    Y = arctan(
        x1**(
            1 / 2)) * (
                - arctan(
                    50 * (
                        x2 + 0.2)**4 * (
                            x1 - 0.3 - 0.2 * x2)) + 2 * x1 * (
                                1 - x1) + (
                                    1 + x2) * x2)
    return Y

def aero2(X):
    X = (X + 1) / 2.0
    x1 = X[:, 0]
    x2 = X[:, 1]
    Y = (arctan(x1**(1 / 2)) * (- arctan(80 * (x2 + 0.2)**3 * (x1 - 0.3 - 0.3 * x2)) + 2 * x1 * (1 - x1) +
        (1 + x2) * x2)) * (1 + 0.2 * sin(5 * pi * (1 - x1) * x2) + 0.2 * sin(x1)) * (1 + 0.5 * arctan(4 * (x2 - 0.5)) - x2)
    return Y

def aero3(X):
    X = (X + 1) / 2.0
    x1 = X[:, 0]
    x2 = X[:, 1]
    x1 = x1 * (1 + 0.3 * sin(pi * (1 - x1) * x2))
    x2 = x2 * (1 + 0.15 * sin(pi * x1 * (1 - x2)))
    phi = 0.003 * pi * (0.1 + (x1 - 0.7)**2 + (x2 - 0.6)**2)**(- 2)
    x1 = 0.7 + (x1 - 0.7) * cos(phi) - (x2 - 0.6) * sin(phi)
    x2 = 0.6 + (x1 - 0.7) * sin(phi) + (x2 - 0.6) * cos(phi)
    Y = (1 + 0.2 * sin(2 * pi * x1 * x2)) * (x1 + 0.2 * sin(- pi * x2**6) + 0.4 * x1 * (1 - x2) * ((cos(x1**2 + x2**2 - 1))**100) +
                                             0.3 * exp(-30 * (absolute(x1 - 0.6))**1.9 - 150 * (x2 - 0.7)**2) +
                                             0.2 * exp(-50 * (x1 - 0.3)**2 - 200 * (absolute(x2 - 0.8))**2.4) +
                                             0.3 * exp(-500 * (x1 - 0.9)**2 - 400 * (x2 - 0.75)**2) +
                                             0.13 * exp(-700 * (x1 - 0.85)**2 - 500 * (x2 - 0.65)**2))
    return Y


def func2D():
    return {fn_.__name__: fn_ for fn_ in [
                heaviside, f1, f2, f3, f4, f5, branin, mystery,
                michalewicz, sqmichalewicz, rosenbrock, ellipsoidal,
                rastrigin, schwefel, himmelblau, six_hump_camel_back,
                GoldsteinPrice, Kink, CurvedKink, StraightDiscont,
                RoundDiscont, TwoHumps, aero1, aero2, aero3, ]}
