#!/usr/bin/env python3

####################################################
# filename: optimize_interventions.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: helper code for optimizing
# intervention parameters
####################################################

from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
from scipy.optimize import root
import numpy as np
from scipy.special import lambertw
from scipy.integrate import odeint

def I_of_S(S, R0, S0 = 1, I0 = 0):
    return (I0 + S0 - (1/R0) * np.log(S0) -
            S + (1/R0) * np.log(S))

def S_of_I(I, R0, S0 = 1, I0 = 0):

    Ipeak = I_of_S(1/R0, R0, S0 = S0, I0 = I0)
    
    if I > Ipeak:
        raise ValueError("I must be smaller "
                         "than I_peak = {}".format(Ipeak))
    exponent = R0 * (I - I0 - S0)
    result = np.real(-lambertw(-S0 * np.exp(exponent) * R0) / R0)
    return result

def S_of_t(t, R0, gamma, fineness = 100000, I0 = 1e-6):
    times = np.linspace(0, t, fineness)
    def deriv(state, time):
        S, I = state
        return np.array([-R0 * gamma * S * I,
                         R0 * gamma * S * I - gamma * I])
    return (odeint(deriv, [1 - I0, I0],
                   times)[-1:,0])

def I_max(S, I, R0):
    return (S + I - (1/R0) * np.log(S) -
            (1/R0) + (1/R0) * np.log(1/R0))

def I_max_after_opt_intervention(S_i, I_i, f, R0, gammatau):
    S_f = S_i - I_i * gammatau * f
    I_f = I_i * np.exp(-gammatau * (1 - f))
    return I_max(S_f, I_f, R0)
    
def I_max_opt_of_S_i(S_i, f, R0, gammatau):
    I_i = I_of_S(S_i, R0)
    return max(
        I_max_after_opt_intervention(S_i, I_i, f, R0,
                                     gammatau),
        I_i)


def log_Imax_opt_S(S, f, R0, gammatau):
    I = I_of_S(S, R0)
    return np.log(I_max_opt_of_S_i(S, f, R0, gammatau) + 1)

def log_Imax_opt_f(f, S, R0, gammatau):
    I = I_of_S(S, R0)
    return np.log(I_max_opt_of_S_i(S, f, R0, gammatau) + 1)

def log_Imax_vec(Sf_vec, R0, gammatau):
    S, f = Sf_vec
    return np.log(I_max_opt_of_S_i(S, f, R0, gammatau))


def constrain_Scrit(Sf_vec, R0, gammatau):
    S, f = Sf_vec
    I_i = I_of_S(S, R0)
    Scrit = (1/R0) + 1e-5
    return np.atleast_1d(S - gammatau * f * I_i - Scrit)

def constrain_I_i(Sf_vec, R0, gammatau):
    S, f = Sf_vec
    I_i = I_of_S(S, R0)
    return np.atleast_1d(I_max_opt_of_S_i(S, f, R0, gammatau) - I_i)


def calc_Sf_opt(R0,
                gammatau,
                method = None,
                func = log_Imax_vec):
    S_target = 1
    f_target = max(0, 1 - 1/gammatau)
    guess_vec = [S_target , f_target]
    minima = minimize(func, guess_vec, (R0, gammatau),
                      method = method,
                      constraints = [{"fun": constrain_Scrit,
                                      "type": "ineq",
                                      "args": (R0, gammatau)},
                                      {"fun": constrain_I_i,
                                       "type": "eq",
                                       "args": (R0, gammatau)}],
                      bounds = ([0, 1], [0, 1]))
    if minima.success:
        return minima.x
    else:
        return np.array([np.nan, np.nan])


def tau_crash(R0):

    def to_solve(x):
        return x * (np.log(x) + x * np.log(R0 * x) - R0 * x) + 1
    
    x = root(to_solve, 0.5)
    if x.success:
        return -np.log(1 - float(x.x))
    else:
        return np.nan

def full_sup_asymptote(R0):
    return 0.5 + (1 / (R0 * 2)) * (np.log(1/R0) - 1)

def calc_Sf_opt_brute(
        R0,
        gammatau,
        f_guess_init = 0.5,
        n_refinements = 5):

    f_guess = f_guess_init
    f_tol = 1
    for i_ref in range(n_refinements):
        fs = np.linspace(max(0, f_guess - f_tol),
                         min(1, f_guess + f_tol),
                         100)
        vals = [I_max_opt_of_S_i(
            calc_S_var_opt(R0, gammatau, f),
            f, R0, gammatau)
                for f in fs]
        f_guess = fs[np.argmin(vals)]
        f_tol = fs[1] - fs[0]
    S_i = calc_S_var_opt(R0, gammatau, f_guess)
    return np.array([S_i, f_guess])


def calc_f_opt(R0,
               gammatau,
               gamma = None,
               S_i = None,
               t_i = None,
               method = "bounded"):
    if t_i is not None and gamma is not None:
        S_i = S_of_t(S, R0, gamma)
    elif S_i is None:
        raise ValueError("Must provide either S_i or "
                         "t_i and gamma")
    I_i = I_of_S(S_i, R0)
    if S_i <= 1/R0:
        raise ValueError("S must be greater than "
                         "Scrit = 1/R0")
    max_f = (S_i - 1/R0) / (gammatau * I_i)
    minima = minimize_scalar(log_Imax_opt,
                             args = (S_i, R0, gammatau),
                             method = method,
                             bounds = ([0, max_f]))
    if minima.success:
        return minima.x
    else:
        return np.nan

    I_i = I_of_S(S_i, R0)
    if S_i <= 1/R0:
        raise ValueError("S must be greater than "
                         "Scrit = 1/R0")
    max_f = (S_i - 1/R0) / (gammatau * I_i)
    minima = minimize_scalar(log_Imax_opt_f,
                             args = (S_i, R0, gammatau),
                             method = method,
                             bounds = ([0, max_f]))
    if minima.success:
        return minima.x
    else:
        return np.nan


def min_S(R0, f, gammatau):
    """
    largest possible S_i 
    that does not result
    in dipping below Scrit.
    """
    if f <= 0 or gammatau <= 0:
        return 1/R0
    fg = f * gammatau
    Rfgo = R0 * (fg + 1)
    candidate = np.real(
        -fg * lambertw(-np.exp(-(1/fg) - R0) * Rfgo / fg,
                               k = -1) / Rfgo)
    return max(1/R0, candidate)


def calc_S_var_opt(R0,
                   gammatau,
                   f,
                   method = "bounded"):

    min_S_val = min_S(R0, f, gammatau)

    minima = minimize_scalar(log_Imax_opt_S,
                             args = (f, R0, gammatau),
                             method = method,
                             bounds = ([min_S_val, 1]))
    if minima.success:
        return minima.x
    else:
        return np.nan

    I_i = I_of_S(S_i, R0)
    if S_i <= 1/R0:
        raise ValueError("S must be greater than "
                         "Scrit = 1/R0")
    max_f = (S_i - 1/R0) / (gammatau * I_i)
    minima = minimize_scalar(log_Imax_opt,
                             args = (S_i, R0, gammatau),
                             method = method,
                             bounds = ([0, max_f]))
    if minima.success:
        return minima.x
    else:
        return np.nan


def t_of_S(S, R0, gamma, I0 = 1e-6, Rec0 = 0):
    S0 = 1 - I0 - Rec0
    def deriv(t, S_val):
        I = I_of_S(S_val, R0, S0 = S0, I0 = I0)
        return -1 / (R0 * gamma * S_val * I)
    return odeint(deriv, 0, np.linspace(S0, S, 2))[1]


def Imax_of_S_i_b(S_i,
                  b,
                  R0,
                  gamma,
                  tau):
    I_i = I_of_S(S_i, R0)

    def deriv(state, time):
        beta = R0 * gamma * b
        S, I = state
        dS = -beta * S * I
        dI = beta * S * I - gamma * I
        return np.array([dS, dI])
    
    intervention = odeint(deriv, [S_i, I_i],
                          np.linspace(0, tau,
                                      (1 + int(tau)) * 1000))
    I_max_interv = np.max(intervention[:, 1])
    S_f, I_f = intervention[-1]
    I_max_f = I_max(S_f, I_f, R0)
    return np.max([I_i, I_max_f, I_max_interv])

def calc_Sb_opt(R0,
                gamma,
                tau,
                verbose = False,
                method = None,
                S_guess = 0.97,
                b_guess = 0.5,
                n_max_tries = 10):
    raw_peak = I_max(1, 0, R0)
        
    def func(Sb_vec):
        S_i, b = Sb_vec
        return Imax_of_S_i_b(S_i,
                             b,
                             R0,
                             gamma,
                             tau)
    Imax = 1
    result = np.array([np.nan, np.nan])
    for k in range(n_max_tries):
        b_guess = max(0, min(np.random.normal(b_guess, 0.05), 1))
        guess_vec = np.array(
            [S_guess, b_guess])
        minima = minimize(func, guess_vec,
                          method = method,
                          bounds = ([1/R0, 1],
                                    [0, 1]))
        if minima.success:
            if minima.fun < Imax:
                result = minima.x
                Imax = minima.fun
    return result
