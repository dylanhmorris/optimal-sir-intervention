#!/usr/bin/env python3

####################################################
# filename: InterventionSIR.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: specifies the SIR models
# with and without or peak-minimizing
# intervention according to our
# model of intervention functions b(t)
####################################################

import numpy as np 
from scipy.integrate import odeint

import optimize_interventions as oi

class InterventionSIR():
    """
    class for handling our SIR
    model with interventions
    """

    def __init__(self,
                 b_func = None,
                 R0 = None,
                 gamma = None,
                 inits = None
    ):
        if b_func is None:
            b_func = Intervention()
        self.b_func = b_func
        self.R0 = R0
        self.gamma = gamma
        self.inits = inits
        self.reset()
        
    def reset(self):
        self.state = self.inits
        self.time = 0
        self.time_ts = np.array([])
        self.state_ts = np.array([[], [], []]).reshape((-1, 3))

    def deriv(self, state, time):
        S, I, R = state
        beta = self.R0 * self.gamma
        b = self.b_func(time, beta, self.gamma, S, I)
        
        dS = -b * beta * S * I
        dI = b * beta * S * I - self.gamma * I
        dR = self.gamma * I

        return np.array([dS, dI, dR])

    def deriv_null(self, state, time):
        S, I, R = state
        beta = self.R0 * self.gamma
        dS = -beta * S * I
        dI = beta * S * I - self.gamma * I
        dR = self.gamma * I
        return np.array([dS, dI, dR])
    
    def integrate(self, final_time, fineness = 10000):
        times = np.linspace(self.time,
                            final_time,
                            fineness)
        results = odeint(self.deriv, self.state, times)
        self.state = results[-1]
        self.time = final_time
        self.time_ts = np.concatenate([self.time_ts,
                                       times])
        self.state_ts = np.concatenate([self.state_ts,
                                        results])

        return (times, results)

    def integrate_null(self, final_time, fineness = 10000):
        times = np.linspace(self.time,
                            final_time,
                            fineness)
        results = odeint(self.deriv_null, self.state, times)
        return (times, results)
    
    def I_max_SI(self, S_x, I_x):
        """
        get the maximum value of I(t) 
        in the window from t s.t. S = S_x,
        I = I_x to t = infinity
        """
        return (S_x + I_x - 
                (1/self.R0) * np.log(S_x) - 
                (1/self.R0) + 
                (1/self.R0) * np.log(1/self.R0))

    def I_of_S(self, S):
        S0, I0, Rec0 = self.inits
        return (I0 + S0 - (1/self.R0) * np.log(S0) -
                S + (1/self.R0) * np.log(S))

    def t_of_S(self, S_target):
        S0, I0, Rec0 = self.inits
        if np.isnan(S_target):
            raise ValueError("Cannot find time "
                             "for non-numeric/nan S\n\n"
                             "check that S is being "
                             "calculated correctly")
        def deriv(t, S_val):
            I = self.I_of_S(S_val)
            return -1 / (self.R0 * self.gamma * S_val * I)
        return odeint(deriv, 0, np.linspace(S0, S_target, 2))[-1]

    def get_I_max(self,
                  allow_boundary_max = False):
        last_timestep_error = (
            "Max at last timestep. "
            "You likely need to "
            "increase integration "
            "max time. If this was expected, "
            "set allow_boundary_max = True")
        first_timestep_error = (
            "Max at first timestep. "
            "Your model may be misspecified. "
            "If this was expected, "
            "set allow_boundary_max = True")
        ## check that we didn't get a boundary soln
        wheremax = np.argmax(self.state_ts[:, 1])
        if not allow_boundary_max:
            if wheremax == self.state_ts[:, 1].size: 
                raise ValueError(last_timestep_error)
            elif wheremax == 0:
                raise ValueError(first_timestep_error)
        return self.state_ts[wheremax, 1]
    
    def get_t_peak(self):
        return self.t_of_S(1 / self.R0)

    def __repr__(self):
        return ("InterventionSIR with R0 = {}, "
                "gamma = {}, and an intervention "
                "function {}".format(
                    self.R0,
                    self.gamma,
                    self.b_func))



class Intervention():
    """
    class for defining intervention
    functions b(t)
    """

    def __init__(self,
                 tau = None,
                 t_i = None,
                 sigma = None,
                 f = None,
                 S_i_expected = None,
                 I_i_expected = None,
                 strategy = None):
        self.tau = tau
        self.t_i = t_i
        self.sigma = sigma
        self.f = f
        self.S_i_expected = S_i_expected
        self.I_i_expected = I_i_expected
        self.strategy = strategy
    
        self.repertoire = {
            "fixed": self.fixed_b,
            "maintain-suppress-time": self.maintain_suppress_time,
            "maintain-suppress-state": self.maintain_suppress_state,
            "full-suppression": self.fixed_b}

    def __call__(self,
                 time,
                 beta,
                 gamma,
                 S,
                 I):
        return self.repertoire[self.strategy](
            time,
            beta,
            gamma,
            S,
            I)

    def fixed_b(self,
                time,
                beta,
                gamma,
                S,
                I):
        """
        Fixed intervention of strictness
        sigma
        """
        if time >= self.t_i and time < self.t_i + self.tau:
            result = self.sigma
        else:
            result = 1
        return result

    def maintain_suppress_time(self,
                              time,
                              beta,
                              gamma,
                              S,
                              I):
        """
        Variable maintain/suppress
        intervention tuned by 
        current time
        """
        if time >= self.t_i and time < self.t_i + self.tau * self.f:
            S_expected = (self.S_i_expected -
                          gamma * (time - self.t_i) *
                          self.I_i_expected)
            result = gamma / (beta * S_expected)
        elif (time >= self.t_i + self.tau * self.f and
              time < self.t_i + self.tau):
            result = 0
        else:
            result = 1
        return result

    def maintain_suppress_state(self,
                               time,
                               beta,
                               gamma,
                               S,
                               I):
        """
        Variable maintain/suppress
        intervention tuned by 
        current state of the system
        (S(t), I(t))
        """
        if time >= self.t_i and time < self.t_i + self.tau * self.f:
            result = gamma / (beta * S)
        elif (time >= self.t_i + self.tau * self.f and
              time < self.t_i + self.tau):
            result = 0
        else:
            result = 1
        return result
        



## helper functions for the above:
def get_SI_expected(
        t_i_expected,
        R0,
        gamma,
        I0 = 1e-6,
        Rec0 = 0,
        integration_fineness = 1000000):
    """
    What S_i and I_i do we
    expect in a time tuned
    intervention that we
    plan to start at a time
    t_i_expected?
    """
    S0 = 1 - I0 - Rec0
    def deriv(state, time):
        beta = R0 * gamma
        S, I = state
        dS = -beta * S * I
        dI = beta * S * I - gamma * I
        return np.array([dS, dI])

    state = odeint(deriv, [S0, I0],
                   np.linspace(0, t_i_expected, integration_fineness))
    expected_S_i, expected_I_i = state[-1]
    return (expected_S_i, expected_I_i)

def make_state_tuned_variable_b_func(tau, t_i, f):
    """
    create a function to execute
    the variable-b intervention
    with parameters t_i, tau and f
    where we operate based on the current
    state
    """
    return Intervention(
        tau = tau,
        t_i = t_i,
        f = f,
        strategy = "maintain-suppress-state")

def make_time_tuned_variable_b_func(tau,
                                    t_i,
                                    f,
                                    S_i_expected,
                                    I_i_expected):
    """
    create a function to execute
    the variable-b intervention
    with parameters t_i, tau and f
    where we operate based on the current
    time
    """                                          
    return Intervention(
        tau = tau,
        f = f,
        t_i = t_i,
        S_i_expected = S_i_expected,
        I_i_expected = I_i_expected,
        strategy = "maintain-suppress-time")



def make_fixed_b_func(tau, t_i, sigma):
    """
    create a function to execute
    the fixed intervention
    with parameters t_i, tau and sigma 
    """
    return Intervention(
        tau = tau,
        t_i = t_i,
        sigma = sigma,
        strategy = "fixed")


def check_gamma(gamma,
                model,
                offset,
                strategy = None,
                R0 = None,
                tau = None,
                verbose = False):
    """
    Check that a gamma 
    value produces a viable
    t_i given the model
    and its intervention
    function
    """
    if strategy is None:
        strategy = model.b_func.strategy
    if R0 is None:
        R0 = model.R0
    if tau is None:
        tau = model.b_func.tau

    model.gamma = gamma
        
    if strategy in ["maintain-suppress-time", "maintain-suppress-state"]:
        opt_S_i, opt_f = oi.calc_Sf_opt(
            R0,
            gamma * tau)
        t_i = model.t_of_S(opt_S_i)[0]
    elif strategy in ["fixed"]:
        fixed_S_i, fixed_sigma = oi.calc_Sb_opt(
            R0,
            gamma,
            tau)
        t_i = model.t_of_S(fixed_S_i)[0]
    elif strategy in ["full-suppression"]:
        full_supp_S_i = oi.calc_S_var_opt(
            R0,
            gamma * tau,
            0)
        t_i = model.t_of_S(full_supp_S_i)[0]
    if verbose:
        print("strategy:", strategy)
        print("gamma:", gamma)
        print("t_i:", t_i)
        print("offset:", offset)
    return min(t_i + offset, t_i) > 0

