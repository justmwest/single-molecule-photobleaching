#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 13:08:07 2023

@author: justin

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import R

def f_dimer_model(pl, Keq):
    return ((1 + 4 * pl * Keq) - np.sqrt(1 + 8 * pl * Keq)) / (4 * pl * Keq)


def find_equilibrium_constant(pl_list, f_dimer_list, initial_guess=[0.003]):
    """
    Initial guess is for the Keq.
    """
    fit_params, fit_covariance = curve_fit(f_dimer_model, pl_list, f_dimer_list, p0=initial_guess)
    k_eq = fit_params[0]

    # Assess the success of the fit using the fit covariance
    if np.any(np.diag(fit_covariance) < 0):
        print("Warning: Poor fit detected. Covariance matrix indicates fit issues.")

    return k_eq, fit_params


def calc_delta_g(k_eq, temperature=298):
    """
    Temperature is in Kelvin.
    """
    gas_constant = 1.98720425864083e-3  # kcal K^-1 mol^-1
    delta_g = -gas_constant * temperature * np.log(k_eq)  # kcal mol^-1
    return delta_g

def plot_keq_fit(pl_list, f_dimer_list, fit_params, delta_g):

    plt.figure(figsize=(9, 6))
    plt.loglog(pl_list, f_dimer_list, 'ko', label='Data', markersize=8)
    plt.loglog(pl_list, f_dimer_model(pl_list, *fit_params), 'k-', label='Fit')

    # Set plot labels and style
    plt.xlabel('Peptide to Lipid Molar Ratio (P/L)')
    plt.ylabel('Fraction Dimer (calculated)')
    plt.legend()

    # Display ΔG value on the plot
    plt.text(max(pl_list)-0.01, min(f_dimer_list), f'ΔG = {delta_g:.2f} kcal/mol', fontsize=14)

    plt.show()



def main():
    # Inputs
    peptide_concentrations = np.array([1.5, 0.75, 0.46875, 0.375, 0.1875]) * 1e-6
    lipid_concentration = 75 * 1e-6  # M
    pls = peptide_concentrations / lipid_concentration

    # f_dimer data determined by fitting simulations to experimental data.
    # see fdimer.py. Below is synthetic data.
    f_dimer = np.array([.5, .4, .3, .25, .20])

    k_eq, fit_params = find_equilibrium_constant(pls, f_dimer)

    delta_g = calc_delta_g(k_eq)

    plot_keq_fit(pls, f_dimer, fit_params, delta_g)












if __name__ == "__main__":
    main()