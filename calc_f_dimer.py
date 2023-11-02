#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 16:17:19 2023

@author: justin


See Chadda & Robertson, 2016. "Measuring membrane protein dimerization equilibrium..."
Section 5.2 (equation 8)

"""

import numpy as np
from scipy.optimize import minimize

def calculate_R2(FDimer, pl, pe, pm, pd):
    # Calculate R2 for the given FDimer
    R2 = 0
    for n in range(1, 4):
        R2 += (pe[n-1] - ((1 - FDimer) * pm[n-1] + FDimer * pd[n-1]))**2
    return R2

def fit_FDimer(pl, pe, pm, pd):
    # Define the objective function for least-squares minimization
    def objective_func(FDimer):
        return calculate_R2(FDimer, pl, pe, pm, pd)

    # Use minimize function to find the value of FDimer that minimizes R2
    result = minimize(objective_func, x0=0.5, bounds=[(0, 1)], method='SLSQP')  # Initial guess: 0.5

    return result.x[0]
    


def main():
    # User input
    peptide_concentrations = np.array([1.5, 0.75, 0.46875, 0.375, 0.1875, 0.09375])*1e-6
    lipid_concentration = 200e-6
    pl = peptide_concentrations/lipid_concentration
   
    monomer = np.array([24.45141066, 45.80645161, 61.49425287, 68.6, 72.92225201])
    dimer = np.array([51.41065831, 41.61290323, 31.6091954, 26.4, 23.05630027])
    trimer_plus = np.array([24.13793103, 12.58064516, 6.896551724, 5, 4.021447721])
    
    # Create a 2D array 'pe' by stacking the individual arrays
    pe = np.vstack((monomer, dimer, trimer_plus)).T.tolist()
    
    # Convert to a list of lists
    pe = [list(pe_i) for pe_i in pe]
    
    # Now 'pe' is a list of lists with the desired structure
    print(pe)
    
    # These values were generated from sim_occ_over_pls
    # List of p1, p2, p3+ for noninteracting monomers
    pm = [[0.3151, 0.5246, 0.27280000000000004], 
          [0.3509, 0.3151, 0.0701], 
          [0.3097, 0.1929, 0.022699999999999998], 
          [0.2701, 0.1506, 0.0145], 
          [0.1729, 0.0664, 0.0021],
          [0.0965, 0.0296, 0.0003]]
    
    # List of p1, p2, p3+ for noninteracting dimers
    pd = [[0.0004, 0.0006, 1.997999999999991],
          [0.0009, 0.0023, 1.994299999999997],
          [0.0025, 0.0057, 1.9864999999999986],
          [0.0036, 0.0091, 1.9788000000000006],
          [0.0179, 0.0409, 1.9054999999999993],
          [0.1128, 0.2557, 1.4496]]  # List of pd values for different pl values

    # Example usage to fit FDimer as a function of pl
    FDimer_values = []
    for pe_i, pm_i, pd_i in zip(pe, pm, pd):
        FDimer_value = fit_FDimer(pl, pe_i, pm_i, pd_i)
        FDimer_values.append(FDimer_value)
    
    print("FDimer values:", FDimer_values)



    
    
    
    
    pass

if __name__ == "__main__":
    main()

