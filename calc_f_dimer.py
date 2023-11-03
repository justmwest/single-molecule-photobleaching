#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 16:17:19 2023

@author: justin


See Chadda & Robertson, 2016. "Measuring membrane protein dimerization equilibrium..."
Section 5.2 (equation 8)

"""

import numpy as np
import matplotlib.pyplot as plt


def fit_frac_dimer(pe, pm, pd):
    """ 
    This function is made to calculate the fraction of dimer
    by minimization of the sum of squared residuals between simulated
    ideal monomers and dimers 
    Adapted from Chadda & Robertson JGP 2017 MatLab app.
    
    pe, pm, pd: each is a list of 3 probabilities in order. [p1, p2, p3+]
    pe = Experimentally determined probabilities
    pm = Simulated probabilities for noninteracting monomers
    pd = Simulated probabilities for noninteracting dimers
    """
    # Create an array of dimer fractions (0 to 1) with a step of 0.01
    FD = np.arange(0, 1.01, 0.01)
    nFD = FD.shape  # Get the size of the FD array
    
    N = 3  # number of data points
    
    # Calculation of the sum of squared residuals (R2)
    R2 = []
    for i in range(nFD[0]):
        total_sum = 0
        for k in range(N):
            Pfit = (1 - FD[i]) * pm[k] + FD[i] * pd[k]
            total_sum += (pe[k] - Pfit) ** 2
        R2.append(total_sum)
    
    # Plotting the results
    plt.semilogy(FD, R2)
    plt.xlabel('Dimer Fraction (FD)')
    plt.ylabel('Sum of Squared Residuals (R2)')
    plt.title('Least-Squares Analysis')
    plt.show()
    
    # Minimum R2 value and Fdimer prediction
    min_R2 = min(R2)
    min_R2_index = R2.index(min_R2)
    Fdimer_prediction = FD[min_R2_index]
    
    print(f"Minimum R2 value: {min_R2}")
    print(f"Fdimer prediction: {Fdimer_prediction}")

    



def main():
    # User input
    peptide_concentrations = np.array([1.5, 0.75, 0.46875, 0.375, 0.1875, 0.09375])*1e-6
    lipid_concentration = 75e-6
    pl = peptide_concentrations/lipid_concentration
   
    # List of p1, p2, p3+ from experimental data
    # No value for peptide concentration 0.09375e-6
    pe = [[0.024451410660000002, 0.05141065831, 0.02413793103], 
          [0.04580645161000001, 0.041612903230000005, 0.012580645159999999], 
          [0.061494252869999995, 0.0316091954, 0.0068965517240000005], 
          [0.0686, 0.0264, 0.005], 
          [0.07292225201, 0.023056300270000003, 0.004021447721]]
    
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

    fit_frac_dimer(pe[1], pm[1], pd[1])

    pe_1=[sublist[0] for sublist in pe]
    pe_2=[sublist[1] for sublist in pe]
    pe_3=[sublist[2] for sublist in pe]
    
    print(pe_1)
    print(pe_2)
    print(pe_3)
    
    
    
    pass

if __name__ == "__main__":
    main()

