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


def least_squares_plot(FD, R2, Fdimer_prediction):
    """
    Plots the results of least-squares analysis.
    """
    plt.semilogy(FD, R2)
    plt.axvline(x=Fdimer_prediction, color='red', linestyle='--')
    inlay_text = Fdimer_prediction + 0.65 # Offset
    plt.text(inlay_text, min(R2), f'Predicted FD = {Fdimer_prediction}', color='red', horizontalalignment='left')
    plt.xlabel('Dimer Fraction (FD)')
    plt.ylabel('Sum of Squared Residuals (R2)')
    plt.title("Least-Squares Analysis")
    plt.show()

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
    
    least_squares_plot(FD, R2, Fdimer_prediction)
    
    print(f"Minimum R2 value: {min_R2}")
    print(f"Fdimer prediction: {Fdimer_prediction}")



def main():

   
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

if __name__ == "__main__":
    main()


# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.stats import chi2

# # Assuming you have the following data
# N = 3  # Replace with the actual number of data points
# PMX = np.random.rand(N)  # Replace with your actual monomer probability data
# PDX = np.random.rand(N)  # Replace with your actual dimer probability data
# Pexpt = np.random.rand(N)  # Replace with your actual experimental data

# # Define FD, I, and X (if not already defined)
# FD = 0.5  # Replace with the desired FD value
# I = 1  # Replace with the desired index
# X = np.arange(1, N + 1)

# # Plot resulting Pfit
# Pfit = [(1 - FD) * PMX[0] + FD * PDX[0], (1 - FD) * PMX[1] + FD * PDX[1], (1 - FD) * PMX[2] + FD * PDX[2]

# # Outputs Pfit
# # Replace the following lines with your actual output code (e.g., GUI updates or print statements)
# print(f"Pfit1 Value: {Pfit[0]}")
# print(f"Pfit2 Value: {Pfit[1]}")
# print(f"Pfit3 Value: {Pfit[2]}")

# # Plot Pfit
# plt.bar(X, Pfit)
# plt.xlabel('X')
# plt.ylabel('Pfit')
# plt.title('Pfit Distribution')
# plt.show()

# # Chi-squared analysis
# # H0: the expt and fitted distributions are the same
# # H1: the expt and fitted distributions are different
# chi2_value = 0
# for k in range(N):
#     chi2_value += ((Pexpt[k] - Pfit[k]) ** 2) / Pfit[k]

# # Calculate the p-value using the chi-squared distribution
# p_value = 1 - chi2.cdf(chi2_value, N - 1)

# # Replace the following lines with your actual output code (e.g., GUI updates or print statements)
# print(f"Chi-squared Value: {chi2_value}")
# print(f"P-value: {p_value}")

# # Check if p-value is greater than 0.05
# if p_value > 0.05:
#     print("The distributions are different")
# else:
#     print("The distributions are the same")
