#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 19:38:22 2023

@author: justin
"""

# See Chadda & Robertson, 2016. "Measuring membrane protein dimerization equilibrium..."
# section 5.2

# Inputs

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.stats as stats

# User input
peptide_concentrations = np.array([1.5, 0.75, 0.46875, 0.375, 0.1875, 0.09375]) * 1e-6
lipid_concentration = 75 * 1e-6  # M

# Assuming normality, this is the radius of the hydrophobic core of the smalp (Jamshad 2015)
core_radius_mean = 38  # Refers to the lipid core of the SMALP.
core_radius_std_dev = 2  

def probability_of_radius(value, core_radius_mean, core_radius_std_dev):
    # Create a normal distribution object
    distribution = stats.norm(loc=core_radius_mean, scale=core_radius_std_dev)
    
    # Calculate the probability for the user-defined value
    probability = distribution.pdf(value)
    
    return probability

# Test probability_of_radius function
user_defined_value = 39  # The value for which you want to calculate the probability
result = probability_of_radius(user_defined_value, core_radius_mean, core_radius_std_dev)
print(f"Probability of finding {user_defined_value} in the normal distribution: {result:.4f}")

def plot_probability_of_radius(core_radius_mean, core_radius_std_dev):
    """ Only used for display. Test function. """
    # Create an array of values from -1 to 1
    x = np.linspace(core_radius_mean-10, core_radius_mean+10, 1000)
    
    # Calculate PRadius for each value of x
    y = [probability_of_radius(r, core_radius_mean, core_radius_std_dev) for r in x]
    
    # Create the plot
    plt.plot(x, y)
    plt.title('PRadius vs. Range')
    plt.xlabel('Range')
    plt.ylabel('PRadius')
    plt.grid(True)
    plt.show()

# Number of SMALPs.
lipids_per_smalp = 132 
volume_of_solution = 200 * 10e-6 # L
moles_lipid = (lipid_concentration*volume_of_solution) # M
avogadros_number = 6.02214076 * 10e23
n_lipids = moles_lipid * avogadros_number
n_smapls = n_lipids / lipids_per_smalp 

def n_smalps_at_radius(radius):
    """ returns the number of smalps given the probability distribution
    radius is int. requires: n_smalps and probability_of_radius (func) """
    return np.round( n_smalps * probability_of_radius(radius) )

def lipid_area(radius):
    """ gives the lipid area of a smalp with this radius. """
    np.pi * np.square(radius)
    


pls = np.flip(peptide_concentrations / lipid_concentration)

# Pull in experimental data

monomer = np.array([24.45141066, 45.80645161, 61.49425287, 68.6, 72.92225201])
dimer = np.array([51.41065831, 41.61290323, 31.6091954, 26.4, 23.05630027])
trimer_plus = np.array([24.13793103, 12.58064516, 6.896551724, 5, 4.021447721])

experimental_occupancy_curve = np.array([monomer, dimer, trimer_plus])
experimental_data = np.array([list(zip(pls[1:], experimental_occupancy_curve[i])) for i in range(3)])

# Fit the data

def dimer_model(pl, keq):
    return ((1 + 4 * pl * keq) - np.sqrt(1 + 8 * pl * keq)) / (4 * pl * keq)

dimer_data = np.array(list(zip(pls[:-1], dimer / 100)))

# Curve fitting
params, _ = curve_fit(dimer_model, dimer_data[:, 0], dimer_data[:, 1], p0=[0.003])

keq = params[0]  # equilibrium constant in mole fraction units (unitless)

R = 1.98720425864083 * 1e-3  # kcal K^-1 mol^-1
T = 298  # K
delta_G = -R * T * np.log(keq)  # kcal mol^-1

xticks = np.arange(0, 0.021, 0.005)
frame_ticks = [xticks, None]
frame_label = ["Peptide to lipid molar ratio", "Fraction Dimer"]
font_size = 18

plt.figure(figsize=(10, 6))
plt.plot(pls[:-1], dimer_data[:, 1], 'o', markersize=10, color='black', label='Data')
plt.plot(pls[:-1], dimer_model(pls[:-1], keq) * 100, color='black', label='Fit')
plt.xlabel(frame_label[0], fontsize=font_size)
plt.ylabel(frame_label[1], fontsize=font_size)
plt.xticks(xticks, [f"{x:.3f}" for x in xticks])
plt.tick_params(axis='both', which='major', labelsize=font_size)
plt.legend(fontsize=font_size)
plt.show()

print(f"ΔG = {delta_G:.3f} kcal mol^-1")
