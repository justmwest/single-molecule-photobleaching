#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 07:57:27 2023

@author: justin
"""

import numpy as np
from scipy.stats import norm

# Reference values
n_smalps = 100000
peptides_per_smalp = 2.125

# Reference functions
smalp_core_radius_mean = 38
smalp_core_radius_std_dev = 2
smalp_core_radius_distribution = norm(loc=smalp_core_radius_mean, scale=smalp_core_radius_std_dev)

def pradius(radius):
    """ Takes a radius and returns a probability of finding that radius """
    return smalp_core_radius_distribution.pdf(radius)

# Create the simulation

# First define the number of bins in order to divide up the distribution
# For simplicity I only want whole numbers.
largest_radius_to_consider = smalp_core_radius_mean + smalp_core_radius_std_dev * 5
smallest_radius_to_consider = smalp_core_radius_mean - smalp_core_radius_std_dev * 5
radius_range = largest_radius_to_consider - smallest_radius_to_consider
n_bins = radius_range # This can also be manually defined.

# Now define the list of bin boundaries
bin_edges = np.linspace(smallest_radius_to_consider, largest_radius_to_consider, n_bins+1)
# bin_width = bin_edges[1] - bin_edges[0]

def calc_n_rows(bin_start, bin_end):
    """Using the start and end radii of the bin, and the probability of finding
    a smalp within that size range, calcuclate how many rows to make the matrix."""
    bin_width = bin_end - bin_start
    bin_center = (bin_start + bin_end) / 2
    probability_density = pradius(bin_center) * bin_width
    n_rows = round(n_smalps * probability_density)
    return n_rows

def create_matrix(n_rows, n_columns=3):
    """ 
    For ideal monomers we use 3 columns: empty, one-fluorophore, and two-fluorophores.
    We assume non-specific labeling can occur only once making max fluorophores
    per peptide equal to two. 
    n_columns must be int.
    """
    return np.zeros((n_rows, n_columns))

n_rows_cumulative = 0
for i, _ in enumerate(bin_edges):
    if i > 0:
        n_rows = calc_n_rows(bin_edges[i-1], bin_edges[i])
        n_rows_cumulative += n_rows
        print(f"Bin starting at {bin_edges[i-1]} gets {n_rows} rows ")
print(f"n_smalps: {n_smalps}")
print(f"n_rows_cumulative: {n_rows_cumulative}")
    
    