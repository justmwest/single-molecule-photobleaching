#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 20:31:58 2023

@author: justin

"""

import numpy as np
from scipy.stats import norm

smalp_core_radius_mean = 38
smalp_core_radius_std_dev = 2
smalp_core_radius_distribution = norm(loc=smalp_core_radius_mean, scale=smalp_core_radius_std_dev)

def pradius(radius):
    """ Takes a radius and returns a probability of finding that radius """
    return smalp_core_radius_distribution.pdf(radius)

def test_pradius():
    """ Tests pradius at a range of values """
    num_bins = 10000
    
    smallest_radius_to_consider = float(smalp_core_radius_mean - (smalp_core_radius_std_dev * 5))
    largest_radius_to_consider = float(smalp_core_radius_mean + (smalp_core_radius_std_dev * 5))
    
    bin_ranges = np.linspace(smallest_radius_to_consider, largest_radius_to_consider, num_bins + 1)
    bin_width = bin_ranges[1] - bin_ranges[0]
    
    probabilities = [pradius(i) for i in bin_ranges]
    probability_densities = np.array(probabilities) * bin_width
    print(f"Sum of probabilities: {np.sum(probability_densities)}")
    
    plt.plot(bin_ranges, probabilities)
    plt.xlabel("Radius (Å)")
    plt.ylabel("Probability")
    plt.show()
    
if __name__ == "__main__":
    test_pradius()