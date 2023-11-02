#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 20:31:58 2023

@author: justin

This is to test the Pradius function.
Do not pull this function for use elsewhere.
The reason is that we cannot modify the distribution in the function.

"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt


def pradius(radius, distribution):
    """ 
    Takes a radius (float or int) and a distribution and returns a probability of finding that radius.
    Requires that the distribution is defined beforehand.
    Tested with distribution made from scipy.stats.norm. Needs to be a callable
    probability distribution function. 
    """
    return distribution.pdf(radius)

def test_pradius():
    """ Tests pradius at a range of values """
    num_bins = 10000
    
    smallest_radius_to_consider = float(smalp_core_radius_mean - (smalp_core_radius_std_dev * 5))
    largest_radius_to_consider = float(smalp_core_radius_mean + (smalp_core_radius_std_dev * 5))
    
    bin_ranges = np.linspace(smallest_radius_to_consider, largest_radius_to_consider, num_bins + 1)
    bin_width = bin_ranges[1] - bin_ranges[0]
    
    probabilities = [pradius(i, smalp_core_radius_distribution) for i in bin_ranges]
    probability_densities = np.array(probabilities) * bin_width
    print(f"Sum of probabilities: {np.sum(probability_densities)}")
    
    plt.plot(bin_ranges, probabilities)
    plt.xlabel("Radius (Å)")
    plt.ylabel("Probability")
    plt.show()
    
if __name__ == "__main__":
    smalp_core_radius_mean = 38
    smalp_core_radius_std_dev = 2
    smalp_core_radius_distribution = norm(loc=smalp_core_radius_mean, scale=smalp_core_radius_std_dev)
    test_pradius()