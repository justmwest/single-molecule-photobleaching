#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 20:43:08 2023

@author: justin

This script shows that calculating the number of smalps of a normal distribution
can be approximated by just dividing the total system area by the mean. It's
not that far off. This likely won't work for a different distribution. 
"""

import numpy as np
from pradius import pradius
from scipy.constants import N_A # Avogadro's number
import matplotlib.pyplot as plt


def concentration_to_number(concentration, volume):
    """ Concentration in M, volume in L """
    return concentration * volume * N_A

def radius_to_area(radius):
    """ Calculates the area of a circle. """
    return (np.pi * radius ** 2)

def gen_radius_bins(n_bins, distribution):
    """
    Takes a distribution and the number of bins (int)
    and returns a list of numbers within 5 standard deviations
    of the mean of the distribution. 
    """
    smallest = distribution.mean() - (distribution.std() * 5)
    largest = distribution.mean() + (distribution.std() * 5)
    radius_list = np.linspace(smallest, largest, n_bins + 1)
    return radius_list
    

def calc_number_of_smalps(total_system_area, distribution, n_bins=25000):
    """ Calculates number of SMALPs from the calculated total lipid area 
    and the radius distribution of SMALPs. Will greatly affect results. """
    
    radius_list = gen_radius_bins(n_bins, distribution)
    
    number_of_smalps = 0
    cumulative_probability = 0
    
    # np.linspace makes the width the same between any two points.
    bin_width = radius_list[1] - radius_list[0] 
    
    for i, _ in enumerate(radius_list):
        if i > 0:
            bin_start_radius = radius_list[i-1]
            bin_end_radius = radius_list[i]
            
            # Calculate the radius in the center of the bin
            bin_average_radius = (bin_start_radius + bin_end_radius) / 2
            
            # Calculate the lipid surface area taken up by a SMALP of a given radius
            # Multiplied by two because two leaflets.
            bin_average_area = (radius_to_area(bin_average_radius) * 2) # Å^2 * 2
    
            # It doesn't matter that we use the probability from the radius
            # distribution since the probability is unitless.
            probability = pradius(bin_average_radius, distribution)
            
            # This is the estimated area under the area distribution curve.
            probability_density = probability * bin_width
            cumulative_probability += probability_density
            
            area_within_bin = total_system_area * probability_density
            
            number_of_smalps_within_bin = area_within_bin / bin_average_area 
            
            number_of_smalps += number_of_smalps_within_bin
    print(f"Cumulative probability: {cumulative_probability}")
    return number_of_smalps
        
    
########################################################################
########################## TESTING FUNCTIONS ###########################
########################################################################

    
def test_calc_number_of_smalps(peptide_concentration, volume_of_solution, peptide_area, total_lipid_area, distribution):
    """ This function shows that there is a asymptotic relationship between
    the number of bins and the number of SMALPs calculated. The choice of 
    number of bins is thus a trade-off between accuracy and computational 
    cost. I'm not sure why small bin numbers bias in one direction though. """

    n_peptides = concentration_to_number(peptide_concentration, volume_of_solution) 
    print(f"n_peptides: {n_peptides}")
    
    # We multiply this number by two because the peptide occupies both leaflets.
    total_peptide_area = (peptide_area*2) * n_peptides
    print(f"Total peptide area: {total_peptide_area} Å^2")
    
    total_system_area = total_lipid_area + total_peptide_area
    print(f"Total system area: {total_system_area} Å^2")
    
    num_bins_list = [1000 * i for i in range(10, 50, 4)]
    
    n_smalps_list =[]
    for num_bins in num_bins_list:
        print(f"\nNumber of bins: {num_bins}")
        
        n_smalps = calc_number_of_smalps(total_system_area, distribution, num_bins)
        print(f"Number of SMALPs: {np.format_float_scientific(n_smalps)}")
        
        n_smalps_list.append(n_smalps)
    
    plt.plot(num_bins_list, n_smalps_list)
    plt.xlabel("Number of bins")
    plt.ylabel("Number of SMALPs calculated")
    plt.show()
    

def calc_number_of_smalps_simple(radius, peptide_concentration, volume_of_solution, peptide_area, total_lipid_area):
    
    n_peptides = concentration_to_number(peptide_concentration, volume_of_solution) 
    print(f"\nn_peptides: {np.format_float_scientific(n_peptides)}")
    
    # We multiply this number by two because the peptide occupies both leaflets.
    total_peptide_area = (peptide_area*2) * n_peptides
    print(f"Total peptide area: {total_peptide_area} Å^2")
    
    total_system_area = total_lipid_area + total_peptide_area
    print(f"Total system area: {total_system_area} Å^2")
    
    number_of_smalps = total_system_area / (radius_to_area(radius)*2)
    print(f"Number of SMALPs: {np.format_float_scientific(number_of_smalps)}")
    
    print(f"Peptides per SMALP: {n_peptides / number_of_smalps}")
    return number_of_smalps


def main():
    from scipy.stats import norm
    
    # User input
    volume_of_solution = 200e-6  # L
    lipid_concentration = 75 * 1e-6  # M
    test_peptide_concentration = 1.5e-6
    test_radius = 38
    
    # Reference values
    lipid_area = 68.3  # Å^2 at 30°C. Kučerka et al. 2005
    smalp_core_radius_mean = 38  # Å Jamshad et al. 2015
    smalp_core_radius_std_dev = 2  # Å
    peptide_radius = 6  # Å
    
    # Calculations
    n_lipids = concentration_to_number(lipid_concentration, volume_of_solution)
    total_lipid_area = lipid_area * n_lipids  # Å^2
    peptide_area = radius_to_area(peptide_radius)  # Å^2
    smalp_core_radius_distribution = norm(loc=smalp_core_radius_mean, scale=smalp_core_radius_std_dev)
    
    # Test zone
    test_calc_number_of_smalps(test_peptide_concentration, volume_of_solution, peptide_area, total_lipid_area, smalp_core_radius_distribution)
    calc_number_of_smalps_simple(test_radius, test_peptide_concentration, volume_of_solution, peptide_area, total_lipid_area)


if __name__ == "__main__":
    main()
