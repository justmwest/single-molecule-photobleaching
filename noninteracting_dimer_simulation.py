#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 14:21:26 2023

@author: justin

See Chadda & Robertson, 2016. "Measuring membrane protein dimerization equilibrium..."
Section 5.1.2

Reusing a lot of code from the noninteracting monomer simulation.

"""

import numpy as np
from tqdm import tqdm
from scipy.stats import norm
from noninteracting_monomer_simulation import calc_n_rows, create_matrix, n_monomers_in_bin, randomly_assign_occupancy, calc_p_star


def calc_n_labeled_dimer(pfluorophore, pnonspecific, n_peptides):
    """
    Calculate how many dimers get 0, 1, 2, 3, or 4 labels.
    Returns list of numbers of peptides in that order.
    pfluorophore and pnonspecific are probabilities, so are floats between 0 and 1
    """
    
    pspecific = pfluorophore - pnonspecific
    n_dimers = n_peptides / 2
    
    # Calculate monomer label probabilities
    P0 = (1 - pspecific) * (1 - pnonspecific)
    P1 = pspecific * (1 - pnonspecific) + (1 - pspecific) * pnonspecific
    P2 = pspecific * pnonspecific

    # Calculate dimer label probabilities
    P0_dimer = P0 * P0
    P1_dimer = (P0 * P1) + (P1 * P0)
    P2_dimer = (P0 * P2) + (P1 * P1) + (P2 * P0)
    P3_dimer = (P1 * P2) + (P2 * P1)
    P4_dimer = P2 * P2

    # Calculate the expected number of dimers with 0, 1, 2, 3, and 4 labels
    num_dimers_0_labels = round(n_dimers * P0_dimer)
    num_dimers_1_label = round(n_dimers * P1_dimer)
    num_dimers_2_labels = round(n_dimers * P2_dimer)
    num_dimers_3_labels = round(n_dimers * P3_dimer)
    num_dimers_4_labels = round(n_dimers * P4_dimer)

    return num_dimers_0_labels, num_dimers_1_label, num_dimers_2_labels, num_dimers_3_labels, num_dimers_4_labels


def simulate_occupancy_dimer(n_peptides, bin_edges, n_smalps, distribution, 
                             pfluorophore=1.0, pnonspecific=0.0, verbose=False):
    """ 
    Simulate occupancy of monomers. Return list of matrices, one for each
    radius bin considered.
    Section 5.1.1 number 5
    """
    matrix_list = []
    
    # For each radius bin
    for i in tqdm(range(len(bin_edges)), desc="Running dimer siumulation for each radius"):
        if i > 0:
            bin_start = bin_edges[i-1]
            bin_end = bin_edges[i]
            
            # Calculate how many SMALPs have that radius
            n_rows = calc_n_rows(bin_start, bin_end, n_smalps, distribution)
            
            # Create the matrix
            matrix = create_matrix(n_rows, 5)        
            
            # Calculate number of monomers in the bin
            n_peptides_in_bin = n_monomers_in_bin(n_peptides, bin_start, bin_end, bin_edges, distribution)
            if verbose: print(f"Bin ending at: {bin_end}:\t n_rows: {n_rows}.\t n_peptides: {n_peptides_in_bin}")
            
            # Calculate fraction of each labeled species
            n_labeled_list = calc_n_labeled_dimer(n_peptides_in_bin, pfluorophore, pnonspecific)
            
            # Randomly assign occupancy
            occupied_matrix = randomly_assign_occupancy(n_labeled_list, matrix)
            
            matrix_list.append(occupied_matrix)
    
    return matrix_list



def main():
    # User input
    n_smalps = 10000
    peptides_per_smalp = 1
    pfluorophore = 1.0
    pnonspecific = 0.2
    smalp_core_radius_mean = 38
    smalp_core_radius_std_dev = 2

    # Calculate number of peptides
    n_peptides = n_smalps * peptides_per_smalp
    smalp_core_radius_distribution = norm(loc=smalp_core_radius_mean, scale=smalp_core_radius_std_dev)
    
    # Define the number of bins in order to divide up the distribution
    # For simplicity I only want whole numbers.
    largest_radius_to_consider = smalp_core_radius_mean + smalp_core_radius_std_dev * 5
    smallest_radius_to_consider = smalp_core_radius_mean - smalp_core_radius_std_dev * 5
    radius_range = largest_radius_to_consider - smallest_radius_to_consider
    n_bins = radius_range # This can also be manually defined.

    # Now define the list of bin boundaries
    bin_edges = np.linspace(smallest_radius_to_consider, largest_radius_to_consider, n_bins+1)
    
    matrix_list = simulate_occupancy_dimer(n_peptides, bin_edges, n_smalps, smalp_core_radius_distribution, pfluorophore, pnonspecific)
    

    
    
    # # Show probabilities
    p1, p2, p3 = calc_p_star(matrix_list, n_smalps)
    print(f"\nP(1): {p1}")
    print(f"P(2): {p2}")
    print(f"P(3+): {p3}")
    

if __name__ == "__main__":
    main()