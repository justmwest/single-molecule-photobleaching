#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 07:57:27 2023

@author: justin
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from tqdm import tqdm

# Reference values
n_smalps = 100000
peptides_per_smalp = 1

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
    a smalp within that size range, calcuclate how many rows to make the matrix.
    Ref 5.1.1, number 2.
    """
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
    Ref 5.1.1, number 2.
    """
    return np.zeros((n_rows, n_columns))

def print_n_rows():
    """ A short function to print the number of rows. Used to optimize n_smalps """
    n_rows_cumulative = 0
    for i, _ in enumerate(bin_edges):
        if i > 0:
            n_rows = calc_n_rows(bin_edges[i-1], bin_edges[i])
            n_rows_cumulative += n_rows
            print(f"Bin starting at {bin_edges[i-1]} gets {n_rows} rows ")
    print(f"n_smalps: {n_smalps}")
    print(f"n_rows_cumulative: {n_rows_cumulative}")
        
    
# Calculate the number of peptides and the degree of labeling.
# Ref section 5.1.1, number 3. 
n_peptides = n_smalps * peptides_per_smalp
pfluorophore = 1.0 # Probability of labeling
pnonspecific = 0 # Probability of nonspecific labeling


def calc_n_labeled(n_peptides, pfluorophore=1, pnonspecific=0):
    """ calc degree of labeling. Returns list of numbers of peptides. """
    pspecific = pfluorophore - pnonspecific 
    # Calculate the probabilities of different labeling scenarios for a single subunit
    # Calculate the number of monomers with 0, 1, and 2 labels
    n_with_0_labels = round((1 - pspecific) * (1 - pnonspecific) * n_peptides)
    n_with_1_label = round(pspecific * (1 - pnonspecific) * n_peptides)
    n_with_2_labels = round(pspecific * pnonspecific * n_peptides)
    
    return n_with_0_labels, n_with_1_label, n_with_2_labels


def calc_fraction_unoccupied():
    """ 
    f0 is the fraction of unoccupied smalps.
    As protein concentration increases, f0 decreases.
    However, very small smalps may be inaccesible to dimers.
    Authors found that 40% of liposomes were empty for ClC-ec1 dimers.
    Ref section 4.2.2
    For now, we assume no SMALPs are off-limits, since the smallest size we
    consider, radius 28, has an area of 2463 Å^2 and a peptide dimer has a 
    radius of 226 Å^2.
    Update with empirical data if possible.
    """
    return 0

def weighted_probability(radius, bin_width):
    """ probability density weighted by the square of the radius """
    return radius**2 * pradius(radius) * bin_width

def cumulative_weighted_probability(bin_edges):
    """ calculate the cumulative weighted probability.
    Ref section 5.1.1 number 4, denominator. """
    bin_width = bin_edges[1] - bin_edges[0]
    cumulative_weighted_probability = 0
    for i, _ in enumerate(bin_edges):
        if i > 0:
            bin_start = bin_edges[i-1]
            bin_end = bin_edges[i]
            bin_center = (bin_start + bin_end) / 2
            prob = weighted_probability(bin_center, bin_width)
            cumulative_weighted_probability += prob
    return cumulative_weighted_probability

def n_monomers_in_bin(n_monomers, bin_start, bin_end, bin_edges):
    """
    "The proportion of each monomer species to be inserted into the liposomes
    is calculated using the fractional surface area from the accessible
    liposome population, Pradius,M(r) determined from the measurements
    of F0 (Section 4.2.2)."
    Ref section 5.1.1 number 4
    We assume Pradius,M(r) = Pradius(r), since we assume F0 to be 0. 
    """
    bin_width = bin_end - bin_start
    bin_center = (bin_start + bin_end) / 2
    numerator = weighted_probability(bin_center, bin_width)
    denominator = cumulative_weighted_probability(bin_edges)
    n_monomers_in_bin = round((numerator / denominator) * n_monomers)
    return n_monomers_in_bin

def randomly_assign_occupancy(n_labeled_list, matrix):
    """ 
    n_labeled is a list of 3 numbers of peptides. Position 0 has no labels, 
    Position 1 has 1 label. Position 2 has 2 labels.
    Modifies the matrix and returns it.
    """
    # For the monomer simulation it's convenient that the j from enumerate
    # is the same as both the column number and the number of labels.
    # This might be confusing. Should refactor for dimer simulation.
    if len(matrix) > 0:
        for j, n_labeled in enumerate(n_labeled_list):
            if n_labeled > 0:
                # Choose random positions
                random_positions = np.random.choice(len(matrix), n_labeled, replace=True)
                
                # Modify random positions
                # Could do it in one line like this:
                #   matrix[:,j][random_positions] += j
                # However, that way doesn't modify the same smalp twice.
                for position in random_positions:
                    matrix[:,j][position] += j
    
    return matrix
    

def matrix_to_histogram(matrix):
    num_columns = matrix.shape[1]

    for col in range(num_columns):
        if matrix[:, col].max() > 0:
            n_bins = int(matrix[:, col].max())
        else:
            n_bins = 3
        plt.figure()
        plt.hist(matrix[:, col], bins=n_bins, color='skyblue', edgecolor='black', alpha=0.9)
        plt.ylabel('Frequency')
        plt.xlabel('Number of fluorophores')
        plt.title(f'Distribution of peptides with {col} label(s)')

    plt.show()



def simulation(n_peptides, bin_edges):
    """ 
    asdf 
    Section 5.1.1 number 5
    """
    # bin_width = bin_edges[1] - bin_edges[0]
    matrix_list = []
    
    # For each radius bin
    for i in tqdm(range(len(bin_edges))):
        if i > 0:
            bin_start = bin_edges[i-1]
            bin_end = bin_edges[i]
            bin_center = (bin_start + bin_end) / 2 # Average radius in bin
            
            # Calculate how many SMALPs have that radius
            n_rows = calc_n_rows(bin_start, bin_end)
            
            # Create the matrix
            matrix = create_matrix(n_rows)
            
            # Calculate number of monomers in the bin
            n_peptides_in_bin = n_monomers_in_bin(n_peptides, bin_start, bin_end, bin_edges)
            
            # Calculate fraction of each labeled species
            n_labeled_list = calc_n_labeled(n_peptides_in_bin)
            
            # Randomly assign occupancy
            occupied_matrix = randomly_assign_occupancy(n_labeled_list, matrix)
            
            matrix_list.append(occupied_matrix)

            if i < 5:
                print(f"\nRadius: {bin_center}")
                print(f"N peptides in bin: {n_peptides_in_bin}")
                print(f"Number of rows in the matrix: {n_rows}")
                print(f"N_labeled_list: {n_labeled_list}")
                print(occupied_matrix)
    
    return matrix_list
                
                
            
            
            
    
matrix_list = simulation(n_peptides, bin_edges)

matrix_to_histogram(matrix_list[10])
    



            
    
    