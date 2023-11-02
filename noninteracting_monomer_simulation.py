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
from pradius import pradius

def calc_n_rows(bin_start, bin_end, n_smalps, distribution):
    """Using the start and end radii of the bin, and the probability of finding
    a smalp within that size range, calcuclate how many rows to make the matrix.
    Ref 5.1.1, number 2.
    """
    bin_width = bin_end - bin_start
    bin_center = (bin_start + bin_end) / 2
    probability_density = pradius(bin_center, distribution) * bin_width
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


def calc_n_labeled(n_peptides, pfluorophore=1, pnonspecific=0):
    """ Calculate how many peptides get 0, 1 or 2 labels.
    Returns list of numbers of peptides in that order.
    pfluorophore and pnonspecific are probabilities, so are floats between 0 and 1"""
    pspecific = pfluorophore - pnonspecific 
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

def weighted_probability(radius, bin_width, distribution):
    """
    Probability density weighted by the square of the radius.
    Ref section 5.1.1 number 4 (Equation 7)
    """
    return radius**2 * pradius(radius, distribution) * bin_width

def cumulative_weighted_probability(bin_edges, distribution):
    """ calculate the cumulative weighted probability.
    Ref section 5.1.1 number 4, denominator. """
    bin_width = bin_edges[1] - bin_edges[0]
    cumulative_weighted_probability = 0
    for i, _ in enumerate(bin_edges):
        if i > 0:
            bin_start = bin_edges[i-1]
            bin_end = bin_edges[i]
            bin_center = (bin_start + bin_end) / 2
            prob = weighted_probability(bin_center, bin_width, distribution)
            cumulative_weighted_probability += prob
    return cumulative_weighted_probability

def n_monomers_in_bin(n_monomers, bin_start, bin_end, bin_edges, distribution):
    """"
    Calculates the number of peptides to be assigned in a single radius bin.
    
    We assume Pradius,M(r) = Pradius(r), since we assume F0 to be 0. 
    
    "The proportion of each monomer species to be inserted into the liposomes
    is calculated using the fractional surface area from the accessible
    liposome population, Pradius,M(r) determined from the measurements
    of F0 (Section 4.2.2)."
    Ref section 5.1.1 number 4 (Equation 7)
    """
    bin_width = bin_end - bin_start
    bin_center = (bin_start + bin_end) / 2
    numerator = weighted_probability(bin_center, bin_width, distribution)
    denominator = cumulative_weighted_probability(bin_edges, distribution)
    n_monomers_in_bin = round((numerator / denominator) * n_monomers)
    return n_monomers_in_bin

def randomly_assign_occupancy(n_labeled_list, matrix):
    """ 
    n_labeled_list is a list of numbers of peptides. Position 0 is the number
    of peptides with no labels. Position 1 is the number with one label and so on. 
    Modifies the input matrix and returns it.
    """
    # For the monomer simulation it's convenient that the i from enumerate
    # is the same as both the column number and the number of labels to add.
    # This might be confusing. Should refactor for dimer simulation.
    if len(matrix) > 0:
        for i, n_labeled in enumerate(n_labeled_list):
            if n_labeled > 0:
                # Choose random positions
                random_positions = np.random.choice(len(matrix), n_labeled, replace=True)
                
                # Modify random positions
                # Could do it in one line like this:
                #   matrix[:,j][random_positions] += j
                # However, that way doesn't modify the same smalp twice.
                for position in random_positions:
                    matrix[:, i][position] += i
    
    return matrix

def simulate_occupancy(n_peptides, bin_edges, n_smalps, distribution, pfluorophore, pnonspecific):
    """ 
    Simulate occupancy of monomers. Return list of matrices, one for each
    radius bin considered.
    Section 5.1.1 number 5
    """
    matrix_list = []
    
    # For each radius bin
    for i in tqdm(range(len(bin_edges)), desc="Running siumulation for each radius"):
        if i > 0:
            bin_start = bin_edges[i-1]
            bin_end = bin_edges[i]
            
            # Calculate how many SMALPs have that radius
            n_rows = calc_n_rows(bin_start, bin_end, n_smalps, distribution)
            
            # Create the matrix
            matrix = create_matrix(n_rows)
            
            # Calculate number of monomers in the bin
            n_peptides_in_bin = n_monomers_in_bin(n_peptides, bin_start, bin_end, bin_edges, distribution)
            
            # Calculate fraction of each labeled species
            n_labeled_list = calc_n_labeled(n_peptides_in_bin, pfluorophore, pnonspecific)
            
            # Randomly assign occupancy
            occupied_matrix = randomly_assign_occupancy(n_labeled_list, matrix)
            
            matrix_list.append(occupied_matrix)
    
    return matrix_list


def calc_p_star(matrix_list, n_smalps):
    """ Calculate the probability of observing 1, 2, or 3+ steps. 
    Ref: section 5.1.1 numbers 6 and 7 """
    # Stack the matrix into one tall array of 3 columns.
    combined_matrix = np.vstack(matrix_list)
    
    # Count up the results
    uniques, counts = np.unique(combined_matrix, return_counts=True)
    
    frequencies = counts / n_smalps
    
    # This is normally not zero, but as explained in the function, we assume.
    # It is not taken from the simulation results.
    p0 = calc_fraction_unoccupied()

    p1 = frequencies[1] / (1 - p0)
    
    p2 = frequencies[2] / (1 - p0)
    
    p3plus = sum(frequencies[3:]) / (1 - p0)
    
    return p1, p2, p3plus

########################################################################
########################## TESTING FUNCTIONS ###########################
########################################################################

def print_n_rows(bin_edges, n_smalps):
    """ Print the number of rows. Used to optimize n_smalps """
    n_rows_cumulative = 0
    for i, _ in enumerate(bin_edges):
        if i > 0:
            n_rows = calc_n_rows(bin_edges[i-1], bin_edges[i])
            n_rows_cumulative += n_rows
            print(f"Bin starting at {bin_edges[i-1]} gets {n_rows} rows ")
    print(f"n_smalps: {n_smalps}")
    print(f"n_rows_cumulative: {n_rows_cumulative}")


def print_sim_results(matrix_list, n_peptides):
    """ Prints counts of different fluorophore numbers """
    
    # Stack the matrix into one tall array of 3 columns.
    combined_matrix = np.vstack(matrix_list)
    
    # Count up the results
    uniques, counts = np.unique(combined_matrix, return_counts=True)
    
    total_spots = 0
    print("")
    for i, (unique, count) in enumerate(zip(uniques, counts)):
        frequency = count / n_peptides
        print(f"{unique}: {count} / {n_peptides} = {frequency}")
        
        # Add up the number of spots observed, skipping over 0-steps
        if i > 0:
            total_spots += count
            
    print(f"Sum of the entire matrix (number of fluorophores): {combined_matrix.sum()}")
    print(f"Total number of peptides: {n_peptides}")
    print(f"Total spots observed: {total_spots}")
    
def matrix_to_histogram(matrix):
    """ Show the results of the simulation in histogram format
    Ref section 5.1.1 number 6, first sentence.
    """
    num_columns = matrix.shape[1]

    for col in range(num_columns):
        if matrix[:, col].max() > 0:
            n_bins = int(matrix[:, col].max())
        else:
            n_bins = 3
        plt.figure()
        plt.hist(matrix[:, col], bins=n_bins, color='skyblue', edgecolor='black', alpha=0.9)
        plt.ylabel('Frequency')
        plt.xlabel('Theoretical number of fluorophores observed')
        plt.title(f'Distribution of peptides with {col} label(s)')

    plt.show()


def main():
    """
    This function is for testing the simulation.
    Uncomment tests that you want to run.
    """
    
    # Reference values
    n_smalps = 100000
    peptides_per_smalp = 1.15
    smalp_core_radius_mean = 38
    smalp_core_radius_std_dev = 2
    smalp_core_radius_distribution = norm(loc=smalp_core_radius_mean, scale=smalp_core_radius_std_dev)
    
    # First define the number of bins in order to divide up the distribution
    # For simplicity I only want whole numbers.
    largest_radius_to_consider = smalp_core_radius_mean + smalp_core_radius_std_dev * 5
    smallest_radius_to_consider = smalp_core_radius_mean - smalp_core_radius_std_dev * 5
    radius_range = largest_radius_to_consider - smallest_radius_to_consider
    n_bins = radius_range # This can also be manually defined.

    # Now define the list of bin boundaries
    bin_edges = np.linspace(smallest_radius_to_consider, largest_radius_to_consider, n_bins+1)

    # Calculate the number of peptides and the degree of labeling.
    # Ref section 5.1.1, number 3. 
    n_peptides = n_smalps * peptides_per_smalp
    pfluorophore = 1.0 # Probability of labeling
    pnonspecific = 0.0 # Probability of nonspecific labeling
    
    # Run the simulation
    matrix_list = simulate_occupancy(n_peptides, bin_edges, n_smalps, smalp_core_radius_distribution, pfluorophore, pnonspecific)
    
    # Below are various tests that you can uncomment to show.
    
    # # Print the number of rows
    # print_n_rows(bin_edges, n_smalps)
    
    # # Print the results to the python console.
    # print_sim_results(matrix_list, n_peptides)
    
    # # Show the results as a histogram
    # combined_matrix = np.vstack(matrix_list)
    # matrix_to_histogram(combined_matrix)
    
    # # Show probabilities
    p1, p2, p3 = calc_p_star(matrix_list, n_smalps)
    print(f"\nP(1): {p1}")
    print(f"P(2): {p2}")
    print(f"P(3+): {p3}")
    
if __name__ == "__main__":
    main()


            
    
    