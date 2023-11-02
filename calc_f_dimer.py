#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 16:17:19 2023

@author: justin


See Chadda & Robertson, 2016. "Measuring membrane protein dimerization equilibrium..."
Section 5.2 (equation 8)

"""
import numpy as np
from scipy.stats import norm
from scipy.constants import N_A
from noninteracting_monomer_simulation import simulate_occupancy, calc_p_star
from noninteracting_dimer_simulation import simulate_occupancy_dimer
from pradius import pradius



def main():
    # Experimental probabilities of 1, 2, or 3 photobleaching step for a single mole fraction.
    # p1expt = 0.3
    # p2expt = 0.5
    # p3expt = 0.2 # probability of 3 or more steps
     
    # # User input
    # n_smalps = 10000
    # peptides_per_smalp = 1
    # pfluorophore = 1.0
    # pnonspecific = 0.2
    # smalp_core_radius_mean = 38
    # smalp_core_radius_std_dev = 2

    # # Calculate number of peptides
    # n_peptides = n_smalps * peptides_per_smalp
    # smalp_core_radius_distribution = norm(loc=smalp_core_radius_mean, scale=smalp_core_radius_std_dev)
    
    # # Define the number of bins in order to divide up the distribution
    # # For simplicity I only want whole numbers.
    # largest_radius_to_consider = smalp_core_radius_mean + smalp_core_radius_std_dev * 5
    # smallest_radius_to_consider = smalp_core_radius_mean - smalp_core_radius_std_dev * 5
    # radius_range = largest_radius_to_consider - smallest_radius_to_consider
    # n_bins = radius_range # This can also be manually defined.

    # # Now define the list of bin boundaries
    # bin_edges = np.linspace(smallest_radius_to_consider, largest_radius_to_consider, n_bins+1)
    
    # # Run the simulations
    # monomer_matrix_list = simulate_occupancy(n_peptides, bin_edges, n_smalps, smalp_core_radius_distribution, pfluorophore, pnonspecific)
    # dimer_matrix_list = simulate_occupancy_dimer(n_peptides, bin_edges, n_smalps, smalp_core_radius_distribution, pfluorophore, pnonspecific)

    # # Calculate the probabilities from the simulation results
    # p1m, p2m, p3m = calc_p_star(monomer_matrix_list, n_smalps)
    # p1d, p2d, p3d = calc_p_star(dimer_matrix_list, n_smalps)
    
    import numpy as np
    from scipy.optimize import minimize
    
    def calculate_R2(FDimer, pl, pe, pm, pd):
        # Calculate R2 for the given FDimer
        R2 = 0
        for n in range(1, 4):
            R2 += (pe[n-1] - ((1 - FDimer) * pm[n-1] + FDimer * pd[n-1]))**2
        return R2
    
    def fit_FDimer(pl, pe, pm, pd):
        # Define the objective function for least-squares minimization
        def objective_func(FDimer):
            return calculate_R2(FDimer, pl, pe, pm, pd)
    
        # Use minimize function to find the value of FDimer that minimizes R2
        result = minimize(objective_func, x0=0.5, bounds=[(0, 1)], method='SLSQP')  # Initial guess: 0.5
    
        return result.x[0]
    

    pl = np.array([1.5, 0.75, 0.46875, 0.375, 0.1875, 0.09375])*1e-6
   
    monomer = np.array([24.45141066, 45.80645161, 61.49425287, 68.6, 72.92225201])
    dimer = np.array([51.41065831, 41.61290323, 31.6091954, 26.4, 23.05630027])
    trimer_plus = np.array([24.13793103, 12.58064516, 6.896551724, 5, 4.021447721])
    
    # Create a 2D array 'pe' by stacking the individual arrays
    pe = np.vstack((monomer, dimer, trimer_plus)).T.tolist()
    
    # Convert to a list of lists
    pe = [list(pe_i) for pe_i in pe]
    
    # Now 'pe' is a list of lists with the desired structure
    print(pe)
    
    pm = [[0.4, 0.5, 0.6], [0.4, 0.5, 0.6], [0.4, 0.5, 0.6], [0.4, 0.5, 0.6], [0.4, 0.5, 0.6]]  # List of pm values for different pl values
    pd = [[0.6, 0.7, 0.8], [0.6, 0.7, 0.8], [0.6, 0.7, 0.8], [0.6, 0.7, 0.8], [0.6, 0.7, 0.8]]  # List of pd values for different pl values
    

    
    # Example usage to fit FDimer as a function of pl
    FDimer_values = []
    for pe_i, pm_i, pd_i in zip(pe, pm, pd):
        FDimer_value = fit_FDimer(pl, pe_i, pm_i, pd_i)
        FDimer_values.append(FDimer_value)
    
    print("FDimer values:", FDimer_values)



    
    
    
    
    pass

if __name__ == "__main__":
    main()

