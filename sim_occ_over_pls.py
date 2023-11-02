#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 17:38:52 2023

@author: justin

This is a script to run the noninteracting monomer and dimer simulations for all
values of P/L. 

Ref: Chadda & Robertson 2016: "Measuring Membrane Protein Dimerization Equilibrium in Lipid Bilayers by Single-Molecule Fluorescence Microscopy"

"""

import numpy as np
from scipy.stats import norm
from calc_n_smalps import calc_number_of_smalps, concentration_to_number, radius_to_area
from noninteracting_monomer_simulation import simulate_occupancy, calc_p_star
from noninteracting_dimer_simulation import simulate_occupancy_dimer

def get_pm_pd(peptide_concentration, volume_of_solution, peptide_area,
              total_lipid_area, distribution, bin_edges, pfluorophore, pnonspecific):

    n_peptides = concentration_to_number(peptide_concentration, volume_of_solution) 
    
    total_peptide_area = (peptide_area*2) * n_peptides

    total_system_area = total_lipid_area + total_peptide_area
    
    n_smalps = calc_number_of_smalps(total_system_area, distribution, n_bins=25000)
    
    peptides_per_smalp = n_peptides / n_smalps
    
    # Redefine n_smalps to a smaller number for the simulations
    n_smalps = 10000
    
    # Use the peptide / smalp ratio to calculate the new number of peptides.
    n_peptides = n_smalps * peptides_per_smalp

    # Run the simulations
    monomer_matrix_list = simulate_occupancy(n_peptides, bin_edges, n_smalps, distribution, pfluorophore, pnonspecific)
    dimer_matrix_list = simulate_occupancy_dimer(n_peptides, bin_edges, n_smalps, distribution, pfluorophore, pnonspecific)
    
    pm = calc_p_star(monomer_matrix_list, n_smalps)
    pd = calc_p_star(dimer_matrix_list, n_smalps)
    
    return [pm, pd]

def main():
    # User input
    volume_of_solution = 200e-6  # L
    lipid_concentration = 75 * 1e-6  # M
    pfluorophore = 1.0
    pnonspecific = 0.2
    smalp_core_radius_mean = 38
    smalp_core_radius_std_dev = 2

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
    
    # Define the number of bins in order to divide up the distribution
    # For simplicity I only want whole numbers.
    largest_radius_to_consider = smalp_core_radius_mean + smalp_core_radius_std_dev * 5
    smallest_radius_to_consider = smalp_core_radius_mean - smalp_core_radius_std_dev * 5
    radius_range = largest_radius_to_consider - smallest_radius_to_consider
    n_bins = radius_range # This can also be manually defined.

    # Now define the list of bin boundaries
    # This can probably be integrated somewhere.
    bin_edges = np.linspace(smallest_radius_to_consider, largest_radius_to_consider, n_bins+1)
    
    # Finally import the peptide concentrations
    peptide_concentrations = np.array([1.5, 0.75, 0.46875, 0.375, 0.1875, 0.09375]) * 1e-6
    
    pm_list = []
    pd_list = []
    # Get values
    for peptide_concentration in peptide_concentrations:
        pm, pd = get_pm_pd(peptide_concentration, volume_of_solution, peptide_area, total_lipid_area, smalp_core_radius_distribution, bin_edges, pfluorophore, pnonspecific)
        pm_list.append(pm)
        pd_list.append(pd)
    
    print("")
    for conc, pm, pd in zip(peptide_concentrations, pm_list, pd_list):
        print(f"{conc}, {pm}, {pd}")



if __name__ == "__main__":
    main()