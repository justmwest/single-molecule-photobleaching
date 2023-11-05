#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 13:08:33 2023

@author: justin

This is a script to perform the "full analysis." First, use the experimental
parameters to calculate the nubmer of SMALPs. Then simulate the noninteracting
monomer and dimer populations. Then use those simulations to calculate Fdimer.
Finally, fit Fdimer to find Keq and thus deltaG.

Pay attention to the inputs.

Translated to python and modified for SMALPs from:
    Chadda R, Robertson JL. Measuring Membrane Protein Dimerization
    Equilibrium in Lipid Bilayers by Single-Molecule Fluorescence
    Microscopy. Methods Enzymol. 2016;581:53-82.
    doi: 10.1016/bs.mie.2016.08.025. Epub 2016 Oct 11.
    PMID: 27793292; PMCID: PMC5568537.

"""

# Import standard libraries. If these fail, you need to install them. e.g.:
#   conda install numpy
import numpy as np
from scipy.constants import N_A
from scipy.stats import norm

# Import from other scripts that should be present in the same directory.
from calc_n_smalps import concentration_to_number, radius_to_area, calc_number_of_smalps
from noninteracting_monomer_simulation import simulate_occupancy, calc_p_star
from noninteracting_dimer_simulation import simulate_occupancy_dimer
from calc_f_dimer import fit_frac_dimer
from plot_probs import plot_probability_curves, plot_probability_comparison
from calc_delta_g import find_equilibrium_constant, calc_delta_g, plot_keq_fit


##############################################################################
###                             USER INPUTS                                ###
##############################################################################

# Information from experiments
volume_of_solution = 200e-6  # L
lipid_concentration = 75e-6  # M

pfluorophore = 1 # Probability of labeling
pnonspecific = 0 # Probability of nonspecific labeling

# These are concentrations in molar units. Must be same length as data.
peptide_concentrations = np.array([1.5, 0.75, 0.46875, 0.375, 0.1875]) * 1e-6

# Experimentally determined probabilities of finding 1, 2, or 3+ photobleaching
# steps. Must be in the same order as peptide concentrations.
p_exp_list = [[0.024451410660000002, 0.05141065831, 0.02413793103],
              [0.04580645161000001, 0.041612903230000005, 0.012580645159999999],
              [0.061494252869999995, 0.0316091954, 0.0068965517240000005],
              [0.0686, 0.0264, 0.005],
              [0.07292225201, 0.023056300270000003, 0.004021447721]]

# Reference values
lipid_area = 68.3  # Å^2 for POPC at 30°C. Kučerka et al. 2005
smalp_core_radius_mean = 38  # Å Jamshad et al. 2015
smalp_core_radius_std_dev = 2  # Å
peptide_radius = 6  # Å, if helix diameter is 12 Å.

# How many SMALPs to simulate. Smaller is faster. Larger, more accurate.
n_simulated_smalps = 1000000
n_simulated_radius_bins = 20

##############################################################################
###                             MAIN FUNCTION                              ###
##############################################################################

# Constants
avogadros_number = N_A # See imports

# Calculations that apply to each iteration below
n_lipids = concentration_to_number(lipid_concentration, volume_of_solution)
total_lipid_area = lipid_area * n_lipids  # Å^2
peptide_area = radius_to_area(peptide_radius)  # Å^2
smalp_core_radius_distribution = norm(loc=smalp_core_radius_mean,
                                      scale=smalp_core_radius_std_dev)

# plot experimental data

p_monomer_list = []
p_dimer_list = []
f_dimer_list = []
for peptide_concentration, pe in zip(peptide_concentrations, p_exp_list):

    print(f"\nStarting iteration for P/L: {peptide_concentration/lipid_concentration:.3f}")

    # Calculate total surface area of all of the lipids and peptides
    n_peptides = concentration_to_number(peptide_concentration,
                                         volume_of_solution)
    total_peptide_area = (peptide_area*2) * n_peptides
    total_system_area = total_lipid_area + total_peptide_area

    # Calculate the number of SMALPs from lipid conc. and the area per lipid.
    n_smalps = calc_number_of_smalps(total_system_area,
                                     smalp_core_radius_distribution,
                                     n_bins=25000)

    print(f"Calculated that {int(n_peptides):.2e} peptides and {int(n_lipids):.2e} "
          f"lipids fit into {int(n_smalps):.2e} SMALPs of radius "
          f"{smalp_core_radius_mean}.")

    peptides_per_smalp = n_peptides / n_smalps
    n_peptides_simulated = round(n_simulated_smalps * peptides_per_smalp)

    print(f"With {peptides_per_smalp:.2f} peptides per SMALP, simulating "
          f"occupancy of {n_peptides_simulated:.2e} peptides in "
          f"{n_simulated_smalps:.2e} SMALPs.")

    # Simulate occupancy of SMALPs with noninteracting monomers.
    monomer_matrix_list = simulate_occupancy(n_peptides_simulated,
                                             n_simulated_smalps,
                                             smalp_core_radius_distribution,
                                             n_simulated_radius_bins,
                                             pfluorophore, pnonspecific)
    # Add up all the occupancies from the simulation
    pm = calc_p_star(monomer_matrix_list, n_simulated_smalps)
    p_monomer_list.append(pm)

    # Simulate occupancy of SMALPs with noninteracting dimers.
    dimer_matrix_list = simulate_occupancy_dimer(n_peptides_simulated,
                                                 n_simulated_smalps,
                                                 smalp_core_radius_distribution,
                                                 n_simulated_radius_bins,
                                                 pfluorophore, pnonspecific)
    # Add up all the occupancies from the simulation
    pd = calc_p_star(dimer_matrix_list, n_simulated_smalps)
    p_dimer_list.append(pd)

    # Calc f dimer using simulations and experimental data.
    fd_pred = fit_frac_dimer(pe, pm, pd)
    f_dimer_list.append(fd_pred)

# Calculate delta G
k_eq, fit_params = find_equilibrium_constant(pl_list, f_dimer_list)
delta_g = calc_delta_g(k_eq)

print(f"\nCalculated delta G: {delta_g}")
print(
      """Note: if any of the distributions were statistically different,
      then that means the fitting used to determine FDimer did not go well.
      This means don't trust the delta g value."""
      )

##############################################################################
###                                 PLOTS                                  ###
##############################################################################


# Plot simulation results
pl_list = peptide_concentrations / lipid_concentration
p1e, p2e, p3e = zip(*p_exp_list)
plot_probability_curves(pl_list, p1e, p2e, p3e, title="Experimental data")

p1m, p2m, p3m = zip(*p_monomer_list)
plot_probability_curves(pl_list, p1m, p2m, p3m, title="Monomer simulation")

p1d, p2d, p3d = zip(*p_dimer_list)
plot_probability_curves(pl_list, p1d, p2d, p3d, title="Dimer simulation")

# Plot comparisons
plot_probability_comparison(pl_list, p1e, p1m, p1d, title="$p_1$")
plot_probability_comparison(pl_list, p2e, p2m, p2d, title="$p_2$")
plot_probability_comparison(pl_list, p3e, p3m, p3d, title="$p_{3+}$")

# Plot K_eq fit
plot_keq_fit(pl_list, f_dimer_list, fit_params, delta_g)

