#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 19:37:35 2023

@author: justin
"""

# Values are baseline subtracted
# Can be found in workbook A488-EphA2-A488.xlsx
# Path = "/Users/Justin/Documents/Data/Conjugations/2020-02-24 TM-EphA2-A488/UV-Vis/A488 and TM-EphA2-A488.xlsx"

def calc_p_fluorophore(a280, a_fluorophore, cf_fluorophore, epsilon_280=5500,
                      epsilon_fluorophore=72000, path_length=1,
                      dilution_factor=10):
    """
    Caclulates pfluorophore from absorbance values.
    Returns pfluorophore. Prints calculated protein concentration.
    All inputs can be int or float.

    Values related to the conjugated fluorophore:
      a_fluorophore:        uncorrected absorbance at emission max.
      epsilon_fluorophore:  molar extinction coefficient of emission max.

    Values related to the intrinsic peptide fluorescence:
      a280:         absorbance at 280 of the conjugated peptide
      epsilon_280:  molar extinction coefficient of peptide alone at 280.

    cf_fluorophore:
      A correction factor for fluorophore absorbance at 280. Calculated as
      A280 / A_at_emission_maximum of the unconjugated fluorophore alone.

    Ref:
        Chadda R, Robertson JL. Measuring Membrane Protein Dimerization
        Equilibrium in Lipid Bilayers by Single-Molecule Fluorescence
        Microscopy. Methods Enzymol. 2016;581:53-82.
        doi: 10.1016/bs.mie.2016.08.025. Epub 2016 Oct 11.
        PMID: 27793292; PMCID: PMC5568537.
    """

    # Calculations
    undiluted_a280 = a280 * dilution_factor
    undiluted_a_fluorophore = a_fluorophore * dilution_factor
    correction = undiluted_a_fluorophore * cf_fluorophore

    prot_conc = (undiluted_a280 - correction) / (epsilon_280 * path_length) # Ref. eq. 5
    p_fluorophore = (undiluted_a_fluorophore / (epsilon_fluorophore * path_length)) / prot_conc  # Ref. eq. 6

    # Output
    print(f"Protein Concentration: {prot_conc}")
    print(f"p_fluorophore: {p_fluorophore}")

    return p_fluorophore

def main():
    # Constants
    path_length = 1
    dilution_factor = 10
    a_fluorophore = 0.133015669  # A_max of peptide + fluorophore
    epsilon_fluorophore = 72000
    cf_fluorophore = 0.182312853  # A280 / A_max of fluorophore alone
    a280 = 0.036730266  # A280 of peptide + fluorophore
    epsilon_280 = 5500

    calc_p_fluorophore(a280, a_fluorophore, cf_fluorophore, epsilon_280,
                      epsilon_fluorophore, path_length, dilution_factor)

if __name__ == "__main__":
    main()