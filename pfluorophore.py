#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 19:37:35 2023

@author: justin
"""

# Values are baseline subtracted
# Can be found in workbook A488-EphA2-A488.xlsx
# Path = "/Users/Justin/Documents/Data/Conjugations/2020-02-24 TM-EphA2-A488/UV-Vis/A488 and TM-EphA2-A488.xlsx"

# Constants
dilution_factor = 10
a_fluorophore = 0.133015669  # A488 of peptide + fluorophore
epsilon_fluorophore = 72000
cf_fluorophore = 0.182312853  # A280 / A488 of fluorophore alone
a280 = 0.036730266  # A280 of peptide + fluorophore
epsilon_280 = 5500

# Calculations
prot_conc = ((a280 * dilution_factor) - (a_fluorophore * dilution_factor * cf_fluorophore)) / epsilon_280  # eq. 5
p_fluorophore = (a_fluorophore * dilution_factor) / (prot_conc * epsilon_fluorophore)  # eq. 6

# Output
print(f"Protein Concentration: {prot_conc}")
print(f"Protein Subscript[Fluorophore, P]: {p_fluorophore}")
