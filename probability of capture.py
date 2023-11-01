#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 19:38:22 2023

@author: justin
"""

# See Chadda & Robertson, 2016. "Measuring membrane protein dimerization equilibrium..."
# section 5.2

# Inputs
def old():
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    import scipy.stats as stats
    
    # User input
    peptide_concentrations = np.array([1.5, 0.75, 0.46875, 0.375, 0.1875, 0.09375]) * 1e-6
    lipid_concentration = 75 * 1e-6  # M
    
    # Assuming normality, this is the radius of the hydrophobic core of the smalp (Jamshad 2015)
    core_radius_mean = 38  # Refers to the lipid core of the SMALP.
    core_radius_std_dev = 2  
    
    def probability_of_radius(value, core_radius_mean, core_radius_std_dev):
        # Create a normal distribution object
        distribution = stats.norm(loc=core_radius_mean, scale=core_radius_std_dev)
        
        # Calculate the probability for the user-defined value
        probability = distribution.pdf(value)
        
        return probability
    
    # Test probability_of_radius function
    user_defined_value = 39  # The value for which you want to calculate the probability
    result = probability_of_radius(user_defined_value, core_radius_mean, core_radius_std_dev)
    print(f"Probability of finding {user_defined_value} in the normal distribution: {result:.4f}")
    
    def plot_probability_of_radius(core_radius_mean, core_radius_std_dev):
        """ Only used for display. Test function. """
        # Create an array of values from -1 to 1
        x = np.linspace(core_radius_mean-10, core_radius_mean+10, 1000)
        
        # Calculate PRadius for each value of x
        y = [probability_of_radius(r, core_radius_mean, core_radius_std_dev) for r in x]
        
        # Create the plot
        plt.plot(x, y)
        plt.title('PRadius vs. Range')
        plt.xlabel('Range')
        plt.ylabel('PRadius')
        plt.grid(True)
        plt.show()
    
    # Number of SMALPs.
    lipids_per_smalp = 132 
    volume_of_solution = 200 * 10e-6 # L
    moles_lipid = (lipid_concentration*volume_of_solution) # M
    avogadros_number = 6.02214076e23
    n_lipids = moles_lipid * avogadros_number
    n_smalps = n_lipids / lipids_per_smalp 
    
    def n_smalps_at_radius(radius, n_smalps):
        """ returns the number of smalps given the probability distribution
        radius is int. requires: n_smalps and probability_of_radius (func) """
        return np.round( n_smalps * probability_of_radius(radius) )
    
    def lipid_area(radius):
        """ gives the lipid area of a smalp with this radius. """
        np.pi * np.square(radius)
    
    # Define the experimental mole fraction.
    # See Chadda & Robertson 2016: "Measuring Membrane Protein Dimerization Equilibrium in Lipid Bilayers by Single-Molecule Fluorescence Microscopy"     
    # Section 5.1.1, number 1
    def calc_experimental_mole_fraction():
        """ Calculates experimental mole fraction """
        n_subunits = 1
        n_smalps = 1
        surface_area_per_lipid = 68.3 # Å^2 for POPC at 30ºC, Kučerka et al. 2005
        peptide_diameter = 12 # Å
        peptide_area = np.pi * np.square(peptide_diameter/2) # Å^2. For one monolayer.   
    
    
    
    pls = np.flip(peptide_concentrations / lipid_concentration)
    
    # Pull in experimental data
    
    monomer = np.array([24.45141066, 45.80645161, 61.49425287, 68.6, 72.92225201])
    dimer = np.array([51.41065831, 41.61290323, 31.6091954, 26.4, 23.05630027])
    trimer_plus = np.array([24.13793103, 12.58064516, 6.896551724, 5, 4.021447721])
    
    experimental_occupancy_curve = np.array([monomer, dimer, trimer_plus])
    experimental_data = np.array([list(zip(pls[1:], experimental_occupancy_curve[i])) for i in range(3)])
    
    # Fit the data
    
    def dimer_model(pl, keq):
        return ((1 + 4 * pl * keq) - np.sqrt(1 + 8 * pl * keq)) / (4 * pl * keq)
    
    dimer_data = np.array(list(zip(pls[:-1], dimer / 100)))
    
    # Curve fitting
    params, _ = curve_fit(dimer_model, dimer_data[:, 0], dimer_data[:, 1], p0=[0.003])
    
    keq = params[0]  # equilibrium constant in mole fraction units (unitless)
    
    R = 1.98720425864083 * 1e-3  # kcal K^-1 mol^-1
    T = 298  # K
    delta_G = -R * T * np.log(keq)  # kcal mol^-1
    
    xticks = np.arange(0, 0.021, 0.005)
    frame_ticks = [xticks, None]
    frame_label = ["Peptide to lipid molar ratio", "Fraction Dimer"]
    font_size = 18
    
    plt.figure(figsize=(10, 6))
    plt.plot(pls[:-1], dimer_data[:, 1], 'o', markersize=10, color='black', label='Data')
    plt.plot(pls[:-1], dimer_model(pls[:-1], keq) * 100, color='black', label='Fit')
    plt.xlabel(frame_label[0], fontsize=font_size)
    plt.ylabel(frame_label[1], fontsize=font_size)
    plt.xticks(xticks, [f"{x:.3f}" for x in xticks])
    plt.tick_params(axis='both', which='major', labelsize=font_size)
    plt.legend(fontsize=font_size)
    plt.show()
    
    print(f"ΔG = {delta_G:.3f} kcal mol^-1")


def main():
    import numpy as np
    from scipy.stats import norm
    from scipy.constants import N_A # Avogadro's number
    import matplotlib.pyplot as plt
    from tqdm import tqdm
    
    # Ref: Chadda & Robertson 2016: "Measuring Membrane Protein Dimerization Equilibrium in Lipid Bilayers by Single-Molecule Fluorescence Microscopy"
    
    # Reference values
    lipid_area = 68.3  # Å^2 at 30°C. Kučerka et al. 2005
    smalp_core_radius_mean = 38  # Å Jamshad et al. 2015
    smalp_core_radius_std_dev = 2  # Å
    peptide_radius = 6  # Å
    
    # Constants
    avogadros_number = N_A 
    
    # Inputs
    volume_of_solution = 200e-6  # L
    peptide_concentrations = np.array([1.5, 0.75, 0.46875, 0.375, 0.1875, 0.09375]) * 1e-6
    lipid_concentration = 75 * 1e-6  # M
    
    # Calculation functions
    def concentration_to_number(concentration, volume):
        """ Concentration in M, volume in L """
        return concentration * volume * avogadros_number
    
    def radius_to_area(radius):
        """ Calculates the area of a circle. """
        return (np.pi * radius ** 2)

    
    # Calculations
    pl_ratios = peptide_concentrations / lipid_concentration
    n_lipids = concentration_to_number(lipid_concentration, volume_of_solution)
    total_lipid_area = lipid_area * n_lipids  # Å^2
    peptide_area = radius_to_area(peptide_radius)  # Å^2
    
    # Create a normal distribution for SMALP core radii
    smalp_core_radius_distribution = norm(loc=smalp_core_radius_mean, scale=smalp_core_radius_std_dev)
    
    def pradius(radius):
        """ Takes a radius and returns a probability of finding that radius """
        return smalp_core_radius_distribution.cdf(radius)
    
    def calc_number_of_smalps(radius_list, total_system_area):
        """ Calculates number of SMALPs from the calculated total lipid area 
        and the radius distribution of SMALPs. Will greatly affect results. """
        total_bilayer_area = total_system_area / 2
        
        number_of_smalps = 0
        
        # np.linspace makes the width the same between any two points.
        bin_width = radius_list[1] - radius_list[0] 
        


        for i, _ in enumerate(radius_list):
            if i > 0:
                bin_start_radius = radius_list[i-1]
                bin_end_radius = radius_list[i]
                
                # Calculate the radius in the center of the bin
                bin_average_radius = (bin_start_radius + bin_end_radius) / 2
                
                # Calculate the bin width as area.
                # We need to convert the distance between the bins into units of 
                # area, since a smalp with a given radius will take up 
                # an area amount of bilayer lipid. 
                # Since the area equation involves a square, the difference
                # in area is *not the same* from bin to bin. 
                # You can test this yourself:
                #   radius_to_area(51) - radius_to_area(50) = 317.3
                #   radius_to_area(29) - radius_to_area(28) = 179.1
                # the bins both differ by radius 1, but by different areas.
                bin_start_area = radius_to_area(bin_start_radius)
                bin_end_area = radius_to_area(bin_end_radius)
                bin_width_as_area = bin_end_area - bin_start_area
        
                # It doesn't matter that we use the probability from the radius
                # distribution since the probability is unitless.
                probability = pradius(bin_average_radius)
                
                # This is the estimated area under the area distribution curve.
                probability_density = probability * bin_width_as_area
                
                number = total_bilayer_area * probability_density
                
                number_of_smalps += number
        
        return number_of_smalps
        
    
    def test_calc_number_of_smalps(peptide_concentration):
        """ This function shows that there is a asymptotic relationship between
        the number of bins and the number of SMALPs calculated. The choice of 
        number of bins is thus a trade-off between accuracy and computational 
        cost. I'm not sure why small bin numbers bias in one direction though. """
        # As it stands, this isn't affecting the result at all.
        n_peptides = concentration_to_number(peptide_concentration, volume_of_solution) 
        print(f"n_peptides: {n_peptides}")
        
        # We multiply this number by two because the peptide occupies both leaflets.
        total_peptide_area = (peptide_area*2) * n_peptides
        print(f"Total peptide area: {total_peptide_area} Å^2")
        
        total_system_area = total_lipid_area + total_peptide_area
        print(f"Total system area: {total_system_area} Å^2")
        
        smallest_radius_to_consider = float(smalp_core_radius_mean - (smalp_core_radius_std_dev * 4))
        largest_radius_to_consider = float(smalp_core_radius_mean + (smalp_core_radius_std_dev * 4))
        num_bins_list = [1000 * i for i in range(1,50,5)]
        
        n_smalps_list =[]
        for num_bins in num_bins_list:
            print(f"\nNumber of bins: {num_bins}")
            
            # Calculate the radius of each bin based on the number of bins.
            bin_ranges = np.linspace(smallest_radius_to_consider, largest_radius_to_consider, num_bins + 1)
            
            n_smalps = calc_number_of_smalps(bin_ranges, total_system_area)
            print(f"Number of SMALPs: {np.format_float_scientific(n_smalps)}")
            
            n_smalps_list.append(n_smalps)
        
        plt.plot(num_bins_list, n_smalps_list)
        plt.xlabel("Number of bins")
        plt.ylabel("Number of SMALPs calculated")
        plt.show()
        

        
    def calc_number_of_smalps_simple(radius, peptide_concentration):
        n_peptides = concentration_to_number(peptide_concentration, volume_of_solution) 
        print(f"n_peptides: {np.format_float_scientific(n_peptides)}")
        
        # We multiply this number by two because the peptide occupies both leaflets.
        total_peptide_area = (peptide_area*2) * n_peptides
        print(f"Total peptide area: {total_peptide_area} Å^2")
        
        total_system_area = total_lipid_area + total_peptide_area
        print(f"Total system area: {total_system_area} Å^2")
        
        number_of_smalps = total_system_area / (radius_to_area(radius)*2)
        print(f"Number of SMALPs: {np.format_float_scientific(number_of_smalps)}")
        
        print(f"Peptides per SMALP: {n_peptides / number_of_smalps}")
        return number_of_smalps
    
    
    # Define simulate_occupancy function
    def simulate_occupancy_binned(peptide_concentration, num_simulations=100):
        
        # TODO: fix the simulation function.
        
        # Determine the number of bins. By binning we assume  that a simulation
        # with fewer SMALPs will resemble one with more, as long as we get the
        # ratio of SMALPs to peptides right.
        num_bins = int(1e3)  # Adjust to computation needs. Smaller = faster.
        
        # Calculate the radius of each bin based on the number of bins.
        smallest_radius_to_consider = float(smalp_core_radius_mean - (smalp_core_radius_std_dev * 5))
        largest_radius_to_consider = float(smalp_core_radius_mean + (smalp_core_radius_std_dev * 5))
        bin_ranges = np.linspace(smallest_radius_to_consider, largest_radius_to_consider, num_bins + 1)
        
        # As it stands, this isn't affecting the result at all.
        n_peptides = concentration_to_number(peptide_concentration, volume_of_solution) 
        print(f"n_peptides: {n_peptides}")
        
        # We multiply this number by two because the peptide occupies both leaflets.
        total_peptide_area = (peptide_area*2) * n_peptides
        print(f"Total peptide area: {total_peptide_area} Å^2")
        
        total_system_area = total_lipid_area + total_peptide_area
        print(f"Total system area: {total_system_area} Å^2")
        
        n_smalps = calc_number_of_smalps(bin_ranges, total_system_area)
        
        n_smalps_binned = round(total_system_area / np.sum([pradius(r) * (radius_to_area(r)*2) for r in bin_ranges]))
        print(f"n_smalps_binned: {n_smalps_binned}")
        

        
        # Ref section 5.1.1, number 2
        # Determine the number of rows needed for each matrix based on Pradius(r) * Nliposomes
        # num_rows_per_matrix = [int(pradius(radius) * n_smalps) for radius in bin_ranges]
        
        
        


    
        


        # Calculate the cumulative area of all the SMALPs in each bin
        # This looks like a normal distribution, but it should be skewed
        # to the right because larger SMALPs contain more lipid, even though
        # There are fewer of them (assuming the size is normally distributed).
        bin_areas = np.diff(smalp_core_radius_distribution.cdf(bin_ranges))
        # total_system_area = np.sum(bin_areas)

    
        # Create arrays to track occupancy of each SMALP.
        smalp_bin_occupancies = np.zeros(num_bins)
    
        # Simulate occupancy for the given peptide concentration with tqdm progress bar
        for _ in tqdm(range(num_simulations), desc="Simulating occupancy"):
            random_bin = np.random.choice(num_bins, p=bin_areas / np.sum(bin_areas))
            smalp_bin_occupancies[random_bin] += 1

    
        # Calculate occupancy counts
        empty_count = np.count_nonzero(smalp_bin_occupancies == 0)
        monomer_count = np.count_nonzero(smalp_bin_occupancies == 1)
        dimer_count = np.count_nonzero(smalp_bin_occupancies == 2)
        trimer_plus_count = num_simulations - (empty_count + monomer_count + dimer_count)
    
        return [empty_count, monomer_count, dimer_count, trimer_plus_count]
    
    # # Define to_occupied_percents function
    # def to_occupied_percents(occupancy_counts):
    #     total_occupied = np.sum(occupancy_counts[1:])
    #     return [count / total_occupied * 100 for count in occupancy_counts[1:]]
    
    # # Simulate occupancy for the given concentrations
    # occupancy_curve = [to_occupied_percents(simulate_occupancy_binned(concentration)) for concentration in peptide_concentrations]
    
    # occupancy_curve_fixed = [[],[],[]]
    # for pl_simulation in occupancy_curve:
    #     [occupancy_curve_fixed[i].append(value) for i, value in enumerate(pl_simulation)]
    
    # # Plot the simulated occupancy curves
    # pl_ratios = pl_ratios[::-1]  # Reverse pl_ratios for the plot
    # labels = ['Monomer', 'Dimer', 'Trimer+']
    # plt.figure()
    # for i, label in enumerate(labels):
    #     plt.plot(pl_ratios, occupancy_curve_fixed[i], label=label)
    # plt.xlabel('P/L')
    # plt.ylabel('Percent (%)')
    # plt.title('Simulated Occupancy for Ideal Monomers')
    # plt.legend()
    # plt.show()
    
    # test_calc_number_of_smalps(1.5e-6)
    # calc_number_of_smalps_simple(38,1.5e-6)


if __name__ == "__main__":
    main()