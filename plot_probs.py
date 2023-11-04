#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 18:14:29 2023

@author: justin

This script contains useful functions for plotting probabilities in various layouts.

"""
import numpy as np
import matplotlib.pyplot as plt

def plot_probability_bar(y_values, title=None):
    """
    p is a list of 3 probabilities in the order of p1, p2, p3+
    TODO: Add values inside or outside the bar, depending on the value.
    """
    # Check if there are exactly 3 y values
    if len(y_values) != 3:
        raise ValueError("Input list must contain exactly 3 y values.")

    # Labels for the bars
    category_labels = ["$p_1$", "$p_2$", "$p_{3+}$"]

    # X-axis positions for the bars
    x = range(len(category_labels))

    # Create the bar chart
    plt.bar(x, y_values, tick_label=category_labels, color="skyblue", edgecolor="black")
    plt.ylabel("Probability")

    if title:
        plt.title(title)

    for i, v in enumerate(y_values):
        if v > 0.9*max(y_values):
            plt.text(i, v - 0.02*max(y_values), f'{v:.2f}', ha='center', va='top', fontsize=12)  # Adjust the position
        else:
            plt.text(i, v + 0.01*max(y_values), f'{v:.2f}', ha='center', va='bottom', fontsize=12)  # Default position

    plt.show()


def plot_probability_curves(pl, p1, p2, p3, title=None):
    """
    Plots the p1, p2, and p3 curves for a simulation or experimental data.
    """
    if len(pl) != len(p1) != len(p2) != len(p3):
        raise ValueError("Input lists must have the same length.")

    plt.plot(pl, p1, label='p1', marker="o", color="black")
    plt.plot(pl, p2, label='p2', marker="o", color="red")
    plt.plot(pl, p3, label='p3+', marker="o", color="green")

    plt.xlabel("P/L")
    plt.ylabel("Probability")

    if title:
        plt.title(title)

    plt.legend()
    plt.show()

def plot_probability_comparison(pl, pe, pm, pd, title=None):
    """
    Plots a comparison of one probability among the experimental data and the
    monomer and dimer simulations.
    """

    if len(pl) != len(pe) != len(pm) != len(pd):
        raise ValueError("Input lists must have the same length.")

    plt.plot(pl, pe, label='Experimental', marker="o", color='red')
    plt.plot(pl, pm, label='Monomer sim.', marker=".", color='lightgray')
    plt.plot(pl, pd, label='Dimer sim.', marker=".", color='black')

    # Fill between the lines
    plt.fill_between(pl, pd, color='black')
    plt.fill_between(pl, pd, pm, color='lightgray')

    plt.xlabel("P/L")
    plt.ylabel("Probability")

    if title:
        plt.title(title)

    plt.legend()
    plt.show()

def main():

    # Example usage
    pm = [0.024451410660000002, 0.05141065831, 0.02413793103]
    plot_probability_bar(pm, "Monomer, P/L = 0.02")



    # Example usage
    peptide_concentrations = np.array([1.5, 0.75, 0.46875, 0.375, 0.1875, 0.09375])*1e-6
    lipid_concentration = 75e-6
    pl = peptide_concentrations/lipid_concentration
    pe1=[0.024451410660000002, 0.04580645161000001, 0.061494252869999995, 0.0686, 0.07292225201]
    pe2=[0.05141065831, 0.041612903230000005, 0.0316091954, 0.0264, 0.023056300270000003]
    pe3=[0.02413793103, 0.012580645159999999, 0.0068965517240000005, 0.005, 0.004021447721]

    # Trim the 6th value because it's not in the experimental data.
    pl, pe1, pe2, pe3 = pl[:5], pe1[:5], pe2[:5], pe3[:5]

    plot_probability_curves(pl, pe1, pe2, pe3, title="Experimental")



    # Example usage
    p1e = [0.024451410660000002, 0.04580645161000001, 0.061494252869999995, 0.0686, 0.07292225201]
    p1m = [0.3151, 0.3509, 0.3097, 0.2701, 0.1729, 0.0965]
    p1d = [0.0004, 0.0009, 0.0025, 0.0036, 0.0179, 0.1128]

    # Trim the 6th value because it's not in the experimental data.
    pl, p1e, p1m, p1d = pl[:5], p1e[:5], p1m[:5], p1d[:5]

    plot_probability_comparison(pl, p1e, p1m, p1d, title="p1")

if __name__ == "__main__":
    main()