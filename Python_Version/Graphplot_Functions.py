#!/usr/bin/python3


"""
File contains functions to plot graphs
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def Consensus_plot(Consensus):

    # Take the consensus frequency of each nucleotide from the CSV file and plot a bar graph for all residue positions
    df_cons = pd.read_csv((Consensus + '.csv'), sep=',')  # Place csv file in a dataframe
    cons_x = df_cons['Alignment_Position']  # Get x values
    cons_y = df_cons['Consensus_Frequency%']  # Get y values
    fig, ax = plt.subplots(figsize=[15, 8])  # Generate graph layout and size
    cons_bars = ax.bar(cons_x, cons_y, color='navy')  # Plot bar graph with initialized x- and y-axis values
    ax.set_xticks(np.arange(0, len(cons_x) + 1, 50))  # Set the x-axis tick marks
    ax.margins(x=0.01)  # Remove margin gaps on left and right side of x-axis
    labels = ax.get_xticks().tolist()  # Put the x-axis tick values in a variable
    ax.set_xticklabels(labels, rotation=90)  # Rotate the x-axis tick values by 90 degrees
    for bar in cons_bars:  # Loop over each bar in the graph and remove gaps in between
        bar.set_width(1)
    ax.set_xlabel('Residue Position', fontsize=12, weight='bold')  # Set x-axis label
    ax.set_ylabel('Consensus Residue Frequency', fontsize=12, weight='bold')  # Set y-axis label
    ax.set_title('Frequency of Consensus residues in sequence\n(as %)', fontsize=14, weight='bold')  # Set graph title
    plt.savefig((Consensus + '.png'))  # Save the graph as a PNG under the same name as CSV file




def RE_plot(Entropy, Position, Max_RE, Cumulative_RE):

    # Generate the plot for the Max and Cumulative Relative entropy at each position across all organisms/strains
    RE_x = Position   # Store x-axis values
    RE_y1 = Max_RE  # Store primary y-axis values
    RE_y2 = Cumulative_RE  # Store secondary y-axis values
    fig, ax1 = plt.subplots(figsize=[15, 8])  # Generate graph layout, size with axes
    ax1.plot(RE_x, RE_y1, color='r', linewidth=0.5)  # Plot x and primary y-axis values, y is a line graph
    ax1.set_ylim(bottom=0)  # Ensure that the primary y-axis 0 is at the x-axis
    ax1.set_xticks(np.arange(0, len(RE_x) + 1, 50))  # Set the x-axis tick marks
    labels = ax1.get_xticks().tolist()  # Put the x-axis tick values in a variable
    ax1.set_xticklabels(labels, rotation=90)  # Rotate the x-axis tick values by 90 degrees
    ax2 = ax1.twinx()  # Generate graph layout for secondary y-axis, having the same x-axis
    bars = ax2.bar(RE_x, RE_y2, color='grey')  # Plot x and secondary y-axis values, y is a bar graph
    for bar in bars:  # Loop over each bar in the graph and remove gaps in between
        bar.set_width(1)
    ax2.set_ylim(bottom=0)  # Make sure the y-tick 0 coincides with x-axis
    ax1.margins(x=0.01)  # Remove margin gaps on left and right side of x-axis
    ax1.set_xlabel('Residue Position', fontsize=12, weight='bold')  # Set x-axis label
    ax1.set_ylabel('Maximum Relative Entropy', fontsize=12, weight='bold')  # Set primary y-axis label
    ax2.set_ylabel('Cumulative Relative Entropy', fontsize=12, weight='bold', rotation=270)  # Set y-axis label
    ax2.yaxis.set_label_coords(1.05, 0.5)
    ax1.set_title('Maximum and Cumulative Relative Entropy\nat each position', fontsize=14, weight='bold')  # Set graph title
    plt.savefig((Entropy + '.png'))  # Save the graph as a PNG under the same name as Relative Entropy file
