import open3d as o3d
import numpy as np
import math as m
import matplotlib.pyplot as plt
from statistics import mean
import pandas as pd
import os
import sys
import random
from pySankey.sankey import sankey

def get_subfolders(directory):
    subfolders = [f.name for f in os.scandir(directory) if f.is_dir()]
    return subfolders

# PARAMETERS TO MODIFY
directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data/"
save_directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data_plots/Sankey_diagram/"
printInFolder = True
mintrac = 24

# CODE THAT SHOULD NOT BE MODIFIED
folders = get_subfolders(directory)
BraNames = ['T-', 'Tmid', 'T+']

for folder in folders:
    df = pd.DataFrame(columns=['Initial_Ch1', 'End_Ch1', 'Units'])
    C = np.zeros((3, 3))
    subfolders = get_subfolders(directory + folder + '/')
    for subfolder in subfolders:
        print('Folder ' + folder + ': ' + subfolder)

        # LOAD DATA FROM SERVER
        print('Loading data')
        csv = np.genfromtxt(directory + folder + '/' + subfolder + "/cells_no_div_no_jumps.csv", delimiter=",")
        n_size = np.size(csv, 0)
        bra = np.genfromtxt(directory + folder + '/' + subfolder + "/Tmin_Tplus_threshold.csv", delimiter=",")
        thr1 = bra[0]
        thr2 = bra[1]
        thr = (thr1 + thr2) / 2.0
        print('Data loaded')

        # MEASURE THE ABSOLUTE AND RELATIVE NUMBERS OF BRA/T
        max_cell = int(max(csv[1:n_size, 1]))

        bra = np.genfromtxt(directory + folder + '/' + subfolder + "/Tmin_Tplus_threshold.csv", delimiter=",")
        for cell in range(max_cell + 1):
            times = np.where(csv[:, 1] == cell)[0]
            if len(times) >= mintrac:
                ch1i = np.mean(csv[times[0:6], 9])
                ch1f = np.mean(csv[times[-6:], 9])

                if ch1i < bra[0]:
                    ini = 0
                elif ch1i > bra[1]:
                    ini = 2
                else:
                    ini = 1

                if ch1f < bra[0]:
                    end = 0
                elif ch1f > bra[1]:
                    end = 2
                else:
                    end = 1

                C[ini, end] += 1

    for i in range(3):
        for j in range(3):
            if C[i, j] > 0.0:
                df_tmp = pd.DataFrame([(BraNames[i], BraNames[j], int(C[i, j]))], columns=['Initial_Ch1', 'End_Ch1', 'Units'])
                df = pd.concat([df, df_tmp], ignore_index=True)
    print(df)

    # Create a Sankey diagram
    s = sankey(left=df['Initial_Ch1'], right=df['End_Ch1'], leftWeight=df['Units'], rightWeight=df['Units'],
               fontsize=12)

    plt.gcf().set_size_inches((6, 6))
    plt.title('Bra/T changes')

    # Estimate node counts based on the input data
    node_counts_ini = {node: sum(df['Units'][df['Initial_Ch1'] == node]) for node in set(df['Initial_Ch1'])}
    node_counts_end = {node: sum(df['Units'][df['End_Ch1'] == node]) for node in set(df['End_Ch1'])}

    # Print node counts
    ini_counts = 0
    for node, count in node_counts_ini.items():
        ini_counts += count

    end_counts = 0
    for node, count in node_counts_end.items():
        end_counts += count

    plt.text(2, (node_counts_ini['T-']*0.5 - 10), node_counts_ini['T-'], fontsize=12)
    plt.text(2, (node_counts_ini['T-'] + node_counts_ini['Tmid']*0.5 ), node_counts_ini['Tmid'], fontsize=12)
    plt.text(2, (node_counts_ini['T-'] + node_counts_ini['Tmid'] + node_counts_ini['T+'] * 0.5 + 10), node_counts_ini['T+'], fontsize=12)

    plt.text(end_counts/4.2, (node_counts_end['T-'] * 0.5 - 10), node_counts_end['T-'], fontsize=12)
    plt.text(end_counts/4.2, (node_counts_end['T-'] + node_counts_end['Tmid'] * 0.5), node_counts_end['Tmid'], fontsize=12)
    plt.text(end_counts/4.2, (node_counts_end['T-'] + node_counts_end['Tmid'] + node_counts_end[ 'T+'] * 0.5 + 10), node_counts_end['T+'], fontsize=12)

    if printInFolder:
        plt.savefig(save_directory + folder + '_sankey_diagram_total_cells.pdf')
        plt.savefig(save_directory + folder + '_sankey_diagram_total_cells.png', dpi=350)
    else:
        plt.show()
    plt.close()