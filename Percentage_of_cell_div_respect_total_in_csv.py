import open3d as o3d
import numpy as np
import math as m
import matplotlib.pyplot as plt
from statistics import mean
import pandas as pd
import os
import sys

def track_division(times, chan1, csv, edges, cell0):
    acc_times_1 = []
    acc_chan1_1 = []
    acc_times_2 = []
    acc_chan1_2 = []

    acc_times_1.append(csv[cell0, 5])
    acc_chan1_1.append(csv[cell0, 6])
    acc_times_2.append(csv[cell0, 5])
    acc_chan1_2.append(csv[cell0, 6])

    c0 = np.where(edges[:, 0] == csv[cell0, 1])[0]
    cell1 = int(edges[c0[0], 1])
    cell2 = int(edges[c0[1], 1])

    # CELL 1
    if csv[cell1, 2] == 2.0:
        while csv[cell1, 2] == 2.0:
            acc_times_1.append(csv[cell1, 5])
            acc_chan1_1.append(csv[cell1, 6])
            c1 = np.where(edges[:, 0] == csv[cell1, 1])[0][0]
            cell1 = int(edges[c1, 1])
        acc_times_1.append(csv[cell1, 5])
        acc_chan1_1.append(csv[cell1, 6])
        times.append(acc_times_1)
        chan1.append(acc_chan1_1)
    else:
        acc_times_1.append(csv[cell1, 5])
        acc_chan1_1.append(csv[cell1, 6])
        times.append(acc_times_1)
        chan1.append(acc_chan1_1)

    # CELL 2
    if csv[cell2, 2] == 2.0:
        while csv[cell2, 2] == 2.0:
            acc_times_2.append(csv[cell2, 5])
            acc_chan1_2.append(csv[cell2, 6])
            c1 = np.where(edges[:, 0] == csv[cell2, 1])[0][0]
            cell2 = int(edges[c1, 1])
        acc_times_2.append(csv[cell2, 5])
        acc_chan1_2.append(csv[cell2, 6])
        times.append(acc_times_2)
        chan1.append(acc_chan1_2)
    else:
        acc_times_2.append(csv[cell2, 5])
        acc_chan1_2.append(csv[cell2, 6])
        times.append(acc_times_2)
        chan1.append(acc_chan1_2)
def get_subfolders(directory):
    subfolders = [f.name for f in os.scandir(directory) if f.is_dir()]
    return subfolders

# PARAMETERS TO MODIFY
directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data/"
save_directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data_plots/Cell_Division_Percentage/"

# CODE THAT SHOULD NOT BE MODIFIED
# Define column names
column_names = ['Mean', 'Median']

# Create an empty DataFrame with predefined column names
df = pd.DataFrame(columns=column_names)

folders = get_subfolders(directory)

for folder in folders:
    subfolders = get_subfolders(directory + folder + '/')
    for subfolder in subfolders:
        print('Folder ' + folder + ': ' + subfolder)

        # LOAD DATA FROM SERVER
        print('Loading data')

        csv_list = []

        if os.path.exists(directory + folder + '/' + subfolder + "/FeatureAndTagTable-vertices.csv"):
            with open(directory + folder + '/' + subfolder +"/FeatureAndTagTable-vertices.csv", "r") as my_input_file:
                idx = 0
                for line in my_input_file:
                    line = line.split(",")
                    if idx > 2:
                        line_list = []
                        for i in range(len(line)):
                            linei = line[i].replace('"', '')
                            if i == len(line)-1:
                                line_list.append(float(linei[0:-1]))
                            else:
                                line_list.append(float(linei))
                        csv_list.append(line_list)
                    idx += 1
            csv = np.array(csv_list)
            csv_size = np.size(csv, 0)

            edge_list = []

            with open(directory + folder + '/' + subfolder +"/FeatureAndTagTable-edges.csv", "r") as my_input_file:
                idx = 0
                for line in my_input_file:
                    line = line.split(",")
                    if idx > 2:
                        line_list = []
                        for i in range(len(line)):
                            linei = line[i].replace('"', '')
                            if i == 0:
                                splt = str(linei).split(' ? ')
                                line_list.append(float(splt[0]))
                                line_list.append(float(splt[1]))
                            elif i == len(line) - 1:
                                line_list.append(float(str(linei[0:-1])))
                            else:
                                line_list.append(float(str(linei)))
                        edge_list.append(line_list)
                    idx += 1
            edges = np.array(edge_list)
            edge_size = np.size(edges, 0)

            print('Data loaded')

            # ANALYZE CELLS BRA/T AFTER DIVISION
            max_time = int(np.max(csv[:, 5]))
            percent = np.zeros((max_time + 1, 2))

            divs = np.where(csv[:, 2] == 3.0)[0]

            for time in range(max_time+1):
                print('Time ' + str(time) + ' of ' + str(max_time))
                tot_cells = np.where(csv[:, 5] == float(time))[0]
                div_cells = np.where(csv[divs, 5] == float(time))[0]
                # print(len(div_cells))
                # print(len(tot_cells))
                percent[time, 0] = time*5/60
                percent[time, 1] = float(len(div_cells))/float(len(tot_cells))

            # Create a dictionary with data for each row
            data = {
                'Mean': np.mean(percent[:, 1]),
                'Median': np.median(percent[:, 1]),
            }

            # Append the data to the DataFrame
            df.loc[folder] = data
        else:
            print(directory + folder + '/' + subfolder + "/FeatureAndTagTable-vertices.csv not found")

print(df)
df.to_csv(save_directory + 'cell_division_percentage_mean_and_median.csv')