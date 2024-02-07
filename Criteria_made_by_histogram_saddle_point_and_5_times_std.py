import numpy as np
import math as m
import matplotlib.pyplot as plt
import statistics
from statistics import mean
import pandas as pd
import os
import sys

def get_subfolders(directory):
    subfolders = [f.name for f in os.scandir(directory) if f.is_dir()]
    return subfolders

# PARAMETERS TO MODIFY
directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data/"
printInFolder = True

# CODE THAT SHOULD NOT BE MODIFIED
folders = get_subfolders(directory)

for folder in folders:
    subfolders = get_subfolders(directory + folder + '/')
    for subfolder in subfolders:
        print('Folder ' + folder + ': ' + subfolder)
        fig, ax = plt.subplots()

        # LOAD DATA FROM SERVER
        print('Loading data')
        csv = np.genfromtxt(directory + folder + "/" + subfolder + "/cells_no_div_no_jumps.csv", delimiter=",")
        n_size = np.size(csv, 0)
        print('Data loaded')

        # ANALYZE TIME
        max_cell = int(max(csv[1:n_size, 1]))
        max_spot = int(max(csv[1:n_size, 2]))
        max_time = int(max(csv[1:n_size, 7]))

        ch1 = csv[1:n_size, 9]
        meanCh1 = np.mean(ch1)
        minCh1 = np.min(ch1)*1.01
        maxCh1 = 2.0*meanCh1 - minCh1

        weights = 4.0*np.ones_like(ch1) / float(len(ch1))
        #ax.hist(ch1, bins=200, weights=weights)
        histog = ax.hist(ch1, bins=200, cumulative=True, density=True)
        lh = len(histog[0])
        H = np.zeros((lh, 2))
        for i in range(lh):
            H[i, 1] = histog[0][i]
            H[i, 0] = (histog[1][i+1] + histog[1][i])/2.0
        H2D = np.zeros((lh-2, 2))
        for i in range(1, lh-1):
            H2D[i-1, 0] = H[i, 0]
            H2D[i-1, 1] = H[i+1, 1] + H[i-1, 1] - 2.0*H[i, 1]
        for i in range(0, len(H2D)-1):
            if i == 0:
                I = 0
                min_sl = H2D[i+1, 1] - H2D[i, 1]
            else:
                sl = H2D[i+1, 1] - H2D[i, 1]
                if sl < min_sl:
                    I = i
                    min_sl = sl

        Tmin = (H2D[I+1, 0] + H2D[I, 0])/2.0
        lam = np.mean(ch1)
        #Tmax = Tmin + 3.0*m.sqrt(lam)
        Tmax = Tmin + 5.0 * m.sqrt(lam)
        print(Tmin)
        print(Tmax)
        if printInFolder:
            np.savetxt(directory + folder + "/" + subfolder + '/Tmin_Tplus_threshold.csv', np.array([Tmin, Tmax]))

