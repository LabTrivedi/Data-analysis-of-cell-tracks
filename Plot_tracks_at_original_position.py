import open3d as o3d
import numpy as np
import math as m
import matplotlib.pyplot as plt
from statistics import mean
import pandas as pd
import os
import sys

def classify_cell(csv, cell, thr, mintrac=24, tol=2):
    cat = 0
    time = np.where(csv[:, 1] == cell)[0]
    if len(time) >= mintrac:
        if csv[time[0], 9] < thr:
            pureMin = True
            outPure = 0
            moveMin = False
            outMove = 0
            for ch1 in csv[time[1:], 9]:
                if pureMin or moveMin:
                    if ch1 >= thr:
                        outPure += 1
                        if outPure > tol: # If more than 'tol' points are over thr1, the track cannot be pure
                            pureMin = False
                            moveMin = True
                    else:
                        if moveMin:
                            outMove += 1
                            if outMove > tol:
                                moveMin = False
            if np.max(csv[time[1:-tol], 9]) < thr:
                pureMin = True
                moveMin = False
            if pureMin:
                cat = 4
            if moveMin:
                cat = 2
        elif csv[time[0], 9] > thr:
            pureMax = True
            outPure = 0
            moveMax = False
            outMove = 0
            for ch1 in csv[time[1:], 9]:
                if pureMax or moveMax:
                    if ch1 <= thr:
                        outPure += 1
                        if outPure > tol:  # If more than 'tol' points are over thr1, the track cannot be pure
                            pureMax = False
                            moveMax = True
                    else:
                        if moveMax:
                            outMove += 1
                            if outMove > tol:
                                moveMax = False
            if np.min(csv[time[1:-tol], 9]) > thr:
                pureMax = True
                moveMax = False
            if pureMax:
                cat = 3
            if moveMax:
                cat = 1
    return cat
def get_subfolders(directory):
    subfolders = [f.name for f in os.scandir(directory) if f.is_dir()]
    return subfolders

# PARAMETERS TO MODIFY
directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data/"
save_directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data_plots/Tracks/"
printInFolder = True
mintrac = 60

# CODE THAT SHOULD NOT BE MODIFIED
Cats = [0, 1, 2, 3, 4]
colors = ["blue", "orange", "green", "red", "purple"]
folders = get_subfolders(directory)

for folder in folders:
    subfolders = get_subfolders(directory + folder + '/')
    if not os.path.isdir(save_directory + folder):
        os.mkdir(save_directory + folder)
    for subfolder in subfolders:
        print('Folder ' + folder + ': ' + subfolder)

        # LOAD DATA FROM SERVER
        print('Loading data')

        csv = np.genfromtxt(directory + folder + "/" + subfolder + "/cells_no_div_no_jumps.csv", delimiter=",")
        n_size = np.size(csv, 0)
        bra = np.genfromtxt(directory + folder + "/" + subfolder + "/Tmin_Tplus_threshold.csv", delimiter=",")
        thr1 = bra[0]
        thr2 = bra[1]
        thr = (thr1 + thr2) / 2.0
        print('Data loaded')

        # MEASURE THE ABSOLUTE AND RELATIVE NUMBERS OF BRA/T
        max_cell = int(max(csv[1:n_size, 1]))

        fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(8, 12))
        plt.subplots_adjust(wspace=0.3, hspace=0.3)
        xx = 0
        yy = 0
        for Cat in Cats:
            yy += 1
            if yy == 2:
                yy = 0
                xx += 1
            for cell in range(max_cell+1):
                time = np.where(csv[:, 1] == cell)[0]
                if len(time) >= mintrac:
                    cat = classify_cell(csv, cell, thr)
                    time = np.where(csv[:, 1] == cell)[0]
                    if Cat != cat:
                        ax[xx, yy].plot(csv[time[:-1], 4], csv[time[:-1], 5], color='black', alpha=0.25)
                        ax[xx, yy].arrow(csv[time[-2], 4], csv[time[-2], 5],
                                         csv[time[-1], 4] - csv[time[-2], 4],
                                         csv[time[-1], 5] - csv[time[-2], 5],
                                         color='black', alpha=0.25, head_width=1.0)
                    else:
                        ax[xx, yy].plot(csv[time[:-1], 4], csv[time[:-1], 5], color=colors[Cat])
                        ax[xx, yy].arrow(csv[time[-2], 4], csv[time[-2], 5],
                                         csv[time[-1], 4] - csv[time[-2], 4],
                                         csv[time[-1], 5] - csv[time[-2], 5],
                                         color=colors[Cat], head_width=1.0)
            ax[xx, yy].set_aspect('equal')
        ax[0, 0].axis('off')

        if printInFolder:
            plt.savefig(save_directory + folder + '/' + subfolder + '_arrows_at_origin.png', dpi=350)
            plt.savefig(save_directory + folder + '/' + subfolder + '_arrows_at_origin.pdf')
        else:
            plt.show()
        plt.close(fig=fig)

