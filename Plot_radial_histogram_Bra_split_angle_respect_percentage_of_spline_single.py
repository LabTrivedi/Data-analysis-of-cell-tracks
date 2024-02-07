import open3d as o3d
import numpy as np
import math as m
import matplotlib.pyplot as plt
from statistics import mean
import pandas as pd
import os
import sys
import json
def circular_hist(ax, x, bins=16, density=True, offset=0, gaps=True, alpha=1, color='blue'):
    """
    Produce a circular histogram of angles on ax.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.PolarAxesSubplot
        axis instance created with subplot_kw=dict(projection='polar').

    x : array
        Angles to plot, expected in units of radians.

    bins : int, optional
        Defines the number of equal-width bins in the range. The default is 16.

    density : bool, optional
        If True plot frequency proportional to area. If False plot frequency
        proportional to radius. The default is True.

    offset : float, optional
        Sets the offset for the location of the 0 direction in units of
        radians. The default is 0.

    gaps : bool, optional
        Whether to allow gaps between bins. When gaps = False the bins are
        forced to partition the entire [-pi, pi] range. The default is True.

    Returns
    -------
    n : array or list of arrays
        The number of values in each bin.

    bins : array
        The edges of the bins.

    patches : `.BarContainer` or list of a single `.Polygon`
        Container of individual artists used to create the histogram
        or list of such containers if there are multiple input datasets.
    """
    # Wrap angles to [-pi, pi)
    x = (x+np.pi) % (2*np.pi) - np.pi

    # Force bins to partition entire circle
    if not gaps:
        bins = np.linspace(-np.pi, np.pi, num=bins+1)

    # Bin data and record counts
    n, bins = np.histogram(x, bins=bins)

    # Compute width of each bin
    widths = np.diff(bins)

    # By default plot frequency proportional to area
    if density:
        # Area to assign each bin
        area = n / x.size
        # Calculate corresponding bin radius
        radius = (area/np.pi) ** .5
    # Otherwise plot frequency proportional to radius
    else:
        radius = n

    # Plot data on ax
    #patches = ax.bar(bins[:-1], radius, zorder=1, align='edge', width=widths, edgecolor='C0', fill=True, linewidth=1, alpha=alpha, color=color)
    patches = ax.bar(bins[:-1], radius, zorder=1, align='edge', width=widths, edgecolor=color, fill=True, linewidth=1, alpha=alpha, color=color)

    # Set the direction of the zero angle
    ax.set_theta_offset(offset)

    # Remove ylabels for area plots (they are mostly obstructive)
    # if density:
    #     maxrad = max(radius)
    #     print(maxrad)
    #     ax.set_yticks([maxrad*0.25, maxrad*0.5, maxrad*0.75])

    return n, bins, patches, max(radius)
def closestNode(x, y, midline):
    lm = len(midline)
    for i in range(lm):
        dist = (midline[i][0] - x)**2 + (midline[i][1] - y)**2
        if i == 0:
            I = 0
            minDist = dist
        else:
            if dist < minDist:
                I = i
                minDist = dist
    return I
def get_subfolders(directory):
    subfolders = [f.name for f in os.scandir(directory) if f.is_dir()]
    return subfolders

# 0 =
# 1 = cell_id
# 2 = n_spots
# 3 = r
# 4 = x
# 5 = y
# 6 = z
# 7 = t
# 8 = speed
# 9 = mean-ch1
# 10 = mean-ch2

# PARAMETERS TO MODIFY
directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data/"
save_directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data_plots/Radial_Histogram/"
printInFolder = True
dt = 36

# CODE THAT SHOULD NOT BE MODIFIED
aspect_ratio = 0.3467
colors = ['grey', 'lime', 'green']
folders = get_subfolders(directory)

for folder in folders:
    subfolders = get_subfolders(directory + folder + '/')

    fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
    fig.set_figheight(5)
    fig.set_figwidth(5)

    totCellsMin = []
    totCellsMax = []

    for subfolder in subfolders:
        print('Folder ' + folder + ': ' + subfolder)

        if os.path.exists(directory + folder + '/' + subfolder + '/MIP_XY_morpho_params.json'):
            # LOAD DATA FROM SERVER
            print('Loading data')
            csv = np.genfromtxt(directory + folder + '/' + subfolder + "/cells_no_div_no_jumps.csv", delimiter=",")
            n_size = np.size(csv, 0)
            print('Data loaded')

            # MEASURE THE ABSOLUTE AND RELATIVE NUMBERS OF BRA/T
            max_cell = int(max(csv[1:n_size, 1]))
            max_spot = int(max(csv[1:n_size, 2]))
            max_time = int(max(csv[1:n_size, 7]))
            Time = np.array(range(max_time+1))
            # Minimum, Average and Maximum for Bra/T-
            CellsMin = np.zeros((max_time + 1, 3))
            # Minimum, Average and Maximum for Bra/Tmid
            CellsMid = np.zeros((max_time + 1, 3))
            # Minimum, Average and Maximum for Bra/T+
            CellsMax = np.zeros((max_time + 1, 3))

            # PARAMETERS DERIVED FROM INPUTS
            bra = np.genfromtxt(directory + folder + '/' + subfolder + "/Tmin_Tplus_threshold.csv", delimiter=",")
            thr1 = bra[0]
            thr2 = bra[1]
            thr = (thr1+thr2)/2.0

            regMin = np.zeros((max_time + 1, 2))
            regMax = np.zeros((max_time + 1, 2))

            times = np.linspace(0, max_time, max_time+1)

            # Compute vector at extremes
            with open(directory + folder + '/' + subfolder + '/MIP_XY_morpho_params.json', 'r') as fcc_file:
                data = json.load(fcc_file)
                print('Data loaded')
                ld = len(data)

                for time in range(max_time + 1):
                    print('Time ' + str(time) + ' of ' + str(max_time))
                    ld = len(data[time]['midline'])
                    midline = np.zeros((ld, 2))
                    perline = np.zeros(ld)
                    tlength = 0.0
                    startY = data[time]['slice'][0]['start']
                    startX = data[time]['slice'][1]['start']
                    for i in range(ld):
                        midline[i, 1] = (startY + data[time]['midline'][i][0])*aspect_ratio
                        midline[i, 0] = (startX + data[time]['midline'][i][1])*aspect_ratio
                        if i >= 1:
                            segment = np.sqrt((midline[i, 0] - midline[i-1, 0])**2 + (midline[i, 1] - midline[i-1, 1])**2)
                            perline[i] = perline[i-1] + segment
                            tlength += segment
                    perline = perline/tlength

                    cellsMin = []
                    cellsMax = []
                    for cell in range(max_cell+1):
                        times = np.where(csv[:, 1] == cell)[0]
                        if len(times) > dt:
                            x = csv[times[0], 4]
                            y = csv[times[0], 5]
                            I = closestNode(x, y, midline)
                            vcell = np.array([csv[times[dt], 4] - csv[times[0], 4], csv[times[dt], 5] - csv[times[0], 5]])
                            if I == ld - 1:
                                I -= 1
                            vspli = np.array([midline[I+1, 0] - midline[I, 0], midline[I+1, 1] - midline[I, 1]])
                            theta = m.acos(np.dot(vcell, vspli)/(np.linalg.norm(vcell) * np.linalg.norm(vspli)))

                            if csv[times[0], 9] < thr1:
                                cellsMin.append(theta)
                                totCellsMin.append(theta)

                            elif csv[times[0], 9] > thr2:
                                cellsMax.append(theta)
                                totCellsMax.append(theta)

    if len(totCellsMin) > 0 and len(totCellsMax) > 0:
        n1, bins1, patches1, rad1 = circular_hist(ax, np.array(totCellsMin), bins=16, density=True, offset=0, gaps=True, alpha=0.5, color='grey')
        n2, bins2, patches2, rad2 = circular_hist(ax, np.array(totCellsMax), bins=16, density=True, offset=0, gaps=True, alpha=0.5, color='green')
        maxrad = max(rad1, rad2)

        ax.set_xlim(0, m.pi)
        ax.set_ylim(0, maxrad)
        ax.set_yticks([maxrad*0.25, maxrad*0.5, maxrad*0.75])
        ax.set_yticklabels([])
        if printInFolder:
            plt.savefig(save_directory + folder + '_radial_histogram_cell_angle_respect_closest_spline_single.pdf')
            plt.savefig(save_directory + folder + '_radial_histogram_cell_angle_respect_closest_spline_single.png', dpi=350)
        else:
            plt.show()