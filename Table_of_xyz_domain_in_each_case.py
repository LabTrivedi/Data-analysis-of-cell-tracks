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
save_directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data_plots/Domain_xyz/"
mintrac = 2
printInFolder = True

# CODE THAT SHOULD NOT BE MODIFIED
folders = get_subfolders(directory)

data_colums = ['Time', 'Condition', 'X', 'Y', 'Z', 'relX', 'relY', 'relZ', 'diffX', 'diffY', 'diffZ', 'reldiffX', 'reldiffY', 'reldiffZ']
df = pd.DataFrame(columns=data_colums)

if not os.path.exists(directory + "/table_of_xyz_domain.csv"):
    for folder in folders:
        cond = folder.split('_')[0]
        time = folder.split('_')[1]

        if cond == 'a2i':
            cond = '2i'

        subfolders = get_subfolders(directory + folder + '/')
        for subfolder in subfolders:
            print('Folder ' + folder + ': ' + subfolder)

            csv = np.genfromtxt(directory + folder + "/" + subfolder + "/cells_no_div_no_jumps.csv", delimiter=",")
            n_size = np.size(csv, 0)
            bra = np.genfromtxt(directory + folder + "/" + subfolder + "/Tmin_Tplus_threshold.csv", delimiter=",")
            thr1 = bra[0]
            thr2 = bra[1]

            # MEASURE THE ABSOLUTE AND RELATIVE NUMBERS OF BRA/T
            max_cell = int(max(csv[1:n_size, 1]))
            max_spot = int(max(csv[1:n_size, 2]))
            max_time = int(max(csv[1:n_size, 7]))

            devX = []
            devY = []
            devZ = []

            for cell in range(max_cell+1):
                times = np.where(csv[:, 1] == cell)[0]
                if len(times) >= mintrac:
                    meanX = np.mean(csv[times, 4])
                    meanY = np.mean(csv[times, 5])
                    meanZ = np.mean(csv[times, 6])

                    dX = np.linalg.norm(csv[times, 4] - meanX)
                    dY = np.linalg.norm(csv[times, 5] - meanY)
                    dZ = np.linalg.norm(csv[times, 6] - meanZ)

                    devX.append(dX)
                    devY.append(dY)
                    devZ.append(dZ)

            dvX = mean(devX)
            dvY = mean(devY)
            dvZ = mean(devZ)

            totdv = dvX + dvY + dvZ

            max_time = int(max(csv[1:n_size, 7]))
            for t in range(max_time+1):
                times = np.where(csv[:, 7] == t)[0]
                if t == 0:
                    diffX = np.max(csv[times, 4]) - np.min(csv[times, 4])
                    diffY = np.max(csv[times, 5]) - np.min(csv[times, 5])
                    diffZ = np.max(csv[times, 6]) - np.min(csv[times, 6])
                else:
                    diffX = max(np.max(csv[times, 4]) - np.min(csv[times, 4]), diffX)
                    diffY = max(np.max(csv[times, 5]) - np.min(csv[times, 5]), diffY)
                    diffZ = max(np.max(csv[times, 6]) - np.min(csv[times, 6]), diffZ)

            totdiff = diffX + diffY + diffZ

            df_tmp = pd.DataFrame([(time, cond, dvX, dvY, dvZ, dvX/totdv, dvY/totdv, dvZ/totdv, diffX, diffY, diffZ, diffX/totdiff, diffY/totdiff, diffZ/totdiff)], columns=data_colums)
            df = pd.concat([df, df_tmp], ignore_index=True)

    sorted_df = df.sort_values(by=['Time', 'Condition'])
    sorted_df.to_csv(directory + "/table_of_xyz_domain.csv", index=False)

colors = dict()
colors['SL'] = '#838B8B'
colors['1i'] = '#00688B'
colors['2i'] = '#8FBC8F'

shapes = dict()
shapes['48h'] = "o"
shapes['72h'] = "s"
shapes['96h'] = "^"

df = pd.read_csv(directory + 'table_of_xyz_domain.csv')

Time = df['Time'].to_list()
Cond = df['Condition'].to_list()
Shap = [shapes[x] for x in Time]
Cols = [colors[x] for x in Cond]

print(df)

# Relative XYZ
relX = df['relX'].to_numpy()
relY = df['relY'].to_numpy()
relZ = df['relZ'].to_numpy()

idx_48h = df[df['Time'] == '48h'].index.to_numpy()
idx_72h = df[df['Time'] == '72h'].index.to_numpy()
idx_96h = df[df['Time'] == '96h'].index.to_numpy()

col_48h = [colors[df['Condition'][cc]] for cc in idx_48h]
col_72h = [colors[df['Condition'][cc]] for cc in idx_72h]
col_96h = [colors[df['Condition'][cc]] for cc in idx_96h]


print('Print relative positions')

plt.scatter(0*np.ones(len(idx_48h)) - 0.1, relX[idx_48h], color=col_48h, marker="o")
plt.scatter(0*np.ones(len(idx_72h)) + 0.0, relX[idx_72h], color=col_72h, marker="s")
plt.scatter(0*np.ones(len(idx_96h)) + 0.1, relX[idx_96h], color=col_96h, marker="^")
plt.text(-0.05, 0.025, 'relX')

plt.scatter(1*np.ones(len(idx_48h)) - 0.1, relY[idx_48h], color=col_48h, marker="o")
plt.scatter(1*np.ones(len(idx_72h)) + 0.0, relY[idx_72h], color=col_72h, marker="s")
plt.scatter(1*np.ones(len(idx_96h)) + 0.1, relY[idx_96h], color=col_96h, marker="^")
plt.text(0.95, 0.025, 'relY')

plt.scatter(2*np.ones(len(idx_48h)) - 0.1, relZ[idx_48h], color=col_48h, marker="o")
plt.scatter(2*np.ones(len(idx_72h)) + 0.0, relZ[idx_72h], color=col_72h, marker="s")
plt.scatter(2*np.ones(len(idx_96h)) + 0.1, relZ[idx_96h], color=col_96h, marker="^")
plt.text(1.95, 0.025, 'relZ')

plt.xticks([])
plt.yticks([0.05, 0.3, 0.55])
if printInFolder:
    plt.savefig(save_directory + 'scatter_plot_by_time_and_condition_relXYZ.pdf')
    plt.savefig(save_directory + 'scatter_plot_by_time_and_condition_relXYZ.png', dpi=350)
else:
    plt.show()

print('Print difference positions')

# Relative diff XYZ
plt.cla()
plt.clf()

diffrelX = df['reldiffX'].to_numpy()
diffrelY = df['reldiffY'].to_numpy()
diffrelZ = df['reldiffZ'].to_numpy()

plt.scatter(0.0*np.ones(len(idx_48h)) - 0.1, diffrelX[idx_48h], color=col_48h, marker="o")
plt.scatter(0.0*np.ones(len(idx_72h)) + 0.0, diffrelX[idx_72h], color=col_72h, marker="s")
plt.scatter(0.0*np.ones(len(idx_96h)) + 0.1, diffrelX[idx_96h], color=col_96h, marker="^")
plt.text(-0.15, 0.025, 'reldiffX')


plt.scatter(1.0*np.ones(len(idx_48h)) - 0.1, diffrelY[idx_48h], color=col_48h, marker="o")
plt.scatter(1.0*np.ones(len(idx_72h)) + 0.0, diffrelY[idx_72h], color=col_72h, marker="s")
plt.scatter(1.0*np.ones(len(idx_96h)) + 0.1, diffrelY[idx_96h], color=col_96h, marker="^")
plt.text(0.85, 0.025, 'reldiffY')

plt.scatter(2.0*np.ones(len(idx_48h)) - 0.1, diffrelZ[idx_48h], color=col_48h, marker="o")
plt.scatter(2.0*np.ones(len(idx_72h)) + 0.0, diffrelZ[idx_72h], color=col_72h, marker="s")
plt.scatter(2.0*np.ones(len(idx_96h)) + 0.1, diffrelZ[idx_96h], color=col_96h, marker="^")
plt.text(1.85, 0.025, 'reldiffZ')

plt.xticks([])
plt.yticks([0.05, 0.3, 0.55])
if printInFolder:
    plt.savefig(save_directory + 'scatter_plot_by_time_and_condition_reldiffXYZ.pdf')
    plt.savefig(save_directory + 'scatter_plot_by_time_and_condition_reldiffXYZ.png', dpi=350)
else:
    plt.show()