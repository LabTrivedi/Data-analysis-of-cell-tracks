import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

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


aspect_ratio = 0.3467

# Function that transforms the total csv in a list of filtered cells
def create_list_of_filtered_times(csv, mintrac):
    n_size = np.size(csv, 0)
    total_list = []
    cell_list = []
    first_time = True
    for i in range(1, n_size):
        if csv[i, 2] >= mintrac:
            if first_time or csv[i, 1] == csv[i - 1, 1]:
                cell_list.append(csv[i, :])
                first_time = False
            else:
                total_list.append(cell_list)
                cell_list = []
                cell_list.append(csv[i, :])
    return total_list
def msd(cells, c, mintrac, m):
    sum = 0.0
    if m > 0:
        for k in range(0, mintrac-m):
            sum += (cells[c][k+m][4] - cells[c][k][4])**2 + (cells[c][k+m][5] - cells[c][k][5])**2 + (cells[c][k+m][6] - cells[c][k][6])**2
    return sum/float(mintrac-m)
def get_subfolders(directory):
    subfolders = [f.name for f in os.scandir(directory) if f.is_dir()]
    return subfolders

# PARAMETERS TO MODIFY
directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data/"
save_directory = "//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/all_tracking_data_plots/MSD/"
mintrac = 36
printInFolder = True

# CODE THAT SHOULD NOT BE MODIFIED
# Define column names
column_names = ['Total', 'T-', 'Tmid', 'T+', 'Condition']

# Create an empty DataFrame with predefined column names
df = pd.DataFrame(columns=column_names)

folders = get_subfolders(directory)

for folder in folders:
    subfolders = get_subfolders(directory + folder + '/')
    for subfolder in subfolders:
        print('Folder ' + folder + ': ' + subfolder)

        # LOAD DATA FROM SERVER
        print('Loading data')
        csv = np.genfromtxt(directory + folder + "/" + subfolder + "/cells_no_div_no_jumps.csv", delimiter=",")
        n_size = np.size(csv, 0)
        print('Data loaded')

        # MEASURE THE ABSOLUTE AND RELATIVE NUMBERS OF BRA/T
        max_cell = int(max(csv[1:n_size, 1]))
        max_spot = int(max(csv[1:n_size, 2]))
        max_time = int(max(csv[1:n_size, 7]))

        # PARAMETERS DERIVED FROM INPUTS
        bra = np.genfromtxt(directory + folder + "/" + subfolder + "/Tmin_Tplus_threshold.csv", delimiter=",")
        thr1 = bra[0]
        thr2 = bra[1]
        thr = (thr1+thr2)/2.0

        # PLOTS

        # Loop over cells
        cells = create_list_of_filtered_times(csv, mintrac)
        tau = np.linspace(0, mintrac-1, mintrac)
        msd0 = np.zeros(mintrac)
        av_msd = np.zeros(mintrac)
        av_msdmin = np.zeros(mintrac)
        av_msdmid = np.zeros(mintrac)
        av_msdmax = np.zeros(mintrac)
        min_msdmin = np.zeros(mintrac)
        max_msdmin = np.zeros(mintrac)
        min_msdmid = np.zeros(mintrac)
        max_msdmid = np.zeros(mintrac)
        min_msdmax = np.zeros(mintrac)
        max_msdmax = np.zeros(mintrac)
        min_msd = np.zeros(mintrac)
        max_msd = np.zeros(mintrac)
        n_min = 0.0
        n_mid = 0.0
        n_max = 0.0
        intens = np.zeros(len(cells))
        coeffi = np.zeros(len(cells))

        firstMin = np.full(mintrac, True)
        firstMid = np.full(mintrac, True)
        firstMax = np.full(mintrac, True)
        firstCell = np.full(mintrac, True)
        for j in range(len(cells)):
            cellj = np.array(cells[j])
            intens[j] = np.average(cellj[:, 9])

            for i in tau:
                I = int(i)
                msd0[I] = msd(cells, j, mintrac, I)
                av_msd[I] += msd0[I]
                if intens[j] < thr1:
                    av_msdmin[I] += msd0[I]
                    if firstMin[I]:
                        min_msdmin[I] = msd0[I]
                        max_msdmin[I] = msd0[I]
                        firstMin[I] = False
                    else:
                        min_msdmin[I] = min(min_msdmin[I], msd0[I])
                        max_msdmin[I] = max(max_msdmin[I], msd0[I])
                elif thr1 <= intens[j] <= thr2:
                    av_msdmid[I] += msd0[I]
                    if firstMid[I]:
                        min_msdmid[I] = msd0[I]
                        max_msdmid[I] = msd0[I]
                        firstMid[I] = False
                    else:
                        min_msdmid[I] = min(min_msdmid[I], msd0[I])
                        max_msdmid[I] = max(max_msdmid[I], msd0[I])
                else:
                    av_msdmax[I] += msd0[I]
                    if firstMax[I]:
                        min_msdmax[I] = msd0[I]
                        max_msdmax[I] = msd0[I]
                        firstMax[I] = False
                    else:
                        min_msdmax[I] = min(min_msdmax[I], msd0[I])
                        max_msdmax[I] = max(max_msdmax[I], msd0[I])
                if firstCell[I]:
                    min_msd[I] = msd0[I]
                    max_msd[I] = msd0[I]
                    firstCell[I] = False
                else:
                    min_msd[I] = min(min_msd[I], msd0[I])
                    max_msd[I] = max(max_msd[I], msd0[I])
            if intens[j] < thr1:
                n_min += 1.0
            elif thr1 <= intens[j] <= thr2:
                n_mid += 1.0
            else:
                n_max += 1.0

            p = np.polyfit(np.log(tau[1:mintrac]), np.log(msd0[1:mintrac]), 1)
            coeffi[j] = p[0]

        av_msd = av_msd/len(cells)
        av_msdmin = av_msdmin / n_min
        av_msdmid = av_msdmid / n_mid
        av_msdmax = av_msdmax / n_max

        av_p = np.polyfit(np.log(tau[1:mintrac]), np.log(av_msd[1:mintrac]), 1)
        av_pmin = np.polyfit(np.log(tau[1:mintrac]), np.log(av_msdmin[1:mintrac]), 1)
        av_pmid = np.polyfit(np.log(tau[1:mintrac]), np.log(av_msdmid[1:mintrac]), 1)
        av_pmax = np.polyfit(np.log(tau[1:mintrac]), np.log(av_msdmax[1:mintrac]), 1)

        # Create a dictionary with data for each row
        data = [av_p[0], av_pmin[0], av_pmid[0], av_pmax[0], folder.split('_')[0] + '(' + folder.split('_')[1] + ')']

        # Append the data to the DataFrame
        df.loc[-1] = data
        df.index = df.index + 1
        df = df.sort_index()

result_df = df.groupby('Condition').mean().reset_index()
result_df.set_index('Condition', inplace=True)
print(result_df)

plt.pcolor(result_df, cmap='magma')
plt.yticks(np.arange(0.5, len(result_df.index), 1), result_df.index)
plt.xticks(np.arange(0.5, len(result_df.columns), 1), result_df.columns)
plt.title('MSD coefficient')
plt.colorbar()
if printInFolder:
    df.to_csv(save_directory + 'msd_all_dataframe.csv')
    plt.savefig(save_directory + 'msd_heatmap.pdf')
else:
    plt.show()
