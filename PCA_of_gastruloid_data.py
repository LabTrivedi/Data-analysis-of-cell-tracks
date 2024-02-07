import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import json
import sys
import os
import csv
import pandas as pd
from sklearn.decomposition import PCA
from sklearn import preprocessing
import statistics
from pathlib import Path

def split_integral_in_half(points, function_values):
    # Create a sorted list of points and corresponding function values
    sorted_points, sorted_values = zip(*sorted(zip(points, function_values)))

    # Compute the integral of the function
    integral = integrate.cumtrapz(sorted_values, sorted_points, initial=0)

    # Find the point where the integral is split in half
    half_integral = integral[-1] / 2
    split_point = np.interp(half_integral, integral, sorted_points)

    return split_point
def length_of_midline(v):
    sum = 0.0
    for l in range(1, len(v)):
        sum += np.sqrt((v[l][0] - v[l - 1][0]) ** 2 + (v[l][1] - v[l - 1][1]) ** 2)
    return sum
def perform_PCA(data, saveName):
    #########################
    #
    # Perform PCA on the data
    #
    #########################
    # First center and scale the data
    scaled_data = preprocessing.scale(data)

    pca = PCA()  # create a PCA object
    pca.fit(scaled_data)  # do the math
    pca_data = pca.transform(scaled_data)  # get PCA coordinates for scaled_data

    #########################
    #
    # Draw a scree plot and a PCA plot
    #
    #########################

    # The following code constructs the Scree plot
    per_var = np.round(pca.explained_variance_ratio_ * 100, decimals=1)
    labels = ['PC' + str(x) for x in range(1, len(per_var) + 1)]

    plt.bar(x=range(1, len(per_var) + 1), height=per_var, tick_label=labels)
    plt.ylabel('Percentage of Explained Variance')
    plt.xlabel('Principal Component')
    plt.title('Screen Plot')
    if printInFolder:
        plt.savefig(save_directory + 'PCA_percentage_' + saveName + '_all.pdf')
        plt.cla()
        plt.clf()
    else:
        plt.show()

    # the following code makes a fancy looking plot using PC1 and PC2
    pca_df = pd.DataFrame(pca_data, index=list(map(str, data.index)), columns=labels)

    for sample in pca_df.index:
        condi, times, ynchi = sample.split('/')[0].split('_')
        times = int(times.split('h')[0])
        if ynchi == "CHI":
            plt.scatter(pca_df.at[sample, 'PC1'], pca_df.at[sample, 'PC2'], color=dict_condi[condi],
                        marker=dict_times[times - dict_subst[condi]])
        elif ynchi == "noCHI":
            plt.scatter(pca_df.at[sample, 'PC1'], pca_df.at[sample, 'PC2'], edgecolors=dict_condi[condi],
                        marker=dict_times[times - dict_subst[condi]], facecolors='none')
    plt.title('My PCA Graph')
    plt.xlabel('PC1 - {0}%'.format(per_var[0]))
    plt.ylabel('PC2 - {0}%'.format(per_var[1]))

    # for sample in pca_df.index:
    #     plt.annotate(sample, (pca_df.PC1.loc[sample], pca_df.PC2.loc[sample]))

    if printInFolder:
        plt.savefig(save_directory + 'PCA_graph_' + saveName + '_all.pdf')
        plt.cla()
        plt.clf()
    else:
        plt.show()

    #########################
    #
    # Determine which genes had the biggest influence on PC1
    #
    #########################

    ## get the name of the top 10 measurements (genes) that contribute
    ## most to pc1.
    ## first, get the loading scores
    loading_scores = pd.Series(pca.components_[0], index=data_columns)
    ## now sort the loading scores based on their magnitude
    sorted_loading_scores = loading_scores.abs().sort_values(ascending=False)

    # get the names of the top data contribution
    top_data = sorted_loading_scores[0:10].index.values

    ## print the gene names and their scores (and +/- sign)
    loading_scores[top_data].to_csv(save_directory + saveName + "_sorted_PCA_contribution.csv")
def check_files_exist(directory, file_list):
    """
    Check if all files in the given list exist in the specified directory.

    Parameters:
    - directory (str): The directory path where the files should be located.
    - file_list (list): A list of file names to check for existence.

    Returns:
    - True if all files exist, False otherwise.
    """

    for file_name in file_list:
        file_path = os.path.join(directory, file_name)

        if not os.path.exists(file_path):
            print(f"File not found: {file_name}")
            return False

    return True

# PARAMETERS TO MODIFY
directory = '//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/2D_polarization/important_data/'
save_directory = '//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/2D_polarization/important_data_plots/PCA/'
printInFolder = True
shiftTime = False

# CODE THAT SHOULD NOT BE MODIFIED
channel = "ch1_APprofile"

cond_folders = ["SL", "1i", "a2i"]
if shiftTime:
    dict_times = {-24: "P", 0: "o", 24: ",", 48: "^"}
    dict_subst = {'SL': 48, '1i': 72, 'a2i': 96}
    save_directory = save_directory + 'Shift_Time/'
else:
    dict_times = {48: "P", 72: "o", 96: ",", 120: "^"}
    dict_subst = {'SL': 0, '1i': 0, 'a2i': 0}
    save_directory = save_directory + 'No_Shift_Time/'
dict_condi = {"SL": "#838B8B", "1i": "#00688B", "a2i": "#8BBA8B"}

minTime = np.inf
times = []
data_columns = ['Area', 'Midline length', '<Intensity>', 'Polarization', 'Eccentricity']
total_folders = [folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory, folder)) and folder.endswith('CHI')]

data_tot = pd.DataFrame(columns=data_columns)
data_time0 = pd.DataFrame(columns=data_columns)
data_time1 = pd.DataFrame(columns=data_columns)
data_time2 = pd.DataFrame(columns=data_columns)
data_time3 = pd.DataFrame(columns=data_columns)
data_SL = pd.DataFrame(columns=data_columns)
data_1i = pd.DataFrame(columns=data_columns)
data_2i = pd.DataFrame(columns=data_columns)

file_list = ["data_total.csv", "data_SL.csv", "data_1i.csv", "data_2i.csv",
             "data_time0.csv", "data_time1.csv", "data_time2.csv", "data_time3.csv"]

if check_files_exist(save_directory, file_list):
    data_tot = pd.read_csv(save_directory + "data_total.csv", index_col=0)
    print('data_total.csv loaded')
    data_SL = pd.read_csv(save_directory + "data_SL.csv", index_col=0)
    print('data_SL.csv loaded')
    data_1i = pd.read_csv(save_directory + "data_SL.csv", index_col=0)
    print('data_1i.csv loaded')
    data_2i = pd.read_csv(save_directory + "data_2i.csv", index_col=0)
    print('data_2i.csv loaded')
    data_time0 = pd.read_csv(save_directory + "data_time0.csv", index_col=0)
    print('data_time0.csv loaded')
    data_time1 = pd.read_csv(save_directory + "data_time1.csv", index_col=0)
    print('data_time1.csv loaded')
    data_time2 = pd.read_csv(save_directory + "data_time2.csv", index_col=0)
    print('data_time2.csv loaded')
    data_time3 = pd.read_csv(save_directory + "data_time3.csv", index_col=0)
    print('data_time3.csv loaded')
else:
    for condition in cond_folders:
        # Max of polarization value (for relative polarization plots)
        maxp = 1.0
        if os.path.exists(directory + condition + "_maxp.csv"):
            with open(directory + condition + "_maxp.csv", mode='r') as csv_file:
                # Create a CSV reader object
                csv_reader = csv.reader(csv_file)

                # Read the single value from the CSV file
                for row in csv_reader:
                    if len(row) > 0:
                        maxp = float(row[0])
                        print("Read maxp: ", maxp)
                    else:
                        print("Empty CSV file")
        else:
            print(directory + condition + "maxp.csv not found")

        # Initial values of area
        area0 = 0.0
        if os.path.exists(directory + condition + "_area0.csv"):
            with open(directory + condition + "_area0.csv", mode='r') as csv_file:
                # Create a CSV reader object
                csv_reader = csv.reader(csv_file)

                # Read the single value from the CSV file
                for row in csv_reader:
                    if len(row) > 0:
                        area0 = float(row[0])
                        print("Read area0: ", area0)
                    else:
                        print("Empty CSV file")
        else:
            print(directory + condition + "area0.csv not found")

        # Initial values of midline length
        midline0 = 0.0
        if os.path.exists(directory + condition + "_midline0.csv"):
            with open(directory + condition + "_midline0.csv", mode='r') as csv_file:
                # Create a CSV reader object
                csv_reader = csv.reader(csv_file)

                # Read the single value from the CSV file
                for row in csv_reader:
                    if len(row) > 0:
                        midline0 = float(row[0])
                        print("Read midline0: ", midline0)
                    else:
                        print("Empty CSV file")
        else:
            print(directory + condition + "_midline0.csv not found")

        # Initial values of average instensity
        intns0 = 0.0
        if os.path.exists(directory + condition + "_intns0.csv"):
            with open(directory + condition + "_intns0.csv", mode='r') as csv_file:
                # Create a CSV reader object
                csv_reader = csv.reader(csv_file)

                # Read the single value from the CSV file
                for row in csv_reader:
                    if len(row) > 0:
                        intns0 = float(row[0])
                        print("Read intns0: ", intns0)
                    else:
                        print("Empty CSV file")
        else:
            print(directory + condition + "_intns0.csv not found")

        folders = [folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory, folder)) and condition in folder]
        print(folders)
        for folder in folders:
            _, time_folder, chi_folder = folder.split('_')
            print(condition, time_folder, chi_folder)
            path = directory + condition + "_" + time_folder + "_" + chi_folder + "/"

            # Create list of subfolders with .json data
            subfolders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]

            remove_fol = []
            for subfolder in subfolders:
                if not os.path.exists(path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json'):
                    remove_fol.append(subfolder)
            subfolders = [s for s in subfolders if (s not in remove_fol)]
            print(subfolders)
            AMIPE = np.zeros((len(subfolders), 5))
            for s, subfolder in enumerate(subfolders):
                print("Plots: " + path + subfolder)
                with open(path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json', 'r') as fcc_file:
                    with open(path + subfolder + "/result_segmentation/" + subfolder + '_morpho_straight_params.json', 'r') as fcc_file2:
                        with open(path + subfolder + "/result_segmentation/" + subfolder + '_morpho_params.json', 'r') as fcc_file3:
                            dataP = json.load(fcc_file)
                            dataA = json.load(fcc_file2)
                            dataM = json.load(fcc_file3)
                            lp = len(dataP)
                            la = len(dataA)
                            lm = len(dataA)
                            if lp != la or lp != lm or la != lm:
                                print('Inconsistent data!!!!')
                                sys.exit()
                            else:
                                AMIPE_sub = np.zeros((lp, 5))
                                for time in range(lp):
                                    lt = len(dataP[time][channel])
                                    points = np.linspace(0, 1, lt)
                                    function_values = dataP[time][channel]
                                    split_point = split_integral_in_half(points, function_values)
                                    m = length_of_midline(dataM[time]["midline"])
                                    AMIPE_sub[time, 0] = dataA[time]['area'] / area0
                                    AMIPE_sub[time, 1] = m / midline0
                                    AMIPE_sub[time, 2] = np.mean(function_values) / intns0
                                    AMIPE_sub[time, 3] = (max(split_point, 1.0 - split_point) - 0.5) / (maxp - 0.5)
                                    AMIPE_sub[time, 4] = dataA[time]['eccentricity']

                                # Add to data frame
                                # Total
                                data_tot.loc[condition + "_" + time_folder + "_" + chi_folder + "/" + subfolder, data_columns] = np.mean(AMIPE_sub, axis=0)
                                # By condition
                                if condition == "SL":
                                    data_SL.loc[condition + "_" + time_folder + "_" + chi_folder + "/" + subfolder, data_columns] = np.mean(AMIPE_sub, axis=0)
                                elif condition == "1i":
                                    data_1i.loc[condition + "_" + time_folder + "_" + chi_folder + "/" + subfolder, data_columns] = np.mean(AMIPE_sub, axis=0)
                                elif condition == "a2i":
                                    data_2i.loc[condition + "_" + time_folder + "_" + chi_folder + "/" + subfolder, data_columns] = np.mean(AMIPE_sub, axis=0)
                                # By time
                                t = int(time_folder.split('h')[0]) - dict_subst[condition]
                                if (shiftTime and t == -24) or (not shiftTime and t == 48):
                                    data_time0.loc[condition + "_" + time_folder + "_" + chi_folder + "/" + subfolder, data_columns] = np.mean(AMIPE_sub, axis=0)
                                elif (shiftTime and t == 0) or (not shiftTime and t == 72):
                                    data_time1.loc[condition + "_" + time_folder + "_" + chi_folder + "/" + subfolder, data_columns] = np.mean(AMIPE_sub, axis=0)
                                elif (shiftTime and t == 24) or (not shiftTime and t == 96):
                                    data_time2.loc[condition + "_" + time_folder + "_" + chi_folder + "/" + subfolder, data_columns] = np.mean(AMIPE_sub, axis=0)
                                elif (shiftTime and t == 48) or (not shiftTime and t == 120):
                                    data_time3.loc[condition + "_" + time_folder + "_" + chi_folder + "/" + subfolder, data_columns] = np.mean(AMIPE_sub, axis=0)

    data_tot.to_csv(save_directory + "data_total.csv")
    data_SL.to_csv(save_directory + "data_SL.csv")
    data_1i.to_csv(save_directory + "data_1i.csv")
    data_2i.to_csv(save_directory + "data_2i.csv")
    data_time0.to_csv(save_directory + "data_time0.csv")
    data_time1.to_csv(save_directory + "data_time1.csv")
    data_time2.to_csv(save_directory + "data_time2.csv")
    if not shiftTime:
        data_time3.to_csv(save_directory + "data_time3.csv")

print('Print data_tot')
perform_PCA(data_tot, 'data_tot')
print('Print data_SL')
perform_PCA(data_SL, 'data_SL')
print('Print data_1i')
perform_PCA(data_1i, 'data_1i')
print('Print data_2i')
perform_PCA(data_2i, 'data_2i')
print('Print data_time0')
perform_PCA(data_time0, 'data_time0')
print('Print data_time1')
perform_PCA(data_time1, 'data_time1')
print('Print data_time2')
perform_PCA(data_time2, 'data_time2')
if not shiftTime:
    print('Print data_time3')
    perform_PCA(data_time3, 'data_time3')


