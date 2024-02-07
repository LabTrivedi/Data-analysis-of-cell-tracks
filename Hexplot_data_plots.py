import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import integrate
import json
import sys
import os
import csv
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
        sum += np.sqrt((v[l][0] - v[l-1][0])**2 + (v[l][1] - v[l-1][1])**2)
    return sum
def log_AMI(TAMPIEnp):
    TAMIPEnp2 = TAMPIEnp.copy()
    TAMIPEnp2[:, 1] = np.log(TAMPIEnp[:, 1])
    TAMIPEnp2[:, 2] = np.log(TAMPIEnp[:, 2])
    TAMIPEnp2[:, 3] = np.log(TAMPIEnp[:, 3])
    return TAMIPEnp2

# PARAMETERS TO MODIFY
directory = '//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/2D_polarization/important_data/'
save_directory = '//trivedi.embl.es/trivedi/Kerim_Anlas/for_jordi/2D_polarization/important_data_plots/Hexbins/'
printInFolder = True

# CODE THAT SHOULD NOT BE MODIFIED
channel = "ch1_APprofile"
commonMin_list = [False, True]
globalRange_list = [False, True]
whatToPrint_list = ["all", "some"]
shiftTime_list = [False, True]

conditions = ['1i', 'SL', 'a2i']
chi_folders = ["CHI", "noCHI"]

minTime = np.inf
times = []

maxp = 0.0
area0 = 0.0
midline0 = 0.0
intns0 = 0.0

list_some = ['MET', 'MEI', 'MEP', 'TEM', 'TEI', 'TEP', 'TIM', 'TIP', 'TIE']
TAMIPE_names = ["Time (h)", "Area", "Midline length", "<Intensity>", "Polarization", "Eccentricity"]
TAMIPE_colors = ["hot", "winter", "copper", "plasma", "viridis", "cividis"]
TAMIPE_dict = {'T': 0, 'A': 1, 'M': 2, 'I': 3, 'P': 4, 'E': 5}

for globalRange in globalRange_list:
    for shiftTime in shiftTime_list:
        for commonMin in commonMin_list:
            for whatToPrint in whatToPrint_list:
                print('********************')
                print(' shiftTime = ' + str(shiftTime))
                print(' commonMin = ' + str(commonMin))
                print(' globalRange = ' + str(globalRange))
                print(' whatToPrint = ' + str(whatToPrint))
                print('********************')
                if globalRange:
                    globMin_list = [float('inf') for _ in range(3)]
                    globMax_list = [float('-inf') for _ in range(3)]
                    for condition in conditions:
                        print('Condition ' + condition)
                        time_folders = []
                        if condition == 'SL':
                            time_folders = ["48h", "72h"]
                        elif condition == '1i':
                            time_folders = ["72h", "96h"]
                        elif condition == 'a2i':
                            time_folders = ["96h", "120h"]

                        if not commonMin:
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
                                maxp = 0.0
                                for time_folder in time_folders:
                                    for chi_folder in chi_folders:
                                        path = directory + condition + "_" + time_folder + "_" + chi_folder + "/"
                                        # Create list of subfolders with .json data
                                        subfolders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
                                        remove_fol = []
                                        for subfolder in subfolders:
                                            if not os.path.exists(path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json'):
                                                remove_fol.append(subfolder)
                                        subfolders = [s for s in subfolders if (s not in remove_fol)]

                                        for subfolder in subfolders:
                                            print("Compute maxp: " + path + subfolder)
                                            with open(path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json', 'r') as fcc_file:
                                                data = json.load(fcc_file)
                                                ld = len(data)
                                                polarization = np.zeros(ld)
                                                for time in range(ld):
                                                    lt = len(data[time][channel])
                                                    points = np.linspace(0, 1, lt)
                                                    function_values = data[time][channel]
                                                    split_point = split_integral_in_half(points, function_values)
                                                    polarization[time] = max(split_point, 1.0 - split_point)

                                            maxp = max(maxp, max(polarization))
                                # Open the file in write mode
                                with open(directory + condition + "_maxp.csv", mode='w', newline='') as csv_file:
                                    # Create a CSV writer object
                                    csv_writer = csv.writer(csv_file)

                                    # Write the double value to the CSV file
                                    csv_writer.writerow([maxp])

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
                                for time_folder in time_folders:
                                    chi_folder = "CHI"
                                    path = directory + condition + "_" + time_folder + "_" + chi_folder + "/"
                                    # Create list of subfolders with .json data
                                    subfolders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
                                    remove_fol = []
                                    for subfolder in subfolders:
                                        if not os.path.exists(path + subfolder + "/result_segmentation/" + subfolder + '_morpho_straight_params.json'):
                                            remove_fol.append(subfolder)
                                    subfolders = [s for s in subfolders if (s not in remove_fol)]
                                    areas = np.zeros(len(subfolders))

                                    for s, subfolder in enumerate(subfolders):
                                        print("Compute area0: " + path + subfolder)
                                        with open(path + subfolder + "/result_segmentation/" + subfolder + '_morpho_straight_params.json', 'r') as fcc_file:
                                            data = json.load(fcc_file)
                                            areas[s] = data[0]['area']

                                    area0 = np.mean(areas)
                                # Open the file in write mode
                                with open(directory + condition + "_area0.csv", mode='w', newline='') as csv_file:
                                    # Create a CSV writer object
                                    csv_writer = csv.writer(csv_file)

                                    # Write the double value to the CSV file
                                    csv_writer.writerow([area0])
                            print(area0)

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
                                for time_folder in time_folders:
                                    chi_folder = "CHI"
                                    path = directory + condition + "_" + time_folder + "_" + chi_folder + "/"
                                    # Create list of subfolders with .json data
                                    subfolders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
                                    remove_fol = []
                                    for subfolder in subfolders:
                                        if not os.path.exists(
                                                path + subfolder + "/result_segmentation/" + subfolder + '_morpho_params.json'):
                                            remove_fol.append(subfolder)
                                    subfolders = [s for s in subfolders if (s not in remove_fol)]
                                    midlines = np.zeros(len(subfolders))

                                    for s, subfolder in enumerate(subfolders):
                                        print("Compute midline0: " + path + subfolder)
                                        with open(path + subfolder + "/result_segmentation/" + subfolder + '_morpho_params.json',
                                                  'r') as fcc_file:
                                            data = json.load(fcc_file)
                                            midlines[s] = length_of_midline(data[0]["midline"])

                                    midline0 = np.mean(midlines)
                                # Open the file in write mode
                                with open(directory + condition + "_midline0.csv", mode='w', newline='') as csv_file:
                                    # Create a CSV writer object
                                    csv_writer = csv.writer(csv_file)

                                    # Write the double value to the CSV file
                                    csv_writer.writerow([midline0])

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
                                for time_folder in time_folders:
                                    chi_folder = "noCHI"
                                    path = directory + condition + "_" + time_folder + "_" + chi_folder + "/"
                                    # Create list of subfolders with .json data
                                    subfolders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
                                    remove_fol = []
                                    for subfolder in subfolders:
                                        if not os.path.exists(
                                                path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json'):
                                            remove_fol.append(subfolder)
                                    subfolders = [s for s in subfolders if (s not in remove_fol)]
                                    intensities = np.zeros(len(subfolders))

                                    for s, subfolder in enumerate(subfolders):
                                        print("Compute intns0: " + path + subfolder)
                                        with open(path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json',
                                                  'r') as fcc_file:
                                            data = json.load(fcc_file)
                                            intensities[s] = np.mean(data[0][channel])

                                    intns0 = np.mean(intensities)
                                # Open the file in write mode
                                with open(directory + condition + "_intns0.csv", mode='w', newline='') as csv_file:
                                    # Create a CSV writer object
                                    csv_writer = csv.writer(csv_file)

                                    # Write the double value to the CSV file
                                    csv_writer.writerow([intns0])

                        # Perform plots
                        for chi_folder in chi_folders:
                            tamipe_name = condition + "_" + chi_folder
                            if commonMin:
                                tamipe_name = tamipe_name + "_common"
                            if os.path.exists(directory + tamipe_name + ".csv"):
                                TAMIPEnp = np.loadtxt(directory + tamipe_name + ".csv", delimiter=',')
                                print(directory + tamipe_name + ".csv has been loaded")
                            else:
                                EAPIMT = []
                                firstTimeFolder = True
                                for time_folder in time_folders:
                                    path = directory + condition + "_" + time_folder + "_" + chi_folder + "/"
                                    # Create list of subfolders with .json data
                                    subfolders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
                                    remove_fol = []
                                    for subfolder in subfolders:
                                        if not os.path.exists(path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json'):
                                            remove_fol.append(subfolder)
                                    subfolders = [s for s in subfolders if (s not in remove_fol)]

                                    iniTime = int(time_folder.split('h')[0])
                                    if iniTime not in times:
                                        times.append(iniTime)
                                    for subfolder in subfolders:
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
                                                        for time in range(lp):
                                                            lt = len(dataP[time][channel])
                                                            points = np.linspace(0, 1, lt)
                                                            function_values = dataP[time][channel]
                                                            split_point = split_integral_in_half(points, function_values)
                                                            m = length_of_midline(dataM[time]["midline"])
                                                            if shiftTime:
                                                                t = float(time_folder[0:-1]) - float(time_folders[0][0:-1]) + time*10.0/60.0
                                                            else:
                                                                t = float(time_folder[0:-1]) + time * 10.0 / 60.0
                                                            EAPIMT.append([dataA[time]['eccentricity'],
                                                                           dataA[time]['area']/area0,
                                                                           (max(split_point, 1.0 - split_point) - 0.5) / (maxp - 0.5),
                                                                           np.mean(function_values)/intns0,
                                                                           m/midline0,
                                                                           t])

                                EAPIMTnp = np.array(EAPIMT)
                                TAMIPEnp = EAPIMTnp.copy()
                                TAMIPEnp[:, 0] = EAPIMTnp[:, 5]
                                TAMIPEnp[:, 2] = EAPIMTnp[:, 4]
                                TAMIPEnp[:, 4] = EAPIMTnp[:, 2]
                                TAMIPEnp[:, 5] = EAPIMTnp[:, 0]
                                np.savetxt(directory + tamipe_name + ".csv", TAMIPEnp, delimiter=',')

                            for i in range(3):
                                globMin_list[i] = min(np.min(TAMIPEnp[:, i + 1]), globMin_list[i])
                                globMax_list[i] = max(np.max(TAMIPEnp[:, i + 1]), globMax_list[i])

                    print(globMin_list)
                    print(globMax_list)

                if commonMin:
                    for c, condition in enumerate(conditions):
                        print('-----------------------')
                        print(condition)
                        time_folders = []
                        if condition == 'SL':
                            time_folders = ["48h", "72h", "96h"]
                        elif condition == '1i':
                            time_folders = ["72h", "96h", "120h"]
                        elif condition == 'a2i':
                            time_folders = ["72h", "96h", "120h"]

                        # Max of polarization value (for relative polarization plots)
                        if os.path.exists(directory + condition + "_maxp.csv"):
                            with open(directory + condition + "_maxp.csv", mode='r') as csv_file:
                                # Create a CSV reader object
                                csv_reader = csv.reader(csv_file)

                                # Read the single value from the CSV file
                                for row in csv_reader:
                                    if len(row) > 0:
                                        print(row[0])
                                        if c == 0:
                                            maxp = float(row[0])
                                        else:
                                            maxp = max(maxp, float(row[0]))
                                        print("Read maxp: ", maxp)
                                    else:
                                        print("Empty CSV file")
                        else:
                            print(directory + condition + '_maxp.csv not found. Please run first with commonMin as True')
                            sys.exit()

                        # Initial values of area
                        if os.path.exists(directory + condition + "_area0.csv"):
                            with open(directory + condition + "_area0.csv", mode='r') as csv_file:
                                # Create a CSV reader object
                                csv_reader = csv.reader(csv_file)

                                # Read the single value from the CSV file
                                for row in csv_reader:
                                    if len(row) > 0:
                                        print(row[0])
                                        if c == 0:
                                            area0 = float(row[0])
                                        else:
                                            area0 = min(area0, float(row[0]))
                                        print("Read area0: ", area0)
                                    else:
                                        print("Empty CSV file")
                        else:
                            print(directory + condition + '_area0.csv not found. Please run first with commonMin as False')
                            sys.exit()

                        # Initial values of midline length
                        if os.path.exists(directory + condition + "_midline0.csv"):
                            with open(directory + condition + "_midline0.csv", mode='r') as csv_file:
                                # Create a CSV reader object
                                csv_reader = csv.reader(csv_file)

                                # Read the single value from the CSV file
                                for row in csv_reader:
                                    if len(row) > 0:
                                        print(row[0])
                                        if c == 0:
                                            midline0 = float(row[0])
                                        else:
                                            midline0 = min(midline0, float(row[0]))
                                        print("Read midline0: ", midline0)
                                    else:
                                        print("Empty CSV file")
                        else:
                            print(directory + condition + '_midline0.csv not found. Please run first with commonMin as True')
                            sys.exit()

                        # Initial values of average instensity
                        if os.path.exists(directory + condition + "_intns0.csv"):
                            with open(directory + condition + "_intns0.csv", mode='r') as csv_file:
                                # Create a CSV reader object
                                csv_reader = csv.reader(csv_file)

                                # Read the single value from the CSV file
                                for row in csv_reader:
                                    if len(row) > 0:
                                        print(row[0])
                                        if c == 0:
                                            intns0 = float(row[0])
                                        else:
                                            intns0 = min(intns0, float(row[0]))
                                        print("Read intns0: ", intns0)
                                    else:
                                        print("Empty CSV file")
                        else:
                            print(directory + condition + '_intns0.csv not found. Please run first with commonMin as True')
                            sys.exit()
                    print('-----------------------')

                for condition in conditions:
                    print('Condition ' + condition)
                    time_folders = []
                    if condition == 'SL':
                        time_folders = ["48h", "72h", "96h"]
                    elif condition == '1i':
                        time_folders = ["72h", "96h", "120h"]
                    elif condition == 'a2i':
                        time_folders = ["72h", "96h", "120h"]

                    if not commonMin:
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
                            maxp = 0.0
                            for time_folder in time_folders:
                                for chi_folder in chi_folders:
                                    path = directory + condition + "_" + time_folder + "_" + chi_folder + "/"
                                    # Create list of subfolders with .json data
                                    subfolders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
                                    remove_fol = []
                                    for subfolder in subfolders:
                                        if not os.path.exists(path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json'):
                                            remove_fol.append(subfolder)
                                    subfolders = [s for s in subfolders if (s not in remove_fol)]

                                    for subfolder in subfolders:
                                        print("Compute maxp: " + path + subfolder)
                                        with open(path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json', 'r') as fcc_file:
                                            data = json.load(fcc_file)
                                            ld = len(data)
                                            polarization = np.zeros(ld)
                                            for time in range(ld):
                                                lt = len(data[time][channel])
                                                points = np.linspace(0, 1, lt)
                                                function_values = data[time][channel]
                                                split_point = split_integral_in_half(points, function_values)
                                                polarization[time] = max(split_point, 1.0 - split_point)

                                        maxp = max(maxp, max(polarization))
                            # Open the file in write mode
                            with open(directory + condition + "_maxp.csv", mode='w', newline='') as csv_file:
                                # Create a CSV writer object
                                csv_writer = csv.writer(csv_file)

                                # Write the double value to the CSV file
                                csv_writer.writerow([maxp])

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
                            for time_folder in time_folders:
                                chi_folder = "CHI"
                                path = directory + condition + "_" + time_folder + "_" + chi_folder + "/"
                                # Create list of subfolders with .json data
                                subfolders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
                                remove_fol = []
                                for subfolder in subfolders:
                                    if not os.path.exists(path + subfolder + "/result_segmentation/" + subfolder + '_morpho_straight_params.json'):
                                        remove_fol.append(subfolder)
                                subfolders = [s for s in subfolders if (s not in remove_fol)]
                                areas = np.zeros(len(subfolders))

                                for s, subfolder in enumerate(subfolders):
                                    print("Compute area0: " + path + subfolder)
                                    with open(path + subfolder + "/result_segmentation/" + subfolder + '_morpho_straight_params.json', 'r') as fcc_file:
                                        data = json.load(fcc_file)
                                        areas[s] = data[0]['area']

                                area0 = np.mean(areas)
                            # Open the file in write mode
                            with open(directory + condition + "_area0.csv", mode='w', newline='') as csv_file:
                                # Create a CSV writer object
                                csv_writer = csv.writer(csv_file)

                                # Write the double value to the CSV file
                                csv_writer.writerow([area0])

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
                            for time_folder in time_folders:
                                chi_folder = "CHI"
                                path = directory + condition + "_" + time_folder + "_" + chi_folder + "/"
                                # Create list of subfolders with .json data
                                subfolders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
                                remove_fol = []
                                for subfolder in subfolders:
                                    if not os.path.exists(
                                            path + subfolder + "/result_segmentation/" + subfolder + '_morpho_params.json'):
                                        remove_fol.append(subfolder)
                                subfolders = [s for s in subfolders if (s not in remove_fol)]
                                midlines = np.zeros(len(subfolders))

                                for s, subfolder in enumerate(subfolders):
                                    print("Compute midline0: " + path + subfolder)
                                    with open(path + subfolder + "/result_segmentation/" + subfolder + '_morpho_params.json',
                                              'r') as fcc_file:
                                        data = json.load(fcc_file)
                                        midlines[s] = length_of_midline(data[0]["midline"])

                                midline0 = np.mean(midlines)
                            # Open the file in write mode
                            with open(directory + condition + "_midline0.csv", mode='w', newline='') as csv_file:
                                # Create a CSV writer object
                                csv_writer = csv.writer(csv_file)

                                # Write the double value to the CSV file
                                csv_writer.writerow([midline0])

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
                            for time_folder in time_folders:
                                chi_folder = "noCHI"
                                path = directory + condition + "_" + time_folder + "_" + chi_folder + "/"
                                # Create list of subfolders with .json data
                                subfolders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
                                remove_fol = []
                                for subfolder in subfolders:
                                    if not os.path.exists(
                                            path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json'):
                                        remove_fol.append(subfolder)
                                subfolders = [s for s in subfolders if (s not in remove_fol)]
                                intensities = np.zeros(len(subfolders))

                                for s, subfolder in enumerate(subfolders):
                                    print("Compute intns0: " + path + subfolder)
                                    with open(path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json',
                                              'r') as fcc_file:
                                        data = json.load(fcc_file)
                                        intensities[s] = np.mean(data[0][channel])

                                intns0 = np.mean(intensities)
                            # Open the file in write mode
                            with open(directory + condition + "_intns0.csv", mode='w', newline='') as csv_file:
                                # Create a CSV writer object
                                csv_writer = csv.writer(csv_file)

                                # Write the double value to the CSV file
                                csv_writer.writerow([intns0])

                    # Perform plots
                    for chi_folder in chi_folders:
                        tamipe_name = condition + "_" + chi_folder
                        if commonMin:
                            tamipe_name = tamipe_name + "_common"
                        if shiftTime:
                            tamipe_name = tamipe_name + "_time_shifted"
                        if os.path.exists(directory + tamipe_name + ".csv"):
                            TAMIPEnp = np.loadtxt(directory + tamipe_name + ".csv", delimiter=',')
                            print(directory + tamipe_name + ".csv has been loaded")
                        else:
                            EAPIMT = []
                            firstTimeFolder = True
                            for time_folder in time_folders:
                                path = directory + condition + "_" + time_folder + "_" + chi_folder + "/"
                                # Create list of subfolders with .json data
                                subfolders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
                                remove_fol = []
                                for subfolder in subfolders:
                                    if not os.path.exists(path + subfolder + "/result_segmentation/" + subfolder + '_fluo_intensity.json'):
                                        remove_fol.append(subfolder)
                                subfolders = [s for s in subfolders if (s not in remove_fol)]

                                iniTime = int(time_folder.split('h')[0])
                                if iniTime not in times:
                                    times.append(iniTime)
                                for subfolder in subfolders:
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
                                                    for time in range(lp):
                                                        lt = len(dataP[time][channel])
                                                        points = np.linspace(0, 1, lt)
                                                        function_values = dataP[time][channel]
                                                        split_point = split_integral_in_half(points, function_values)
                                                        m = length_of_midline(dataM[time]["midline"])
                                                        if shiftTime:
                                                            t = float(time_folder[0:-1]) - float(
                                                                time_folders[0][0:-1]) + time * 10.0 / 60.0
                                                        else:
                                                            t = float(time_folder[0:-1]) + time * 10.0 / 60.0
                                                        EAPIMT.append([dataA[time]['eccentricity'],
                                                                       dataA[time]['area']/area0,
                                                                       (max(split_point, 1.0 - split_point) - 0.5) / (maxp - 0.5),
                                                                       np.mean(function_values)/intns0,
                                                                       m/midline0,
                                                                       t])

                            EAPIMTnp = np.array(EAPIMT)
                            TAMIPEnp = EAPIMTnp.copy()
                            TAMIPEnp[:, 0] = EAPIMTnp[:, 5]
                            TAMIPEnp[:, 2] = EAPIMTnp[:, 4]
                            TAMIPEnp[:, 4] = EAPIMTnp[:, 2]
                            TAMIPEnp[:, 5] = EAPIMTnp[:, 0]
                            np.savetxt(directory + tamipe_name + ".csv", TAMIPEnp, delimiter=',')

                        TAMIPEnp_lin = TAMIPEnp.copy()
                        TAMIPEnp_log = log_AMI(TAMIPEnp)

                        time_counts = [int(x[:-1]) for x in time_folders]
                        if shiftTime:
                            tmin = 0
                            tmax = max(time_counts) - min(time_counts) + 24
                        else:
                            tmin = min(time_counts)
                            tmax = max(time_counts) + 24

                        if whatToPrint == "some":
                            # nx = int(math.sqrt(len(list_some)))
                            # ny = math.ceil(len(list_some)/nx)
                            nx = 3
                            ny = 3

                            # Plot all formats
                            fig_lin, ax_lin = plt.subplots(nx, ny, figsize=(15, 10))
                            plt.subplots_adjust(wspace=0.3, hspace=0.4)

                            fig_log, ax_log = plt.subplots(nx, ny, figsize=(15, 10))
                            plt.subplots_adjust(wspace=0.3, hspace=0.4)

                            xx = 0
                            yy = 0
                            for ll in range(len(list_some)):
                                i = TAMIPE_dict[list_some[ll][0]]
                                j = TAMIPE_dict[list_some[ll][1]]
                                c = TAMIPE_dict[list_some[ll][2]]

                                ax_lin[xx, yy].set_xlabel(TAMIPE_names[i])
                                ax_lin[xx, yy].set_ylabel(TAMIPE_names[j])
                                ax_lin[xx, yy].set_title(TAMIPE_names[c])
                                ax_log[xx, yy].set_xlabel(TAMIPE_names[i])
                                ax_log[xx, yy].set_ylabel(TAMIPE_names[j])
                                ax_log[xx, yy].set_title(TAMIPE_names[c])
                                if i == 0:
                                    ax_lin[xx, yy].set_xlim([tmin, tmax])
                                    ax_log[xx, yy].set_xlim([tmin, tmax])
                                elif i in [4, 5]:
                                    ax_lin[xx, yy].set_xlim([0.0, 1.0])
                                    ax_log[xx, yy].set_xlim([0.0, 1.0])
                                else:
                                    imin_lin = np.min(TAMIPEnp_lin[:, i]) if not globalRange else globMin_list[i - 1]
                                    imax_lin = np.max(TAMIPEnp_lin[:, i]) if not globalRange else globMax_list[i - 1]
                                    imin_log = np.min(TAMIPEnp_log[:, i]) if not globalRange else np.log(globMin_list[i - 1])
                                    imax_log = np.max(TAMIPEnp_log[:, i]) if not globalRange else np.log(globMax_list[i - 1])

                                    ax_lin[xx, yy].set_xlim([imin_lin, imax_lin])
                                    ax_log[xx, yy].set_xlim([imin_log, imax_log])
                                if j == 0:
                                    ax_lin[xx, yy].set_ylim([tmin, tmax])
                                    ax_log[xx, yy].set_ylim([tmin, tmax])
                                elif j in [4, 5]:
                                    ax_lin[xx, yy].set_ylim([0.0, 1.0])
                                    ax_log[xx, yy].set_ylim([0.0, 1.0])
                                else:
                                    jmin_lin = np.min(TAMIPEnp_lin[:, j]) if not globalRange else globMin_list[j - 1]
                                    jmax_lin = np.max(TAMIPEnp_lin[:, j]) if not globalRange else globMax_list[j - 1]
                                    jmin_log = np.min(TAMIPEnp_log[:, j]) if not globalRange else np.log(globMin_list[j - 1])
                                    jmax_log = np.max(TAMIPEnp_log[:, j]) if not globalRange else np.log(globMax_list[j - 1])

                                    ax_lin[xx, yy].set_ylim([jmin_lin, jmax_lin])
                                    ax_log[xx, yy].set_ylim([jmin_log, jmax_log])

                                if c == 0:
                                    im_lin = ax_lin[xx, yy].hexbin(x=TAMIPEnp_lin[:, i], y=TAMIPEnp_lin[:, j],
                                                                 C=TAMIPEnp_lin[:, c], gridsize=50, vmin=tmin,
                                                                 vmax=tmax, cmap=TAMIPE_colors[c])
                                    im_log = ax_log[xx, yy].hexbin(x=TAMIPEnp_log[:, i], y=TAMIPEnp_log[:, j],
                                                                 C=TAMIPEnp_log[:, c], gridsize=50, vmin=tmin,
                                                                 vmax=tmax, cmap=TAMIPE_colors[c])
                                elif c in [4, 5]:
                                    im_lin = ax_lin[xx, yy].hexbin(x=TAMIPEnp_lin[:, i], y=TAMIPEnp_lin[:, j],
                                                                 C=TAMIPEnp_lin[:, c], gridsize=50, vmin=0.0, vmax=1.0,
                                                                 cmap=TAMIPE_colors[c])
                                    im_log = ax_log[xx, yy].hexbin(x=TAMIPEnp_log[:, i], y=TAMIPEnp_log[:, j],
                                                                 C=TAMIPEnp_log[:, c], gridsize=50, vmin=0.0, vmax=1.0,
                                                                 cmap=TAMIPE_colors[c])
                                else:
                                    cmin_lin = np.min(TAMIPEnp_lin[:, c]) if not globalRange else globMin_list[c - 1]
                                    cmax_lin = np.max(TAMIPEnp_lin[:, c]) if not globalRange else globMax_list[c - 1]
                                    cmin_log = np.min(TAMIPEnp_log[:, c]) if not globalRange else np.log(globMin_list[c - 1])
                                    cmax_log = np.max(TAMIPEnp_log[:, c]) if not globalRange else np.log(globMax_list[c - 1])

                                    im_lin = ax_lin[xx, yy].hexbin(x=TAMIPEnp_lin[:, i], y=TAMIPEnp_lin[:, j],
                                                                 C=TAMIPEnp_lin[:, c], gridsize=50, vmin=cmin_lin,
                                                                 vmax=cmax_lin, cmap=TAMIPE_colors[c])
                                    im_log = ax_log[xx, yy].hexbin(x=TAMIPEnp_log[:, i], y=TAMIPEnp_log[:, j],
                                                                 C=TAMIPEnp_log[:, c], gridsize=50, vmin=cmin_log,
                                                                 vmax=cmax_log, cmap=TAMIPE_colors[c])
                                plt.colorbar(im_lin, ax=ax_lin[xx, yy])
                                plt.colorbar(im_log, ax=ax_log[xx, yy])

                                xx += 1
                                if xx == nx:
                                    xx = 0
                                    yy += 1
                        elif whatToPrint == "all":
                            ldata = len(TAMIPE_names)

                            fig_lin, ax_lin = plt.subplots(ldata, int((ldata - 1) * (ldata- 2 ) / 2), figsize=(40, 20))
                            plt.subplots_adjust(wspace=0.5, hspace=0.5)

                            fig_log, ax_log = plt.subplots(ldata, int((ldata - 1) * (ldata - 2) / 2), figsize=(40, 20))
                            plt.subplots_adjust(wspace=0.5, hspace=0.5)
                            for c in range(ldata):
                                k = 0
                                for i in range(ldata):
                                    if i != c:
                                        for j in range(i+1, ldata):
                                            if j != c:
                                                ax_lin[c, k].set_xlabel(TAMIPE_names[i])
                                                ax_lin[c, k].set_ylabel(TAMIPE_names[j])
                                                ax_lin[c, k].set_title(TAMIPE_names[c])
                                                ax_log[c, k].set_xlabel(TAMIPE_names[i])
                                                ax_log[c, k].set_ylabel(TAMIPE_names[j])
                                                ax_log[c, k].set_title(TAMIPE_names[c])
                                                if i == 0:
                                                    ax_lin[c, k].set_xlim([tmin, tmax])
                                                    ax_log[c, k].set_xlim([tmin, tmax])
                                                elif i in [4, 5]:
                                                    ax_lin[c, k].set_xlim([0.0, 1.0])
                                                    ax_log[c, k].set_xlim([0.0, 1.0])
                                                else:
                                                    imin_lin = np.min(TAMIPEnp_lin[:, i]) if not globalRange else globMin_list[i - 1]
                                                    imax_lin = np.max(TAMIPEnp_lin[:, i]) if not globalRange else globMax_list[i - 1]
                                                    imin_log = np.min(TAMIPEnp_log[:, i]) if not globalRange else np.log(globMin_list[i - 1])
                                                    imax_log = np.max(TAMIPEnp_log[:, i]) if not globalRange else np.log(globMax_list[i - 1])
                                                    ax_lin[c, k].set_xlim([imin_lin, imax_lin])
                                                    ax_log[c, k].set_xlim([imin_log, imax_log])
                                                if j == 0:
                                                    ax_lin[c, k].set_ylim([tmin, tmax])
                                                    ax_log[c, k].set_ylim([tmin, tmax])
                                                elif j in [4, 5]:
                                                    ax_lin[c, k].set_ylim([0.0, 1.0])
                                                    ax_log[c, k].set_ylim([0.0, 1.0])
                                                else:
                                                    jmin_lin = np.min(TAMIPEnp_lin[:, j]) if not globalRange else globMin_list[j - 1]
                                                    jmax_lin = np.max(TAMIPEnp_lin[:, j]) if not globalRange else globMax_list[j - 1]
                                                    jmin_log = np.min(TAMIPEnp_log[:, j]) if not globalRange else np.log(globMin_list[j - 1])
                                                    jmax_log = np.max(TAMIPEnp_log[:, j]) if not globalRange else np.log(globMax_list[j - 1])
                                                    ax_lin[c, k].set_ylim([jmin_lin, jmax_lin])
                                                    ax_log[c, k].set_ylim([jmin_log, jmax_log])

                                                if c == 0:
                                                    im_lin = ax_lin[c, k].hexbin(x=TAMIPEnp_lin[:, i], y=TAMIPEnp_lin[:, j],
                                                                                 C=TAMIPEnp_lin[:, c], gridsize=50, vmin=tmin,
                                                                                 vmax=tmax, cmap=TAMIPE_colors[c])
                                                    im_log = ax_log[c, k].hexbin(x=TAMIPEnp_log[:, i], y=TAMIPEnp_log[:, j],
                                                                                 C=TAMIPEnp_log[:, c], gridsize=50, vmin=tmin,
                                                                                 vmax=tmax, cmap=TAMIPE_colors[c])
                                                elif c in [4, 5]:
                                                    im_lin = ax_lin[c, k].hexbin(x=TAMIPEnp_lin[:, i], y=TAMIPEnp_lin[:, j],
                                                                                 C=TAMIPEnp_lin[:, c], gridsize=50, vmin=0.0, vmax=1.0,
                                                                                 cmap=TAMIPE_colors[c])
                                                    im_log = ax_log[c, k].hexbin(x=TAMIPEnp_log[:, i], y=TAMIPEnp_log[:, j],
                                                                                 C=TAMIPEnp_log[:, c], gridsize=50, vmin=0.0, vmax=1.0,
                                                                                 cmap=TAMIPE_colors[c])
                                                else:
                                                    cmin_lin = np.min(TAMIPEnp_lin[:, c]) if not globalRange else globMin_list[c - 1]
                                                    cmax_lin = np.max(TAMIPEnp_lin[:, c]) if not globalRange else globMax_list[c - 1]
                                                    cmin_log = np.min(TAMIPEnp_log[:, c]) if not globalRange else np.log(globMin_list[c - 1])
                                                    cmax_log = np.max(TAMIPEnp_log[:, c]) if not globalRange else np.log(globMax_list[c - 1])
                                                    im_lin = ax_lin[c, k].hexbin(x=TAMIPEnp_lin[:, i], y=TAMIPEnp_lin[:, j],
                                                                                 C=TAMIPEnp_lin[:, c], gridsize=50, vmin=cmin_lin,
                                                                                 vmax=cmax_lin, cmap=TAMIPE_colors[c])
                                                    im_log = ax_log[c, k].hexbin(x=TAMIPEnp_log[:, i], y=TAMIPEnp_log[:, j],
                                                                                 C=TAMIPEnp_log[:, c], gridsize=50, vmin=cmin_log,
                                                                                 vmax=cmax_log, cmap=TAMIPE_colors[c])
                                                plt.colorbar(im_lin, ax=ax_lin[c, k])
                                                plt.colorbar(im_log, ax=ax_log[c, k])
                                                k += 1

                        if printInFolder:
                            file_name = condition + "_" + chi_folder
                            name = whatToPrint
                            if commonMin:
                                name = name + '_common'
                            if globalRange:
                                name = name + '_globalRange'
                            save_folder = 'No_Shift_Time/'
                            if shiftTime:
                                save_folder = 'Shift_Time/'

                            fig_lin.savefig(save_directory + save_folder + file_name + '_' + name + '_linear.pdf')
                            fig_log.savefig(save_directory + save_folder + file_name + '_' + name + '_logarithmic.pdf')
                        else:
                            plt.show()
                        plt.close(fig=fig_lin)
                        plt.close(fig=fig_log)