# Data-analysis-of-cell-tracks
Python pipeline to perfom data analysis of cell tracks

Here, we provide several scripts to study the data related with cell tracks. The scripts are:
- Criteria_made_by_histogram_saddle_point_and_5_times_std: Automatized script to decide the thresholds for T- and T+ Bra/T expression.
- Heatmap_MSD: Heatmap of MSD coefficient depending on condition and time frame.
- Hexplot_data_plots: Set of hexbin plots where we print combinations of triples of data (x-axis, y-axis and color code) for the set of [time, area, midline length, average intensity, polartization, eccentricity]
- PCA_of_gastruloid_data: PCA analysis of [time, area, midline length, average intensity, polartization, eccentricity] data to understand the clusterings over times and conditions.
- Percentage_of_cell_div_respect_total_in_csv: CSV file where we have the ratio of cell that divide respect the total at each time frame.
- Pie_chart_categories_of_Bra_evolution: Quantification of cells according of Bra/T evolution categories.
- Plot_dynamics_arrow_at_origin: Plot of the difference between initial and final track position, divided by the number of tracks between them. Arrows are later centered at the origin. Plots are color coded according to Bra/T evolution categories.
- Plot_radial_histogram_Bra_split_angle_respect_percentage_of_spline_single: Radial histogram of angle between cells and gastuloid main directions at both tips.
- PLot_tracks_at_original_position: Plot of each cell track, color coded according to Bra/T evolution categories.
- Sankey_diagram_cell_Bra_evolution: Sankey diagram between initial and final expression type for cells in each condition and time frame.
- Table_of_xyz_domain_in_each_case: Plot to study the preference of directionality of cell movement (relX, relY, relZ) and to see the lenght of gastuloid over each dimension (reldiffX, reldiffY, reldiffZ). The objective of this script is to compare the direction of each cell with the dimensions of the gastruloid it lies within.

IMPORTANT: THE FIRST CODE YOU NEED TO USE IS Criteria_made_by_histogram_saddle_point_and_5_times_std, SINCE IT CREATES THE THESHOLD FILES THAT THE REST OF SCRIPTS USE FOR THEIR COMPUTATIONS. ALSO, YOU NEED TO PERFORM FIRST Hexplot_data_plots BEFORE PCA_of_gastruloid_data, SINCE THE LATTER USED DATA GENERATED IN THE FORMER.

Scripts are implemented with the mindset to automatize as much as possible, so the number of parameters that the used need to modify are few. The usual parameters to modify in each script are:
- directory: path where data from MOrgAna are stored
- save_directory: path where plots from scripts are saved. We recommend to have different folders for data and for plots
- mintrac: Minimum track length. This is for cases where it requires a minimum amount of time frames for a cell to consider its data significant (an example would be if we study the evolution of its expression but we only have two data frames). Each script has its own - - predetermined mintrac, but user is encouraged to try different cases to see which minimum track is appropiate for their data.
- printInFolder: boolean that saves the data plots at save_directory if True, if False, python will plot it directly to the screen. Generally, the user does not want to change this variables.
