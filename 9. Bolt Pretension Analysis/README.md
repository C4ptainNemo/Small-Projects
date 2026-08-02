These are some of the files for automating analysis of bolt pretension in ANSYS. Not all files are included, as I've only included what is entirely my own work.

The Parametric Pretension Python script runs in ANSYS's script window. It automatically sets simulation steps, boundary conditions, bolt preload for each model. It then solves the model and saves the data to a csv. This data is post-processed in MATLAB an outputs the two images seen.

![combined_plots][Combined_Results_Plots.png]

![Bolt 1 Infinite Fatigue FOS Heatmap][Bolt 1 Infinite Fatigue FOS Heatmap.png]
