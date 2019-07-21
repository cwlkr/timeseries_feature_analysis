# Feature extraction and clustering
As some scripts are dependent on others, always keep them in the same folder.
Examples to use the function can be seen in an if(FALSE) block at the end of the scripts. They are not executed by sourcing the script.
## importing an experiment and extracting features.
### single experiment
For cleaning rawdata and extracting features of a single experiment the data.R script is used.
It's main function is load_data. The default values of the function load experiment 20190408_systII_siPOOLs_plate1_and_2_singlePulse_changed_order/20190408_150050_815 on the ubuntu OS. Default variables should be suited for the most recent experiment descriptions (column names).
As path are different between Operating systems, the first part of the path can be added to the mnt variable and the sub folder to the experiment variable. Alternatively, mnt can be set to "" and the experiment description can contain the full path to the experiment. 

The experiment path should point to the metadata file and the receptor data, plus the raw trajectory data should be located in a subfolder named "output-data".

The cleaned data, plus the extracted features are saved to a csv, in either the data folder in the current working directory or in the folder specified in data.folder. Saving can be deactivated with setting write.all = FALSE.

All Parameters are also explained in the source code.

### multiple experiments
The script load_multiple.R contains a wrapper function for load_data. instead of an experiment string it takes a vector of experiment path, resulting tables are concatenated to one data.frame for data and features, plus features normalized per experiment. Experiments loaded with this function are ready to be clustered.

## hierarchical clustering
Script name. Hierarchical clustering. Used Chapter of the Thesis.

Create a complete hierarchical clustering and plots cluster averages and cluster proportions.

the main function is:

create_clustering()

## One dimensional Clustering.
For one dimensional clustering on multiple experiments, features have to be extracted and normalized per experiment. Load_multiple does this automatically.
If the features are loaded from a csv. They have to be first normalized with normalize_frame(data, -exp.var, -site.var, -track-var) and then concatenated to a single data frame.
SiRNA treatmens with multiple CTRLS and WT, can be pooled by using pool_treatments. This discards the numbering for each CTRL +/- and WT. Additional treatments can be defined in the function.
Functions for 1D clustering are in onedclustering.R

A usual workflow that was used in my thesis is in the file example.R.

## plate comparisons
Plate comparisons are done on the script zscoretoctrl.R. The function experiment_comparison_plots() can be used to create the plots from the thesis, an example call can be found above the function in the comments. It works with all experiments loaded with the load_multiple function.

## feature importance
For feature Importance a random forest is trained in the script feature_importance.R. The workflow should work with all data loaded with the load_multiple or load_data function. Feature importance is not used anymore as we now look at single features.

## Average Plots.
Treatment average plots per group can be done in the plot.R script with the function plot_average_per_group() and are saved as a pdf to the the pdf.filename location.
Alternatively, the plots can be produced with the shiny app.