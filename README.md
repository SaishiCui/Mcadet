# Mcadet


1. The "Mcadet_function.R" file serves as the main R function file for our proposed framework. It contains detailed information about the required dependent packages and parameters used in the function.
2. The other files included in the framework have specific purposes. Some files are responsible for creating working datasets or performing benchmarking analyses.
3. The "PBMC_datasets_preparing.R" and "SPARSim_Simulation.R" files are used to generate 48 working datasets for real-world PBMC datasets and 48 simulated datasets, respectively.
4. The "Creating_HVGs_PBMC.R" file is utilized to create the ground truth of highly variable gene sets for the 48 PBMC datasets.
5. The "Creating_McadetHVGs_PBMC.R" and "Creating_McadetHVGs_Simulation" files are employed to generate the sets of highly variable genes using the Mcadet method on the PBMC and simulated datasets, respectively.
6. The remaining eight files, namely "Bre_FS.R", "Scry_FS.R", "M3Drop_FS.R", "NBDrop_FS.R", "SeuratDisp_FS.R", "SeuratVst_FS.R," "SeuratMvp_FS.R" and "Random_FS.R" represent other feature selection (FS) methods used to obtain highly variable gene sets.
7. The "Load_all_fs_results.R" file is responsible for loading all the required highly variable gene sets obtained by different FS methods.
8. The "UMAP_visualization_PBMC.R" and "UMAP_visualization_Simulation.R" files are used to generate UMAP visualization results.
9. The "Jaccard_results_PBMC.R" and "Jaccard_results_Simulation.R" files are employed to generate Jaccard similarity index comparison results.
10. The "Evaluation_clustering_PBMC.R" and "Evaluation_clustering_Simulation.R" files are responsible for generating clustering performance comparison results using five evaluation metrics.
11. The "Stability_analysis.R" and "Rarecells_analysis.R" files are used to generate results for consistency comparison and performance analysis on rare cell populations, respectively.
12. The "LengthTable.R" and "Line_graph_trend.R" files are employed to generate a table displaying the number of highly variable genes selected by different FS methods and a trend line graph for both PBMC and simulated datasets.
13. The "Expression_analysis.R" and "Running_Time.R" files are responsible for generating results related to expression analysis and running time comparison.

Each file serves a specific purpose within the framework and contributes to the overall analysis and evaluation of the Mcadet method.
