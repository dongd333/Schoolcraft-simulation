This folder had the RW3D model simulating the biodegradation column experiment by Witt et al. 1999.

RW3D_SC701 changes from RW3D_SC6 in declaring random seed through system clock function, so that the code can be run without re-compiling. 

The model simulates the complete reaction (forward, reverse, and transform) of Nitrate with mobile biomass first, then the complete reaction with immobile biomass, and last the reactions between CT and biomass.  It is the same as matlab code used in the 2015 paper. 

The simulation results are the same as presented in the paper. 

The model results are kept in the folder: C:\Users\Dong\Documents\MATLAB\Plot_RW3D.

Part_path_701_3 read in the particle_path outputs and plot them, then calculate the mean value and standard deviation. Because the number of CT particle is small, the curves are not smooth.

The smoothed curved are in the folder: Column_Plots for making the figures.