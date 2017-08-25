
The base model was run in different computers. The inputs are the same. So only one base model is archieved. 

This  base model for the reactive transport runs in office computer, the files are in the folder: Ding_desktop\Documents\New Folder\Research\Final_Feb2017.

RW3D_Dec_27 is the RW3D code modified for the simulation. Random numbers were generated using clock.

To plot the figures, first, the output files from RW3D model runs -- rw3d_BTCs_0xx.dat are read in matlab code plot_initial_Nitrate2.m and plot_initial_CT_a3.m, and save the breakthrough data to individual csv file. Then, the individual data are normalized with the weighting factor (Factor.csv), which is calculated from transmissivity of the wells. This is because the bin size is larger than well radius. 

The final figure is plotted with matlab code in folder Reactive transport plots: plot_nitrate_smooth.m and plot_CT_smooth.m.