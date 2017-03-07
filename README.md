# stmm
This folder contains code supporting the spatiotemporal mixed model. This code was written for Matlab 2015a and may not run on earlier versions. It is recommended that the user start with "ExampleSimulations.m", which simulates data for a parcel from the Gordon parcellation (Gordon, E. M., Laumann, T. O., Adeyemo, B., Huckins, J. F., Kelley, W. M., & Petersen, S. E. (2016). Generation and evaluation of a cortical area parcellation from resting-state correlations. Cerebral cortex, 26(1), 288-303.) Note the design matrices are constructed from data from the HCP theory of mind experiment and are redistributed under the data use terms, see http://humanconnectome.org/data/data-use-terms/.


Code implementing the analysis for the HCP data is provided, including functions to extract cortical surface data from cifti files for input to matlab. Some degree of adaptation is required for use, and indicated in the "Begin user input--->" to "<------End user input" that appears at the beginning of each .m file.

Estimation is divided into four programs: 1) create subject design matrices; 2) extract cortical surface data from the HCP cifti files; 3) conduct first-level analysis using reduced bias AR parameters, conduct cross validation, and smooth temporal variance parameters; 4) fit the mixed model to the dhat and smoothed temporal variance parameters.

Questions or report bugs: Benjamin Risk, benjamin.risk@gmail.com
