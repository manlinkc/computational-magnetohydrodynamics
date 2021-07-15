# computational-magnetohydrodynamics
Code submission for a short research project on the topic 'Computational Magnetohydrodynamics' as part of the MPhil in Scientific Computing, University of Cambridge 2020-2021 program.

### 1D MHD Code 

This directory contains all of the code used for the 1D MHD simulations. 

**1dmhdfunctions.C** – This file contains all of the functions used in the calculations within the 1D MHD solver. This file contains functions relating to the conversion between variables, time step calculations and flux calculations. 

**1dmhdfunctions.h** - This file contains function headers for the file 1dmhdfunctions.C 

**1dmhdmain.C** – This is the main program used for the 1DMHD solver. 

**1dmhdtests.C** – This file contains functions for setting up the initial conditions in the computational domain for each test simulated.  

**1dmhdtests.h** – This file contains function headers for the file 1dmhdtests.C

**inputs.txt** - Textfile containing the settings parameters for the test to be simulated in the 1D MHD solver. 

**master.txt** - Textfile containing the settings parameters for all of the 1D validation tests considered in this written assignment. The settings for a particular test can be copied from the master.txt file to the inputs.txt file when needed. 

**GNUmakefile** – Builds an executable named ‘my1dmhd’ which can be compiled by typing make in the directory.   

### 2D MHD Code 

This directory contains all of the code used for the 2D MHD simulations. 

**2dmhdfunctions.C** – This file contains all of the functions used in the calculations within the 2D MHD solver. This file also contains functions for setting up the initial conditions in the computational domain for each test simulated.  This file contains functions relating to the conversion between variables, time step calculations and flux calculations. 

**2dmhdfunctions.h**  - This file contains function headers for the file 2dmhdfunctions.C 

**2dmhdmain.C** – This is the main program used for the 2DMHD solver. This file also computes the L1 error at regular intervals when no correction is used and outputs the values to a datafile.  

**inputs.txt** - Textfile containing the settings parameters for the test to be simulated in the 2D MHD solver. 

**master.txt** - Textfile containing the settings parameters for all of the 2D validation tests considered in this written assignment. The settings for a particular test can be copied from the master.txt file to the inputs.txt file when needed. 

**GNUmakefile** – Builds an executable named ‘my2dmhd’ which can be compiled by typing make in the directory.   

### Divergence Cleaning Code 

This directory contains all of the code used for the Divergence Cleaning 2D MHD simulations. 

**dcfunctions.C** – This file contains all of the functions used in the calculations within the Divergence Cleaning 2D MHD solver. This file also contains functions for setting up the initial conditions in the computational domain for each test simulated. This file contains functions relating to the conversion between variables, time step calculations and flux calculations. 

**dcfunctions.h** - This file contains function headers for the file dcfunctions.C 

**dcmain.C** – This is the main program used for the Divergence Cleaning 2DMHD solver. 

**inputs.txt** - Textfile containing the settings parameters for the test to be simulated in the Divergence Cleaning 2D MHD solver. 

**master.txt** - Textfile containing the settings parameters for all of the 2D validation tests considered in this written assignment. The settings for a particular test can be copied from the master.txt file to the inputs.txt file when needed. 

**GNUmakefile** – Builds an executable named ‘divmhd’ which can be compiled by typing make in the directory.   

### Exact Solutions for Sods Test

To compute the exact solutions for Sods test to include in the plots in my report, I have created a program which takes in star values for pressure, velocity and density from Table 4.2 in Toro’s textbook ‘Riemann Solvers and Numerical Methods for Fluid Dynamics’ and runs a sampling function to produce to the exact solution.  

**samplingmain.C** – This file reads in inputs from a setting file and reads in x-values from a data file and then samples the exact solution at each of these x-values. The program outputs a data file which contains the x-values and the exact solution at these values. 

**samplingfunctions.C** – This file contains the sampling function and function for reading in parameters from the settings file.

**samplingfunctions.h** – This file contains function headers for the file samplingfunctions.C .

**inputs.txt** – Textfile containing the settings parameters for computing the exact solution for one of Toro’s five tests. 

**master.txt** – Textfile containing the settings parameters for each of Toro’s five tests. The settings for a particular test can be copied from the master.txt file to the inputs.txt file when needed. 

**GNUmakefile** – Builds an executable named ‘mysampling’ which can be compiled by typing make in the directory.   
