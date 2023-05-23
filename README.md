# Authors

- Dr. Enrico Schiassi - PhD (University of Arizona, Tucson, AZ)
- Dr. Roberto Furfaro - Professor, Systems and Industrial Engineering, University of Arizona, Tucson, AZ
- Dr. Jeffrey S. Kargel - Senior Scientist, Planetary Science Institute, Tucson, AZ

# Installing and running `GLAM_BioLithRT`

- Download and install a recent version of Matlab (from R2015a on)
- all the files needed to run the code (database - w/ input spectra * -, Matlab functions, and Matlab script) must be in the same folder (in any location).
- functions and input spectra are into a folder called functions&input-spectra. Add its path in the main.m to run the code
- scripts and functions to perform MCMC sampling are in the folder mcmcstat. The code is available at http://helios.fmi.fi/~lainema/mcmc/

    * The Input Spectra are taken from the data base in the folder DATA available in WASI4 package (http://www.ioccg.org/data/software.html)

# Major files

- `main.m` : script to enter the input and run the code
- `AOP_Rs.m`: function to simulate the Remote sensing reflectance (Rrs) given the required input
- `InvModeBioLithRT_Copt.m`: function to compute the objective function for the constrained optimization for water components concentration retrieval
- `InvModeBioLithRT_Bopt.m`: function to compute the objective function for the Bayseian optimization for water components concentration retrieval

# Notes

- the code is in source format. The user has access and can modify all the scripts and the functions provided according to his/her tasks (where GLAM BioLithRT can accomplish those)
- if any modifications are done in `AOP_Rs.m`, the same modifications must be done in `InvModeBioLithRT_Copt.m`, `InvModeBioLithRT_Bopt.m`, and vice versa 

# Contacts

- Dr. Roberto Furfaro, robertof@email.arizona.edu
- Dr. Enrico Schiassi, enrico.schiassi@gmail.com
- Dr. Jeffrey S. Kargel, jkargel@psi.edu
