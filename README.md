# Protein-DNA Binding Model

## About

The given package provides a statistical-mechanical model for calculation of non-specific protein-DNA binding free energies from experimental data. While the complete description of the Mathematical formulation can be found in the original publication **here**, a brief overview is as follows -

For a DNA fragment with `n` base pairs, and a protein that binds to `p` continuous base pairs on the DNA in one of `m` modes, the model
1. Algorithmically enumerates all possible binding configurations
2. Calculates partition function and therefore average number of bound proteins at a given concentration
3. Fits the binding energy per mode from experimental anisotropy measurements

Certain limitations and constraints currently baked into the model -
1. The protein can bind in maximum of two modes (`m = 1 or 2`).
2. The DNA strand cannot allow more than two proteins to bind to it (`n < 3p`).
3. Proteins cannot partially bind to the DNA fragment; a protein strictly occupies `p` continuous sites when bound.
4. When multiple proteins can bind in two modes, more than one protein cannot bind in the major mode (i.e. the mode with greater binding free energy).

An example overview of the workflow for a simplified system with `(n, p, m) = (5, 2, 2)` is illustrated below -

![Model overview](/images/figure_1.png)

## Installation Instructions

### Software instructions

Following version of softwares are required -
- MATLAB v24.2.0.2773142 (R2024b) Update 2

## Usage

### Overview

The user must provide the correct input data file, and set up model parameters as required, prior to executing a run.

**Input data -** Must be provided in a `.dat` format and placed inside the `./data/` folder. The file must contain three data columns -
1. Concentration (`µmol L^-1`)
2. Experimentally measured anisotropy
3. Standard deviation of measurements

Incorrect input data shape is flagged by the code.

**Model parameters -** These can be modified in lines `22-26` in the `./main.m` script. These are -
1. `pname` - Protein name used for saving output files. Note that the code searches for the input data file at `./data/{pname}.dat`, so user must set `pname` to have the same name as the input data (or equivalently rename the input data of course).
2. `n` - The number of base pairs on the DNA segment
3. `T` - The temperature at which the system is to be simulated (in Kelvin)
4. `sites` - The number of continuous base pairs that the protein binds to
5. `modes` - The number of binding modes of the protein

The user can then run `./main.m` script. Note that the working directory must be the root directory itself (i.e. the folder with the `main.m` script itself).

**Output files -** During the model execution, the output folder is created as `./results/{pname}/`, and the output files are stored in it. The results will be overwritten if another run is started with the same `pname`, or else the new results are stored in a new folder. Three output files are generated -

1. `{pname}.xlsx` - Stores the following calculated fields at the end of model run -
a) Concentration in `µmol L^-1`
b) The calculate anisotropy values at input concentrations
c) Signal contribution from fragments with only 1 protein bound
d) Signal contribution from fragments with 2 proteins bound
e) The average number of proteins bound at that concentration
2. `{pname}.png` - Stores the output plot as an image, which contains 3 plots -
a) Concentration (`µmol L^-1`) vs Anisotropy
b) Concentration (`µmol L^-1`) vs Occupancy (no. of proteins bound)
c) A 3D-plot of the binding probability at each base pair of the DNA, at different concentrations
3. `{pname}.mat` - Stores the following workspace variables for future reference -
a) `pname` - The given protein name for the run
b) `concs` - The input experimental concentrations
c) `exp_S` - The input experimental anisotropy values
d) `err_S` - The input standard deviations
e) `calc_S` - The calculated anisotropy values at the input concentrations

The code for plotting the graphs is between lines `95 - 156`, and the code for saving output files is between lines `157 - 176`, should the user wish to make any changes.

### Examples

Four example input files (`NHP6A.dat`, `NHP6A_S26D.dat`, `NHP6A_T63D.dat` and `NHP6A_TM.dat`) ae pre-loaded into `./data/`. The user can run these examples by setting `pname` and other input parameters accordingly in `./main.m`.

As a rule of thumb, a single binding mode is chosen when the experimental data appears monophasic, and two binding modes are chosen when the experimental data shows biphasic behaviour. However, both modes can be tested, and the better fit can be chosen by the user basis their context. Other parameters are set based on known protein/DNA sequences and their interactions, but, again, the best fits can be manually selected by the user over a range of candidate input parameters.

## References

If using our protein-DNA binding model, please cite the following articles:

1. 

2. 

Additional reading and references:

1. 

2. 

3. 

4. 

5. 
