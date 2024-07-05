# ORBproject
This `ORBproject` repository contains the code and source files of the paper "Addressing Outcome Reporting Bias (ORB) in Meta-analysis: A Selection Model Perspective". The repository, and, overall, the project, can be seem as a refinement of the `ORBmeta` and `Simulations` repositories, i.e., the work on ORB done in the Master Thesis. 

The focus of the `ORBproject` repository is on outcome reporting bias (ORB) adjustment (and associated simulations according to the simulation study protocol) via a selection model framework in the random effects model meta-analysis of a beneficial outcome.

## Repository Organization

Below is a description of the different folders and their respective contents.

### Simulation Folder

1. `simulate_ORB_data.R` contains the function used to simulate a meta-analysis dataset in the presence of ORB

2. `reORBadjust.R` contains the main function used to estimate the treatment effect and heterogeneity variance of the meta-analysis, both unadjusted and adjusted for ORB according to the desired ORB-adjustement method, which is specified via the selection function

3. `SummaryStatsHelper.R` contains the functions used to calculate statistics summary measures in the simulations

4. `ORBSimAdj_15.R` is the script which conducts the simulation study for the parameter setting wherein ORB is simulated with a value of $\gamma=1.5$

5. `AdjustResults3200k_15.rds` is the .rds file of the resulting simulations obtained with the above script (4) and used in the final manuscript of this project

6. `ORBSimAdj_05.R` is the script which conducts the simulation study for the parameter setting wherein ORB is simulated with a value of $\gamma=0.5$

7. `AdjustResults3200k_05.rds` is the .rds file of the resulting simulations obtained with the above script (6) and used in the final manuscript of this project

8. `PlotResults.R` contains the functions and code to be executed, so as to produce the plots of the summary measures obtained from the simulations, for both simulations processes obtained from points (4) and (5) above-mentioned. These plots are for the main parameter of interest, i.e., the treatment effect estimate

9. `PlotResults_Tau2.R` contains the same functions as point (8) above, but for the secondary parameter of interest in our simulation, i.e., the heterogeneity variance


### Paper Folder

The folder contains the file `ORB2.Rnw` file, which can be compiled in R, so as to obtain the `ORB2.pdf` file, also included in this folder. The outputted pdf file is the main file of the manuscript for this project, dynamically producing text, tables, and figures obtained from the simulations. The folder additionally includes the bibliography file `biblio.bib`, as well as the figures and cache in the respective folders.

### Protocol Folder

The folder contains the file `SimulationAnalysisProposalORB.Rnw` file, which can be compiled in R, so as to obtain the `SimulationAnalysisProposalORB.pdf` file, also included in this folder. The outputted pdf file is the simulation protocol file for this project, dynamically producing text, tables, and figures. The simulations conducted in the Simulation Folder, the results of which are used in the Paper Folder, are obtained following this pre-specified protocol. This folder additionally includes the bibliography file `biblio.bib`.

## Editing and Compiling

In order to compile the manuscript 

1. Clone this repository  with `git clone https://github.com/agaiasaracini/ORBproject`

2. Navigate to the main manuscript folder, the Paper folder with `cd ORBproject` and then `cd Paper`

3. Use the Make file to compile the manuscript with `make` or `make -B`

4. Open the pdf files generated with `open ORB2.pdf` for the main manuscript or `open Supplementary.pdf` for the supplementary material

In order to make edits to the manuscript and/or other contents of the repository, e.g., `nano ORB2.Rnw`, after the chages have been made do the following:

`git commit -m "Changes to ORB2.Rnw"` commit changes
`git push origin main` push changes

Then repeat steps (3) and (4) above.

Alternatively, to pull changes and thus update the clones repository contents use `git pull origin main`.

The repeat steps (3) and (4) above.

## Running Simulations

In order to re-run simulations, execute the files `ORBSimAdj_15.R` and `ORBSimAdj_05.R` in R or in the desired server. The simulations will output .rds files with the final data. These files, defined in (5) and (7) of the Simulation Folder above, are called within the `ORB2.Rnw` file to generate the desired plots for the manuscript.



