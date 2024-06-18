# ORBproject
This repository is a refinement of the "ORBmeta" and "Simulations" repositories, i.e., the work done in the Master Thesis. The focus of the "ORBproject" repository is on outcome reporting bias (ORB) adjustment (and associated simulations according to the protocol) via a selection model framework in the random effects model meta-analysis setting without the use of the ORBIT pre-classification system.

## Repository Organization

Below is a description of the different folders and their respective content.

### Simulation Folder

1. simulate_ORB_data.R contains the function used to simulate a meta-analysis dataset in the presence of ORB

2. reORBadjust.R contains the function used to estimate the treatment effect and heterogeneity variance of the meta-analysis, both unadjusted and adjusted for ORB according to the desired ORB-adjustement method, which is specified via the selection function

3. SummaryStatsHelper.R contains the functions used to calculate statistics summary measures in the simulations

4. ORBSimAdj_15.R is the script which conducts the simulation study for the parameter setting wherein ORB is simulated with a value of $\gamma=1.5$

5. AdjustResults3200k_15.rds is the .rds file of the resulting simulations obtained with the above script (4) and used in the final manuscript of this project

6. ORBSimAdj_05.R is the script which conducts the simulation study for the parameter setting wherein ORB is simulated with a value of $\gamma=0.5$

7. AdjustResults3200k_05.rds is the .rds file of the resulting simulations obtained with the above script (6) and used in the final manuscript of this project

8. PlotResults.R contains the functions and code to be executed, so as to produce the plots of the summary measures obtained from the simulations, for both simulations processes obtained from point scripts (4) and (5) above-mentioned. These plots are for the main parameter of interest, i.e., the treatment effect estimate

9. PlotResults_Tau2.R contains the same functions as script (8) above, but for the secondary parameter of interest in our simulation, i.e., the heterogeneity variance


### Paper Folder

The folder contains the file ORB2.Rnw file, which can be compiled in R, so as to obtain the ORB2.pdf file, also included in this folder. The outputted pdf file is the main file of the manuscript for this project, dynamically producing text, tables, and figures obtained from the simulations. The folder additionally includes the bibliography file biblio.bib, as well as the figures and cache in the respective folders.

### Protocol Folder

The folder contains the file SimulationAnalysisProposalORB.Rnw file, which can be compiled in R, so as to obtain the SimulationAnalysisProposalORB.pdf file, also included in this folder. The outputted pdf file is the simulation protocol file for this project, dynamically producing text, tables, and figures. The simulations conducted in the Simulation Folder, the results of which are used in the Paper Folder, are obtained following this pre-specified protocol. This folder additionally includes the bibliography file biblio.bib.

## Editing and Compiling

In order to make edits to the manuscript, download this repository as a .zip file, and open the project in R studio, where all the contents of the repository will be shown. Here, changes to  the main manuscript file can be made directly by editing the ORB2.Rnw file and subsequently compiling it to a pdf from the R studio interface.

## Running Simulations

In order to re-run simulations, execute the files ORBSimAdj_15.R and ORBSimAdj_05.R in R or in the desired server. The simulations will output .rds files with the final data. These files, defined in (5) and (7) of the Simulation Folder above, are called in the ORB2.Rnw file to generate the desired plots for the manuscript.


