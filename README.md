# CyclomicsManuscript README
Data analysis scripts used to process and plot the data for the manuscript of CyclomicsSeq.

## Date
1 July 2020

## Author
Myrthe Jager (contact: m.jager-2@umcutrecht.nl)
Alessio Marcozzi (contact: alessio.marcozzi@gmail.com)

## Description
Here, you can find all code (R scripts and Python notebooks) that can be used to reproduce the results of all the paper Figures and most Supplementary Figures of:
"Accurate detection of circulating tumor DNA using nanopore consensus sequencing".


## Data
Processed data https://doi.org/10.5281/zenodo.3925250


## Instructions for use

### A) Start up (estimated time < 1 hour)
1. Download and install required packages (see below @ Software dependencies, R packages and Python libraries)
2. Download the processed data from Zenodo (DOI: 10.5281/zenodo.3925250)
3. Create output directory '/Results/' in the main downloaded folder 'Cyclomics_manuscript' with subfolders:
- /10reps/
- /10reps/basetype/
- /10reps/extra/
- /10reps/insert/
- /consensus_loop/
- /consensus_loop/backbone/
- /consensus_loop/insert/
- /consensus_loop/all/
- /consensus_loop/all/error/
- /consensus_loop/all/qscore/
- /patients/
- /patients/patienta/
- /patients/patientb/
- /patients/patientc/
4. Add 'dir' for each script used in steps B and C (except functions.R) @ GET STARTED step 2: directory of downloaded processed data (once per script) 



### B) snFP and deletion analyses (estimated time ~ 1 hour)
5. Run 'functions.R' in RStudio 
6. Run 'consensus_loop.R' in RStudio to get & plot (~ 35 minutes runtime): 
- Figure 3a
- Figure 5a
- Supplementary figures: 4a, 7a, 7b, 8a, 8b and 9b
- Analyses of single-nucleoide FP rate and deletion rate per number of repeat for 1-39 and 40+ repeats
- Separate analyses for entire runs and per insert, and for error and qscore
- Coverage of all runs per number of repeat for 1-39 and 40+ repeats
7. Clear RStudio workspace
8. Run 'functions.R' in RStudio
9. Run '10reps.R' in RStudio to get & plot (~ 25 minutes runtime): 
- Figure 3 b,c,d
- Figure 5 b,c
- Supplementary figures: 4b-g, 5b-e, 6, 7c-f and 8c-f
- Analyses of single-nucleotide FP rate and deletion rate per nucleotide position and per run (for 10+ repeats)
- Nucleotide positions which require snFP forward/reverse correction
- Forward/reverse corrected analyses of single-nucleotide FP rate and deletion rate per nucleotide position and per run (for 10+ repeats)
- Percentage of bases with <0.1%, 0.1-1.0 and >1.0% snFP and deletion rates per run
10. Clear RStudio workspace



### C) Patients (estimated time ~ 2 minutes)
11. Run 'functions' in RStudio 
12. Run 'patientA.R' in RStudio to get & plot (~30 seconds runtime)
- Figure 6 a,d
13. Clear RStudio workspace
14. Run 'functions'
15. Run 'patientB.R' in RStudio to get & plot (~30 seconds runtime)
- Figure 6 b,e
16. Clear RStudio workspace
17. Run 'functions'
18. Run 'patientC.R' in RStudio to get & plot (~30 seconds runtime)
- Figure 6 c,f
19. Clear RStudio workspace



### D) CyclomicsSeq run-stats (estimated time < 1 hour)
The notebook `stats_from_structure.ipynb` was used to generate the plots of Figure 1 and the ones in `SuppementaryData.zip`.
1. Start a jupyter notebook server and open `stats_from_structure.ipynb`
2. Set the value of the `data_folder` variable to match the location where you have downloaded the processed data from Zenodo (DOI: 10.5281/zenodo.3925250)
3. Optionally, modify the `samples` variable to include/exclude specific runs. The default value include all the runs published in "Accurate detection of circulating tumor DNA using nanopore consensus sequencing"
4. Specify the value of `save_plots_folder`, which indicates the folder where the plots will be saved
5. Run all the cells



### E) TP53 exons coverage (estimated time < 10 minutes)
The notebook `TP53_panel_coverage.ipynb` was used to generate the plots of Figure 2.
1. Start a jupyter notebook server and open `TP53_panel_coverage.ipynb`
2. Set the value of the `data_folder` variable to match the location where you have downloaded the processed data from Zenodo (DOI: 10.5281/zenodo.3925250)
4. Set the value of the  `output_folder` variable which indicates the folder where the plots will be saved
5. Optionally, modify the `samples` variable to include/exclude specific runs
6. Run all the cells



### F) False positive rates on TP53 COSMIC mutations (estimated time < 10 minutes)
The notebook `COSMIC_analysis.ipynb` was used to generate the plots of Figure 4.
1. Start a jupyter notebook server and open `COSMIC_analysis.ipynb`
2. Set the value of the `data_folder` variable to match the location where you have downloaded the processed data from Zenodo (DOI: 10.5281/zenodo.3925250)
3. Optionally, modify the `samples` variable to include/exclude specific runs
4. Optionally, modify the `save_as` variable to change the name and the extension of the output file
5. Run all the cells



### G) Visualize the structure of a CyclomicsSeq read (estimated time < 1 hour)
The CyclomicsSeq concatemers can be visualized using `plot_read_structure.ipynb`
1. Start a jupyter notebook server and open `plot_read_structure.ipynb`
2. Set the value of the `data_folder` variable to match the location where you have downloaded the processed data from Zenodo (DOI: 10.5281/zenodo.3925250)
3. Set the value of the `sample_name` variable to specify the run to be analyzed
4. Optionally, set the other input parameters specified in the second cell of the notebook
5. Run all the cells




## Software dependencies
R version 3.5.1
RSrudio version 1.1.456

Python version 3.6
jupyterlab version 1.2.6



## R packages
- cowplot		0.9.3
- ggplot2		3.0.0
- gridExtra		2.3
- reshape2		1.4.3
- scales		0.5.0



## Python libraries
- matplotlib	3.1.3
- numpy			1.18.5
- pandas		1.0
- scipy			1.4
- biopython		1.7
