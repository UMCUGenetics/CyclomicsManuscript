# CyclomicsManuscript README
Data analysis scripts used to process and plot the data for the manuscript of CyclomicsSeq.

## Date
1 July 2020 <br />

## Author
Myrthe Jager (contact: m.jager-2@umcutrecht.nl) <br />

## Description
Here, you can find all code that can be used to reproduce the results of Figures 3, 5 and 6 and most Supplementary Figures of: <br />
"Accurate detection of circulating tumor DNA using nanopore consensus sequencing"

## Data
Processed data DOI: 10.5281/zenodo.3925250


## Instructions for use

### A) Start up (estimated time < 1 hr)
1. Download and install required packages (see below @ Software dependencies & R packages)
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


## Software dependencies
R version 3.5.1<br />
RSrudio version 1.1.456 <br />


## R packages
- cowplot	0.9.3
- ggplot2	3.0.0
- gridExtra	2.3
- reshape2	1.4.3
- scales	0.5.0
