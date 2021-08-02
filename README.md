# ToRNAdoTeam-SRP

This repository contains files for the Re-Analysis of Smajic et. al. scRNA-seq study, for the preprocessing and the website. The original article can be found here https://www.medrxiv.org/content/10.1101/2020.09.28.20202812v1.full while their github repository is found here https://github.com/SpielmannLab.

| We wanted to reproduce the results, that the original authors were able to create through using single cell RNA-sequencing technique on 11 brain midsection samples (5 with Idiopathic Parkinsons Disease and 6 controls). The samples were taken from postmortem human brains and were prepared accordingly for usage on scRNA-seq machines. Our goal here is to do an *exact* step by step computational re-analysis. This includes download sequences from the SRA, converting them into appropriate binary files and processing them with CellRanger tools. After the pre-processing part, we used the scripts that the authors provided on their github page for cell clustering, cell sub-clustering, differential gene expression and trajectory analysis. We split the authors pipeline into three phases. For the alternative analysis we used Seurat and ScanPy. Both included manually coded semi-pipelines which takes their UMI count matrix as an input.

> Phase 1: Raw data pre-proccessing

*This was done via the High Performance Cluster that belongs to University of Leicester*
Scripts for this can be found inside the raw.tar.gz archive. All scripts leading up to merging are part of the data pre-processing.

We needed to change some of the Rscripts used by authors to better suite our machines.

> Phase 2: Re-analysis of data using R

All the Rscripts were posted for free re-use by the authors. Scripts such as merging, clustering, etc. inside the raw.tar.gz archive are part of the original recreation of authors analysis. 

The alternative anlysis can be found in the other .tar.gz archive, where there is only one file called r_analysis_byauthor.R and for the ScanPy based it is called scanpy_analysis.py.

For magma there is an Rscript called magma_prepar.R inside the same .tar.gz archive. The script prepares files for the first step of the MAGMA analysis (which yields SNPs annotated to the reference genome hg37). This is where the recreation ended since we could not reproduce steps 2 and 3. The MAGMA_script is basically just a copy of the authors, where for the gene annotation part, you have to change names and directories accordingly. 

> Phase 3: Website Construction

We constructed a website which presents the graphs/outcomes of our re-analysis. The user can compare the results by selecting/choosing which result they want to see and compare with. The results of individual genes can be searched for using our search engine. 

Our website can be found on following this link: https://www406.lamp.le.ac.uk/index.html (You have to login via University of Leicester HPC)
