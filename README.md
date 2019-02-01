## Using Cluster Ensembles to Identify Psychiatric Patient Subgroups
This is the code that belongs to the paper that has been submitted for publication under the above title. Details are in the manuscript. Data cannot be included due to legal and patient privacy constraints. 

The code consists of three files:
1. `cluster_ensemble_modeling.r` An R script, in which selection of number of clusters (using `NbClust` package), cluster ensemble modelling (using `DiceR`) package, and significance testing of clusters (using `sigclust` package) is done. 
2. `cluster_ensemble.ipynb` A Jupyter notebook that consists the loading and preprocessing of YSR data, and computation of hopkins statistic. After this, the R script (see above, 1) should be executed as a sort of subroutine. Then the notebook takes over again, by creating visualizations, and integrating with clinial notes, DSM diagnosis and other clinically relevant variables. 
2. `hopkins_statistic.py` A script to compute Hopkins statistic in Python
