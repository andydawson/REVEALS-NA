# REVEALS-NA
Workflow to estimate Holocene land cover from fossil pollen data using REVEALS.

## Workflow

`1_prepare_pollen.R`

Prepare the pollen data in order for REVEALS.

`2A_run_reveals_parallel.R`

Prepare REVEALS inputs and run REVEALS for all sites. The code uses parellel computation to do this, and generates an output file with REVEALS estimates for each record.

`3_reveals_amalgomate.R`

Amalgomate the individual record REVEALS estimate into one data frame.   

`4_reveals_to_lct.R`

Translate the REVEALS estimates to land cover estimates in the Summergreen Trees and Shrubs (STS), Evergreen Trees and Shrubs (ETS), and Open Vegetated Land (OVL) categories.
