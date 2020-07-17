# Taylor-et-al-2020-demo

Machine learning demo code for Taylor et al "The molecular basis of socially-mediated phenotypic plasticity in a eusocial paper wasp".

# System requirements

Requires R version v3.6+ with packages *tidyverse*, *e1071*, *DESeq2*, *limma* & *pracma* installed.

Tested in Ubuntu v18.04.4 using R v3.6.3 running in RStudio v1.3.959.

# Installation guide

We provide three sets of gene expression data for use with the demo code:

  * *counts_raw* is the full set of gene expression data generated for the experiment
  * *counts_clean* is a 'cleaned' set of genes from which rows representing genes with low expression across all samples have been removed
  * *counts_clean_subsample* is a smaller set of 1000 genes randomly picked from those in *counts_clean*; we recommend using this smaller dataset to reduce computation times when testing the demo.
  
Each dataset is provided both as an RData object which can be loaded directly into R and also as a .csv file for export to other programs. Additionally, to use the demo you will need to load the object *phenotypic_data* which describes which samples belong to which treatment groups in addition to phenotypic data (the latter are not used in the present demo but are included in case they are of interest). 

# Demo & instructions for use

Once a gene expression data matrix and the phenotypic data matrix have been loaded into an R session, the demo code can be run without further intervention. As output, the demo produces:

 * Classifications and error estimates for a full svm model run using the full set of gene expression data provided
 * A rudimentary plot of the results of feature selection performed upon the provided dataset- this plot should possess a characteristic 'hockeystick' shape, with the error rates of the SVM models gradually decreasing at first then rapidly increasing once caste-informative features begin being removed
 * Classifications and error estimates for the optimised svm model produced using just the genes identified by feature selection as being caste-informative
 
 Note that the demo will likely take several hours to run on a desktop computer if the larger datasets provided are used. We recommend using the smaller subsample provided for testing the demo, which should take under an hour to run. This run time can be further reduced if the user modifies the demo to reduce the number of error estimates made in each loop of the feature selection process (set at 20 by default). Setting the number of loops to one will allow the demo to be run in under 10 minutes if using the subsetted data, but do note that the resulting error estimates will exhibit higher rates of stochasticity.  
