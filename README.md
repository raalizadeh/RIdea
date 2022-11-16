# Retention Index based on Cocamide Diethanolamine homologous series (C(n)-DEA)
<img
     src= https://github.com/raalizadeh/RIdea/raw/main/v1.0/Win/doc/toc_DEA.png width=450 />

This repository contains the models and codes for generating experimental and predicted liquid chromatographic retention index values for emerging contaminants. The RI values are indexed tR values that are created based on Cocamide Diethanolamine homologous series (CDEA) and LC-system independent. Moreover, RI bank based on CDEA indexing system is available under “RI bank tab” for more than 3000 emerging contaminants in the app. The proposed RI system can be used to enhance identification confidence in suspect and nontarget screening while facilitating successful comparability of tR data between various LC settings.  

## System requirements
The application runs on 64-bit windows operating systems. The linux version of app is not available now because there is a discrepancy in generation of png format between linux and windows.  

## How to install  
To run the app, you need to:  
Clone this repository and find related version.  
setup R and the third-party R packages dependencies.  
The following R packages should be installed before running the app:  

[shiny](https://cran.r-project.org/web/packages/shiny/index.html), [rJava](https://cran.r-project.org/web/packages/rJava/index.html), [rcdk](https://cran.r-project.org/web/packages/rcdk/index.html), [ChemmineR](https://www.bioconductor.org/packages/release/bioc/html/ChemmineR.html), [DT](https://cran.r-project.org/web/packages/DT/index.html), [rapportools](https://cran.r-project.org/web/packages/rapportools/index.html), [shinyjs](https://cran.r-project.org/web/packages/shinyjs/index.html), [cluster](https://cran.r-project.org/web/packages/cluster/index.html), [mlbench](https://cran.r-project.org/web/packages/mlbench/index.html), [drat](https://cran.r-project.org/web/packages/drat/index.html), [mxnet](https://github.com/apache/mxnet), [EBImage](https://www.bioconductor.org/packages/release/bioc/html/EBImage.html), [pbapply](https://cran.r-project.org/web/packages/pbapply/index.html), [neuralnet](https://cran.r-project.org/web/packages/neuralnet/index.html), [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html), [caret](https://cran.r-project.org/web/packages/caret/index.html), [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html), [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)  

Note: [mxnet](https://github.com/apache/mxnet) R package needs to be rebuilt under your operating system, as pre-built version might be unavailable. Please, kindly follow the [tutorial](https://mxnet.apache.org/get_started/build_from_source) and build mxnet and install it for R. In addition, the app is written and tested using R version 3.6.3.

Last but not least, the app requires “www/checkmol.exe” to complete the generation of fingerprints map. Checkmol was written by Dr. Norbert Haider, which it can be accessed free of charge for academic use from [here](https://homepage.univie.ac.at/norbert.haider/cheminf/cmmm.html). Checkmol is a software that creates fingerprints for various organic functional groups. 


## Quick Tutorial

Access the manual in "About" tab from app, follow the steps and apply the procedure for any compound of interest. The manual is also available in “www/manual.pdf”. In order to be sure that the app is running correctly, please, verify that the calculations are the same for the following three compounds under your operating system:  

Example: SMILES tR  
c1ccc(c(c1)CC(=O)OCC(=O)O)Nc2c(cccc2Cl)Cl 9.2  
COC1=C(OC)C=C(C(C2=C(N)N=C(N)N=C2)=O)C=C1OC 5.49  
CCCCCCCCn1c(=O)c(c(s1)Cl)Cl 12.49  

 <img
     src= https://github.com/raalizadeh/RIdea/raw/main/v1.0/Win/doc/exam_DEA.png width=850 />

## Purchase of RI-CDEA standards mix
The RI calibrants based on CDEA can be purchased from [TrAMS group](http://trams.chem.uoa.gr/).  
Please, contact us for more details.  

## Authors
Reza Aalizadeh, Email: raalizadeh@chem.uoa.gr  
Varvara Nikolopoulou, Email: vnikol@chem.uoa.gr  
Nikolaos S. Thomaidis, Email: ntho@chem.uoa.gr 

Laboratory of Analytical Chemistry,  
Department of Chemistry,  
National and Kapodistrian University of Athens,  
Panepistimiopolis Zografou, 15771, Athens, Greece  

## Publication

Development of Liquid Chromatographic Retention Index Based on Cocamide Diethanolamine Homologous Series (C(n)-DEA), Analytical Chemistry, 2022. DOI: [10.1021/acs.analchem.2c02893](https://pubs.acs.org/doi/10.1021/acs.analchem.2c02893)
