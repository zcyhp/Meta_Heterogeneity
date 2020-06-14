# Meta_Heterogeneity
This repository contains simulation results and R code for the paper "A Comparison of Hypothesis Tests for Homogeneity in Meta-analysis with Focus on Rare Binary Events".
## All simulation results
All simulation results are stored in excel files with each file for one <img src="https://latex.codecogs.com/gif.latex?\tau^2" title="\tau^2" /> setting in {0, 0.1,..., 1}. In each excel file, results for K=10, 20, 50 are stored in 3 seperate sheets.
## Figures
  * Each figure contains 4 heat maps showing adjusted power values for all the tests under the following 4 corresponding settings:
  
    * <img src="https://latex.codecogs.com/gif.latex?\mu=-2.5\text{,&space;large-sample}" title="\mu=-2.5\text{, large-sample}" />
    * <img src="https://latex.codecogs.com/gif.latex?\mu=-5\text{,&space;large-sample}" title="\mu=-5\text{, large-sample}" />
    * <img src="https://latex.codecogs.com/gif.latex?\mu=-2.5\text{,&space;small-sample}" title="\mu=-2.5\text{, small-sample}" />
    * <img src="https://latex.codecogs.com/gif.latex?\mu=-5\text{,&space;small-sample}" title="\mu=-5\text{, small-sample}" />
  
  * We use the combination of parameters <img src="https://latex.codecogs.com/gif.latex?\{K=20,&space;R=1,&space;\theta=0,&space;w=0\}" title="\{K=20, R=1, \theta=0, w=0\}" /> as the reference setting and generated rest figures by changing one parameter at a time. 
## R code
To implement the R code: 
  * R version 3.6.3 or above is needed.
  * The file `LOR_moments_final.txt` needs to be sourced.
