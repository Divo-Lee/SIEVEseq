# SIEVE

### SIEVE: One-stop differential expression, variability, and skewness analyses using RNA-Seq data

Hongxiang Li and Tsung Fei Khang

SIEVE: A statistical method for building a unified framework for the simultaneous testing of differential expression, variability and skewness of genes using RNA-Seq data. This framework adopts a compositional data analysis approach to modelling RNA-Seq count data, applies the centered log-ratio (CLR) transformation to convert them into continuous variables, and uses a skew-normal distribution to model them.


### Installation:
Install SIEVE from local source with

`install.packages("SIEVE_1.0.0.tar.gz", repos=NULL, type="source")`

Install SIEVE from GitHub with

 `if (!"devtools" %in% installed.packages()) {
  install.packages("devtools")}`
  
 `devtools::install_github("Divo-Lee/SIEVE")`
 
 
### Dependencies:
 `SIEVE` `R` package depends on the following packages: `sn`, `stats`, `utils`, `ggplot2`.
