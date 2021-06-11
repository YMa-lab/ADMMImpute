# ADMMImpute
This package perform imputation on scRNAseq dataset using ADMM algorithm. 
## Getting Started
### Prerequisites
You need to install devtools, Rcpp,RcppArmadillo

```
install.packages("devtools")
install.packages("Rcpp")
install.packages("RcppArmadillo")
```

### Installing
```
library(devtools)
devtools::install_github("YingMa1993/Biostat815_Final_Project")
library(ADMMImpute)

```
## Running the tests
Here are the code scripts for the imputation on the example dataset in the package
```
#### load example data, from Chu, et al, Genome Biology, 2016.
data(Chu)
#### Get cell type information, this step can be negligible if there is no cell type information provided. 
celltype = gsub("(.+?)(\\_.*)", "\\1", colnames(Chu))
celltype = as.factor(celltype)
celltype = as.numeric(celltype)

```
```
#### Perform Imputation Analysis
imputed_count = Impute(raw_data = as.matrix(Chu),          ### raw count
                       labeled = TRUE,                     ### if it is false, then celltype information not needed
                       labels = celltype,                  ### cell type information
                       numCluster = 7,                     ### number of cell subpopulations, can be based on prior knwoledge
                       drop_thre = 0.8,                    ### dropout probability threshold set on
                       rho = 10,                           ### step-size
                       max_iter = 1000,                    ### max iteration
                       tol = 1e-04)                        ### tolerance
#### You can use the imputed_count to do the downstream analysis.

```
## Contributing

Please read [ADMMImpute-package.Rd](https://github.com/YingMa0107/Biostat815_Final_Project/man) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Ying Ma** - *815 Final Project* - [Biostat815_Final_Project](https://github.com/YingMa0107/Biostat815_Final_Project)
