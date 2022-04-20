# ADMMImpute
This package perform imputation on scRNAseq dataset using ADMM algorithm. The model was built originally in the method scImpute [https://www.nature.com/articles/s41467-018-03405-7]. Now I reimplemented this model with a faster optimization framework by using ADMM.
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
devtools::install_github("YingMa0107/ADMMImpute")
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
library(peakRAM)
mem1 <- peakRAM({
imputed_count = Impute(raw_data = as.matrix(Chu),          ### raw count
                       labeled = TRUE,                     ### if it is false, then celltype information not needed
                       labels = celltype,                  ### cell type information
                       numCluster = 7,                     ### number of cell subpopulations, can be based on prior knwoledge
                       drop_thre = 0.8,                    ### dropout probability threshold set on
                       rho = 10,                           ### step-size
                       max_iter = 1000,                    ### max iteration
                       tol = 1e-04)                        ### tolerance
})
#### print the computation time abd efficiency that was cost by ADMMImpute
print(mem1$Peak_RAM_Used_MiB)
890.1  ### ADMMImpute use 890.1 MiB!
print(mem1$Elapsed_Time_sec)
67.966 ### ADMMImpute cost 70 seconds!
#### bechmarking the computation time and memory with scImpute
library(devtools)
install_github("Vivianstats/scImpute")
library(scImpute)
mem2 <- peakRAM({
scimpute(count_path = "./chu.csv",infile = "csv",  outfile = "csv",out_dir = "./",    labeled = FALSE,   drop_thre = 0.8, Kcluster = 7,ncores = 1)
})
#### print the computation time abd efficiency that was cost by scimpute
print(mem2$Peak_RAM_Used_MiB)### scImpute use 992.3 MiB
print(mem2$Elapsed_Time_sec)### scImpute use cost seconds, which is 54 times slower than ADMMImpute
#### compare the results of these two methods
scImpute_imputed_count = read.csv("./scimpute_count.csv",row.names = 1)
ADMMImpute_imputed_count = imputed_count
#### print the correlation of these two impited count
cor(as.numeric(as.matrix(ADMMImpute_imputed_count)),as.numeric(as.matrix(scImpute_imputed_count)))
0.9988723
```
## Contributing

Please read [ADMMImpute-package.Rd](https://github.com/YingMa0107/Biostat815_Final_Project/tree/master/man/) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Ying Ma** - *815 Final Project* - [Biostat815_Final_Project](https://github.com/YingMa0107/Biostat815_Final_Project)
