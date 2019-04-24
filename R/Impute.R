#' A R function to apply imputation_model which was implemented in Rcpp functions to perform imputation on scRNA-seq data
#'
#' @param raw_data The raw data in a matrix format you want to perform imputation on. And the rownames of the count matrix should representing the genes and column names should representing the cells.
#' @param labeled A logical value describing whether cell labels information are known.
#' @param labels A character vector describing the cell label information. Only needed when \code{labeled = TRUE}.
#' @param numCluster A number describing the expected or estimated number of subpopulations. This could be determined by prior knowledge or unsupervised clustering.
#' @param drop_thre A number between 0 and 1 describing the dropout rate threshold to determine the dropout genes and the non-dropout genes.
#' @param rho step size in ADMM algorithm, default value is 10.
#' @param max_iter max iterations in ADMM algorithm, default value is 5000
#' @param tol tolerance in ADMM algorithm, default value is 1e-4
#' @return Impute returns imputed count matrix which matches the dimension of the original count matrix
#' @useDynLib ADMMImpute
#' @importFrom Rcpp sourceCpp
#' @importFrom rsvd rpca
#' @importFrom kernlab specc
#' @importFrom utils read.csv read.table write.csv write.table
#' @export
Impute <- function (raw_data, labeled = FALSE,labels = NULL, numCluster, drop_thre, rho = 10,max_iter = 5000,tol = 1e-4)
{	
    print("reading in the raw count matrix ...")
    raw_count = as.matrix(raw_data)
    gene_name = rownames(raw_count)
    cell_name = colnames(raw_count)
    totalC_PerCell = colSums(raw_count)
    totalC_PerCell[totalC_PerCell == 0] = 1
    count_norm  = sweep(raw_count, 2, 10^6/totalC_PerCell, FUN = "*")
    count_norm = log10(count_norm + 1.01)
    if(labeled == FALSE){
         pca = rpca(t(count_norm), k = 100, center = TRUE, scale = FALSE)
         eigvalues = (pca$sdev)^2
         cumVar = cumsum(eigvalues)/sum(eigvalues)
         var_thre = 0.6
         if(max(cumVar) <= var_thre){
             npc = length(cumVar)
         }else{
             npc = which.max(cumVar > var_thre)
             npc = max(npc, numCluster)
         }
         pcscore = t(pca$x[, 1:npc])
         res = specc(t(pcscore),numCluster,kernel = "rbfdot")
         cluster =  res
    }else{
    if(is.null(labels)){
            stop("Please provide the cell lables if 'labeled = TRUE'")}else{
                if(length(labels) != ncol(count_norm)){
                    stop("The number of cell labels does not match the number of cells in raw data!")
                }
                cluster = labels
            }
    }
    print("Imputation starts ...")
    count_norm = as.matrix(count_norm)
    #imputed_count = try(imputation_wlabel_model(count_norm,cluster,log10(1.01),drop_thre,10,1000,1e-3))
    imputed_count = ADMMImpute::imputation_wlabel_model(count_norm,cluster,log10(1.01),drop_thre,rho,max_iter,tol)
    imputed_count = 10^imputed_count - 1.01
    count_imp = sweep(imputed_count, MARGIN = 2, colSums(raw_count)/10^6,FUN = "*")
    count_imp = round(count_imp, digits = 2)
    rownames(count_imp) = gene_name
    colnames(count_imp) = cell_name
    print("Imputation finished")
    return(count_imp)
}




