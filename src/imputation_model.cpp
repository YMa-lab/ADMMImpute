#include <RcppArmadillo.h>
#include <iostream>
#include <Rcpp.h>
#include <cassert>
#include <cmath>
#include <armadillo>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
double func(double alpha, double target){ //root-finding equation
    return(log(alpha) - R::digamma(alpha) - target);
}

// bisection to solve the uniroot problem
double bisection(double a,double b, double target)
{
    double c = a;
    int iter = 0;
    while ((b-a) >= 0.0001220703)
    {
        c = (a+b)/2;
        if (func(c,target) == 0.0){
            break;
        }
        else if (func(c,target)*func(a,target) < 0){
            b = c;
        }
        else{
            a = c;
        }
        iter = iter + 1;
        if (iter > 1000){
            break;
        }
    }
    return(c);
}// end func

// vectorised version of the power function
vec vectorise_pow(vec a, vec b) {
    vec res = arma::zeros(a.size());
    for (int i = 0; i < a.size(); ++i) {
        res(i) = pow(a(i), b(i));
    }
    return res;
}// end func

// Calculate the density of mix distribution for a vector of data and paramters
vec dmix(vec x,vec pars){
    vec prob = pars(0) * pow(x,(pars(1)-1))%exp(-x * pars(2)) * pow(pars(2),pars(1)) / exp(lgamma(pars(1))) + (1 - pars(0)) * normpdf(x, pars(3), pars(4));
    return(prob);
}// end func

// Calculate the pdf of gamma distribution for a vector of data and parameters
vec dgamma(vec x,vec shape, vec scale){
    vec prob = vectorise_pow(x, (shape-1)) % exp(-x % scale) % vectorise_pow(scale,shape) / exp(lgamma(shape));
    return(prob);
}// endfunc

// Calculate the pdf of normal distribution for a vector of data and parameters
vec dnorm(vec x,vec mean, vec sd){
    vec prob = normpdf(x, mean, sd);
    return(prob);
}// endfunc

// calculate cell distances
mat calc_dis(mat A){
    colvec An =  sum(square(A),1);
    mat C = -2 * (A * A.t());
    C.each_col() += An;
    C.each_row() += An.t();
    return(sqrt(C));
}//endfunc

// get quantile
int get_quantile(vec num, double q){
    vec temp = num;
    sort(temp.begin(), temp.end());
    int quantile_idx = temp.size() * q;
    return(quantile_idx);
}// endfunc

// find the different elements of set A to set B
arma::uvec std_setdiff(arma::uvec& x, arma::uvec& y) {
    
    std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
    std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
    std::vector<int> out;
    
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                        std::inserter(out, out.end()));
    
    return arma::conv_to<arma::uvec>::from(out);
}// end func

//*************************************************************************//
//                 Calculate the weight for neighbors                      //
//*************************************************************************//
mat calculate_weight(vec x, vec paramt){
    vec pz1_pz2 = dmix(x,paramt);
    vec pz1 = paramt(0) * pow(x,(paramt(1)-1))%exp(-x * paramt(2)) * pow(paramt(2),paramt(1)) / exp(lgamma(paramt(1)));
    vec pz = pz1/pz1_pz2;
    mat pz_cbind(pz.n_elem, 2);
    pz.replace(datum::nan, 0);
    pz_cbind.col(0) = pz;
    pz_cbind.col(1) = 1 - pz;
    return(pz_cbind);
}// end func

//*************************************************************************//
//                 Update paramters in gamma distribution                  //
//*************************************************************************//
vec update_gmm_pars(vec x, vec wt){ //update parameters in gamma distribution
    double tp_s = sum(wt);
    double tp_t = sum(wt % x);
    double tp_u = sum(wt % log(x));
    double tp_v = -tp_u/tp_s - log(tp_s / tp_t);
    double alpha = 0;
    double alpha0 = 0;
    vec interval = arma::zeros(2);
    colvec res = arma::zeros(2);
    if (tp_v <= 0){
        alpha = 20;
    }else{
        alpha0 = (3 - tp_v + sqrt(pow((tp_v - 3),2) + double(24) * tp_v)) / double(12) / tp_v;
        if (alpha0 >= 20){
            alpha = 20;
        }else{
            alpha = bisection(0.9*alpha0,1.1*alpha0,tp_v);
        }
    }
    double beta = tp_s / tp_t * alpha;
    res(0) = alpha;
    res(1) = beta;
    return(res);
}// end func

//*************************************************************************//
//                    Get mix paramters for one gene                       //
//*************************************************************************//
vec get_mix(vec xdata, double point){
    vec inits = arma::zeros(5);
    vec temp1 = xdata.elem(find(xdata == point));
    inits(0) = double(temp1.size())/double(xdata.size());
    if (inits(0) == 0) {
        inits(0) = 0.01;
    }
    inits(1) = 0.5;
    inits(2) = 1;
    if(xdata.elem(find(xdata > point)).is_empty()){
        inits(3) = R_NaN;
        inits(4) = R_NaN;
    }else{
        inits(3) = mean(xdata.elem(find(xdata > point)));
        inits(4) = stddev(xdata.elem(find(xdata > point)));
    }
    if (R_IsNA(inits(4))) {
        inits(4) = 0;
    }
    vec paramt = inits;
    double eps = 10;
    int iter = 0;
    double loglik_old = 0;
    while(eps > 0.5) {
        mat wt = calculate_weight(xdata, paramt);
        paramt[0] = sum(wt.col(0))/wt.n_rows;
        paramt[3] = sum(wt.col(1) % xdata)/sum(wt.col(1));
        paramt[4] = sqrt(sum(wt.col(1) % pow((xdata - paramt(3)),2))/sum(wt.col(1)));
        paramt(1) = update_gmm_pars(xdata, wt.col(0))(0);
        paramt(2) = update_gmm_pars(xdata, wt.col(0))(1);
        double loglik = sum(log10(dmix(xdata, paramt)));
        eps = pow((loglik - loglik_old),2);
        loglik_old = loglik;
        iter = iter + 1 ;
        if (iter > 100){
            break;
        }
    }
    return(paramt);
}// end func

//*************************************************************************//
//                    Get mix paramters for all genes                      //
//*************************************************************************//
mat get_mix_parameters(mat count, double point){
    mat paramt =  arma::zeros(count.n_rows,5);
    for(int irow = 0; irow < count.n_rows; ++ irow){
        vec xdata = count.row(irow).t();
        vec rowPara = get_mix(xdata, point);
        paramt.row(irow) = rowPara.t();
        double temp = abs(sum(count.row(irow)) - point * double(count.n_cols));
        if(temp < 1e-10){
            paramt.row(irow).fill(datum::nan);
        }
        if(paramt.row(irow).has_nan()){
            paramt.row(irow).fill(datum::nan);
        }
        if(paramt.row(irow)(4) < 1e-10){
            paramt.row(irow).fill(datum::nan);
        }
    }
    return(paramt);
}// end func

//*************************************************************************//
//                        Find highly variable genes                       //
//*************************************************************************//
mat find_hv_genes(mat count){
    vec mu = arma::zeros(count.n_rows);
    vec sdv = arma::zeros(count.n_rows);
    for(int irow = 0; irow < count.n_rows; ++irow){
        vec temp = count.row(irow).t();
        vec temp_sub = temp.elem(find(temp != log10(1.01)));
        if(temp_sub.is_empty()){
            mu[irow] = R_NaN;
            sdv[irow] = R_NaN;
        }else{
            mu[irow] = mean(temp_sub);
            sdv[irow] = stddev(temp_sub);
            
        }
    }
    mu.replace(datum::nan, 0);
    sdv.replace(datum::nan, 0);
    vec cv = sdv/mu;
    cv.replace(datum::nan, 0);
    // obtain quantile
    vec cv_sort = cv;
    sort(cv_sort.begin(), cv_sort.end());
    int quantile_idx = cv_sort.size() * 0.25;
    uvec high_var_genes = find(mu >= 1 && cv >= cv_sort[quantile_idx]);
    // uvec high_var_genes = find(mu >= 1 && cv >= 0.02);
    if(high_var_genes.size() < 500){
        high_var_genes.resize(count.n_rows);
        high_var_genes = find(mu == mu);
    }
    mat count_hv = arma::zeros(high_var_genes.size(),count.n_cols);
    count_hv = count.rows(high_var_genes);
    return(count_hv);
}// end func

//*************************************************************************//
//                        Find neighbors for each cell                     //
//*************************************************************************//
List find_neighbors(mat count_hv,vec cell_labels){
    vec cluster = unique(cell_labels);
    sort(cluster.begin(), cluster.end());
    Rcpp::List neighbors(2);
    Rcpp::List dist_list(cluster.size());
    double npc = 0;
    Rcout << "calculating cell distances ..." << "\n";
    for(int iclust = 0; iclust < cluster.size(); ++iclust){
        uvec cell_inds = find(cell_labels == cluster[iclust]);
        mat count_hv_sub = arma::zeros(count_hv.n_rows,cell_inds.size());
        count_hv_sub = count_hv.cols(cell_inds);
        mat temp = count_hv_sub.each_col() - sum(count_hv_sub,1) / double(count_hv_sub.n_cols);
        mat U;
        vec s;
        mat V; //same as svd$v
        svd_econ(U,s,V,temp.t());
        vec eigs = s % s / double(count_hv_sub.n_cols - 1); //sample as pca$eigs
        vec var_cum = cumsum(eigs)/sum(eigs);
        if(max(var_cum) <= 0.4){ // var_thre = 0.4
            npc = var_cum.size();
        }else{
            uvec greater = find(var_cum>0.4);
            npc = greater(0);
        }
        if (npc < 3){ npc = 3; }
        mat pca_x = temp.t() * V;
        mat dist_cells_mat = calc_dis(pca_x.cols(0,npc));
        dist_list(iclust) = dist_cells_mat;
    }
    neighbors(0) = dist_list;
    neighbors(1) = cell_labels;
    return(neighbors);
}// end func

//*************************************************************************//
//                        Find valid genes for imputation                  //
//*************************************************************************//
uvec find_va_genes(mat paramt, mat subcount){
    uvec va_genes;
    double point = log10(1.01);
    uvec valid_genes1 = find(sum(subcount,1) > point * double(subcount.n_cols));
    uvec Nan_genes = find_nonfinite(paramt.col(0));
    uvec valid_genes = std_setdiff(valid_genes1,Nan_genes);
    if(valid_genes.size() == 0){
        va_genes = valid_genes;
    }else{
        // find out genes that violate assumption
        vec mu = paramt.col(3);
        vec dcheck1 = dgamma(mu+1, paramt.col(1), paramt.col(2));
        vec dcheck2 = dnorm(mu+1, paramt.col(3), paramt.col(4));
        uvec sgene = find((mu <= log10(1+1.01)) || (dcheck1 >= dcheck2 && mu <= 1));
        va_genes = std_setdiff(valid_genes,sgene);
    }
    return(va_genes);
}// end func

//*************************************************************************//
//                  ADMM algorithms to get the optimal weight              //
//*************************************************************************//
vec impute_ADMM(int cellid, mat subcount,
                uvec geneid_drop, uvec geneid_obs,
                uvec neighbors,double rho, int max_iter, double tol, double point){
    vec yimpute = arma::zeros(geneid_drop.size() + geneid_obs.size());
    if (geneid_drop.size() == 0 || geneid_drop.size() == subcount.n_cols) {
        yimpute =  subcount.col(cellid);
        
    }else{
    mat A = subcount(geneid_obs,neighbors); //observed count to provide information
    vec tmp = subcount.col(cellid);
    vec Y = tmp.elem(geneid_obs); // for the specific cell
    vec z = arma::zeros(neighbors.size());        // initial value of z
    vec weight = arma::zeros(neighbors.size());
    vec u = arma::zeros(neighbors.size()); // initial value of u
    mat Imat = eye<mat>(neighbors.size(),neighbors.size());
    mat Invterm = inv(rho * Imat + 2 * A.t() * A);
    mat AtY = A.t() * Y;
    for(int i = 0; i < max_iter; ++i) {
        vec weight_old = weight;
        vec z_old = z;
        weight = Invterm * (2*AtY + rho* (z - u)); // update weight
        z = weight + u;//update z
        z.elem(find(z < point)).fill(point);
        u = u + weight - z;
        vec r = weight - z;
        vec s = -rho * (z - z_old);
        if(norm(r,2) <= tol && norm(s,2) <= tol){
            break;
        }
    }
    mat ximpute = subcount(geneid_drop,neighbors);
    yimpute.elem(geneid_drop) = ximpute * weight;
    yimpute.elem(geneid_obs) = Y;
    vec maxobs = max(subcount, 1);
    yimpute.elem(find(yimpute > maxobs)) = maxobs.elem(find(yimpute > maxobs));
    }
    return(yimpute);
}// end func

//*************************************************************************//
//                        Get Dropout rate                                 //
//*************************************************************************//
mat get_droprate(mat subcount, mat paramt, double drop_thre){
     mat droprate = arma::zeros(subcount.n_rows,subcount.n_cols);
     for(int i = 0; i < subcount.n_rows; ++i){
         mat weightTemp = calculate_weight(subcount.row(i).t(),paramt.row(i).t());
         rowvec dropTemp = weightTemp.col(0).t();
         dropTemp.elem(find(dropTemp > drop_thre && subcount.row(i) > paramt(i,3))).fill(0);
         droprate.row(i) = dropTemp;
     }
     return(droprate);
}// end func

//*************************************************************************//
//                        Imputation model                                 //
//*************************************************************************//
//' A Rcpp function to perform imputation on scRNAseq data based on ADMM algorithm
//' @param count normalized count from raw data
//' @param cell_labels cell labels, it can be estimated or provided by user
//' @param point pseudo zero count
//' @param drop_thre dropout rate threshold
//' @param rho Step size in ADMM algprithm
//' @param max_iter max iteration in ADMM algorithm
//' @param tol tolerance
//' @return A list containing imputed data
//'
//' @export
// [[Rcpp::export]]
arma::mat imputation_wlabel_model(arma::mat count,arma::vec cell_labels, double point, double drop_thre,double rho,int max_iter,double tol){
    arma::mat count_hv = find_hv_genes(count);
    Rcout << "searching candidate neighbors ..." << "\n";
    List neighbors_res = find_neighbors(count_hv,cell_labels = cell_labels);
    List dist_list = neighbors_res(0);
    arma::vec clust = neighbors_res(1); 
    //mixture model
    arma::vec cluster = unique(clust);
    sort(cluster.begin(), cluster.end());
    arma::mat count_imp = count;
    for(int cc = 0; cc < cluster.size(); ++cc){
       Rcout << "estimating dropout probability for type " << cc + 1 << "..." << "\n";
        arma::uvec cell_inds = find(clust == cluster(cc));
        arma::mat paramt = get_mix_parameters(count.cols(cell_inds),point);
        Rcout << "cell size for type " << cc + 1 << " is " << cell_inds.size() << "\n";
        if(cell_inds.size() <= 1) continue;
        Rcout << "searching for valid genes ..." << "\n";
        arma::uvec valid_genes = find_va_genes(paramt, count.cols(cell_inds));
        if(valid_genes.size() <= 10) continue;
        arma::mat subcount = count.submat(valid_genes, cell_inds);
        arma::mat subparamt = paramt.rows(valid_genes);
        arma::mat droprate = get_droprate(subcount,subparamt,drop_thre);
        Rcout << "imputing dropout values for type " << cluster(cc) << "..." << "\n";
        arma::mat subcount_imp = arma::zeros(subcount.n_rows,subcount.n_cols);
        arma::uvec cell_inds_sub = sort_index(subcount.row(0).t());
        sort(cell_inds_sub.begin(),cell_inds_sub.end());
        for(int cellid = 0; cellid < subcount.n_cols; ++cellid){
            arma::uvec neighbors = find(cell_inds_sub != cellid);
            arma::uvec geneid_drop = find(droprate.col(cellid) > drop_thre);
            arma::uvec geneid_obs = find(droprate.col(cellid) <= drop_thre);
            subcount_imp.col(cellid) = impute_ADMM(cellid,subcount,geneid_drop,geneid_obs,neighbors,rho,max_iter,tol,point);
        }
        count_imp.submat(valid_genes,cell_inds) = subcount_imp;
    }
    count_imp.elem(find(count_imp < point)).fill(point);
    return(count_imp);
}// end func


