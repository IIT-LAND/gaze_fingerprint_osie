library(easypackages)
libraries("here","ggplot2","ggridges","reshape2","patchwork","RColorBrewer","plotrix","tidyverse","cluster","ggpackets")

codepath = here("code")
resultpath = here("results")

# make geom_scatterbox for later plotting
geom_scatterbox <- ggpacket() +
  geom_jitter() +
  geom_boxplot(fill = NA, colour = "#000000", outlier.shape = NA)

# make linear model plots 

geom_linear <- ggpacket() + geom_point() + geom_smooth(method = lmrob)


hamming_bin_dist <- function(df4analysis){

  df4analysis = as.matrix(df4analysis)
  res_mat = data.frame(matrix(nrow = dim(df4analysis)[1], ncol = dim(df4analysis)[1]))

  for (irow in 1:nrow(res_mat)){
    for (icol in 1:ncol(res_mat)){
      # print(sprintf("row %d, col %d",irow,icol))
      if (irow==icol){
        res_mat[irow,icol] = 0
      } else{
        mask1 = df4analysis[irow,]==1
        mask2 = df4analysis[icol,]==1
        percentage_overlap = (dim(df4analysis)[2] - sum(mask1 & mask2))/dim(df4analysis)[2]
        res_mat[irow,icol] = percentage_overlap
        res_mat[icol,irow] = percentage_overlap
      } # if (irow==icol)
      
    } # for icol
  } # for irow
  
  return(res_mat)
} # function hamming_bin_dist



# function to use silhouette to determine the best k
best_k_nbclust <- function(df2use, 
                   clust_method = "ward.D2",
                   dist_method = "euclidean", 
                   index = "cindex",
                   min.nc = 2,
                   max.nc = 15){
  
  require(NbClust)
  
  res = NbClust(data = df2use, 
                method = clust_method, 
                distance = dist_method, 
                index = index,
                min.nc = min.nc,
                max.nc = max.nc)
  optimal_k = res$Best.nc[["Number_clusters"]]
  
  return(optimal_k)
} # function best_k_nbclust



# function to use silhouette to determine the best k
best_k <- function(df2use, 
                   k_range = c(2:15), 
                   dist_method = "euclidean", 
                   clust_method = "ward.D2"){
  
  colnames2use = c("k","score")
  silhouette_scores = data.frame(matrix(nrow = length(k_range), ncol = length(colnames2use)))
  colnames(silhouette_scores) = colnames2use
  rownames(silhouette_scores) = k_range
  
  if (dist_method=="hamming") {
    dist_mat = as.dist(hamming_bin_dist(df2use))
  }else{
    dist_mat =  dist(df2use, method = dist_method)
  }
  cres = hclust(dist_mat, method = clust_method)

  for (i in 1:length(k_range)){
    
    k = k_range[i]
    cut_ward_rows =  cutree(cres, k = k)
    
    silhouette_res = silhouette(cut_ward_rows, dist = dist_mat)
    silhouette_scores[i,"k"] = k
    silhouette_scores[i,"score"] = mean(silhouette_res[,"sil_width"])
    
  }
  
  optimal_k = silhouette_scores[silhouette_scores[,"score"]==max(silhouette_scores[,"score"]),"k"]
  return(optimal_k)
} # function best_k


# function to do enrichment tests
enrichmentTest <- function(nFP_given_Semantic, nFP, nSemantic, nTotal){
  
  A = nFP_given_Semantic
  B = nFP - nFP_given_Semantic
  C = nSemantic - nFP_given_Semantic
  D = nTotal - (A+B+C)
  
  if (A==0){
    A = A+1
  }
  
  if (B==0){
    B = B+1
  }
  
  if (C==0){
    C = C+1
  }
  
  if (D==0){
    D = D+1
  }
  
  contingency_table = as.table(t(matrix(c(A,B,C,D), nrow=2, ncol = 2)))
  res = fisher.test(contingency_table, alternative="greater")
  
  colnames2use = c("OR","p","lowCI","highCI")
  res_df = data.frame(matrix(nrow = 1, ncol = length(colnames2use)))
  colnames(res_df) = colnames2use
  
  res_df[,"OR"] = res$estimate
  res_df[,"p"] = res$p.value
  res_df[,"lowCI"] = res$conf.int[1]
  res_df[,"highCI"] = res$conf.int[2]
  
  return(res_df)
} # function enrichmentTest



lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}



cohens_d <- function(x, y, DIM=1, SIGN=TRUE, na.rm=TRUE) {
  #
  # Function will compute cohen's d effect size.
  # Generalized to work on either matrices, data frames, vectors
  #
  # INPUT
  #	x <- matrix or numeric vector
  #	y <- matrix or numeric vector
  #	DIM <- specify the dimension which samples are along
  #
  # Example usage
  #
  # x <- cbind(rnorm(100,0,1), rnorm(100,0,1), rnorm(100,0,1))
  # y <- cbind(rnorm(100,1,1), rnorm(100,2,1), rnorm(100,3,1))
  # d <- cohens_d(x, y, 1)
  #
  # written by mvlombardo - 28.08.2015
  #
  
  library(matrixStats)
  
  # if x and y are vectors, coerce them into matrices
  if (class(x)=="numeric" | class(x)=="integer") {
    x <- as.matrix(x)
  } # if
  
  if (class(y)=="numeric" | class(y)=="integer") {
    y <- as.matrix(y)
  }# if
  
  if (na.rm==TRUE){
    missingValDecision = TRUE
  } else {
    missingValDecision = FALSE
  }
  
  # n-1 for x and y
  lx <- dim(x)[DIM]-1
  ly <- dim(y)[DIM]-1
  
  # if samples are along the rows
  if (DIM==1){
    if (SIGN){
      # mean difference (numerator)
      md <- colMeans(x, na.rm = missingValDecision) - colMeans(y, na.rm = missingValDecision)
    } else{
      # mean difference (numerator)
      md <- abs(colMeans(x, na.rm = missingValDecision) - colMeans(y, na.rm = missingValDecision))
    }# if (SIGN)
    # pooled variance (denominator), but before any sqrt is done
    csd <- (lx * rowVars(t(x),na.rm = missingValDecision)) + (ly * rowVars(t(y), na.rm = missingValDecision))
    
    # else if samples are along the columns
  } else if (DIM==2){
    if (SIGN){
      # mean difference (numerator)
      md <- rowMeans(x, na.rm = missingValDecision) - rowMeans(y, na.rm = missingValDecision)
    } else{
      # mean difference (numerator)
      md <- abs(rowMeans(x, na.rm = missingValDecision) - rowMeans(y, na.rm = missingValDecision))
    }# if (SIGN)
    # pooled variance (denominator), but before any sqrt is done
    csd <- lx * rowVars(x, na.rm = missingValDecision) + ly * rowVars(y, na.rm = missingValDecision)
  }# end if
  
  # divide pooled variance by sum of n-1 for x and y and then square root it
  csd <- sqrt(csd/(lx + ly))
  # compute cohen's d
  cd  <- md/csd
}# end cohens_d <- function(x, y, DIM)

variance_d <- function(d, sample1_n, sample2_n){
  
  #function computes the variance for a given effect size, as calculated by cohen's d 
  #requires three variables: 
  # 1) cohen's d effect size 
  # 2) sample size of first group you are comparing, 
  # 3) sample size of second group you are comparing
  
  #variance is required to statistically compare two effect sizes!
  
  #formula kindly found here:
  #https://stats.stackexchange.com/questions/77269/statistical-comparison-of-2-independent-cohens-ds
  #https://stats.stackexchange.com/questions/144084/variance-of-cohens-d-statistic
  
  #author: S K Crockford 14.06.2023 (sarahcrockford@mac.com)
  
  n1 = 1/(sample1_n)
  n2 = 1/(sample2_n)
  
  d_sqr = d**2
  
  denominator = 2*(sample1_n + sample2_n)
  
  v <- n1 + n2 + (d_sqr/denominator)
  
  return(v)
  
} # end variance_d <- function(cohens_d, n1, n2)

cohendiff <- function(d1, d2, v1, v2){
  
  #function computes the statistical significance between two effect sizes, as calculated by cohen's d 
  #requires four variables: 
  # 1) first cohen's d effect size you are comparing
  # 2) second cohen's d effect size you are comparing
  # 3) variance of first cohen's d effect size you are comparing
  # 4) variance of second cohen's d effect size you are comparing
  
  #variance is required to statistically compare two effect sizes!
  #we will use the R function pnorm to get the CDF/P-value of the results z-statistic
  
  #formula kindly found here:
  #https://stats.stackexchange.com/questions/77269/statistical-comparison-of-2-independent-cohens-ds
  #Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2021). Introduction to meta-analysis. John Wiley & Sons.
  #formula found in Chapter 19: Subgroup Analyses, (pg. 156 of the 2009 version)
  
  #author: S K Crockford 14.06.2023 (sarahcrockford@mac.com)
  
  H0 = "null hypothesis: cohen's d effect size is same for the two groups tested"
  
  cohen_diff = d1 - d2
  root_variance = sqrt(v1 + v2)
  
  z = cohen_diff/root_variance
  
  cdf = pnorm(q = abs(z), lower.tail = TRUE)
  
  p = 2*(1 - cdf)
  
  if (p < 0.05){
    statement = "Cohen's d effect sizes are significantly different at alpha = 0.05, null hypothesis is rejected" }
  else{
    statement = "Cohen's d effect sizes are not significantly different at alpha = 0.05, null hypothesis is retained"
  }
  
  results_list <- list("null hypothesis" = H0, "z_score" = z, "p_value" = p, "result" = statement)
  
  return(results_list)
  
} # end cohendiff <- function(d1, d2, v1, v2)

# get_ggColorHue.R
# 
# Function will spit out color codes for a given n, just like ggplot2. 
# Helpful when you want to manually set colors in plots.
#

get_ggColorHue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
} # function get_ggColorHue