## ============================================================================
## Test the tree structure of cell type transcriptomes
## Main function name: treeness_test
## Dependency: pvalue (bioconductor)
## Preperation functions: treeness_pval; treeness_delta; 
## Input: data frame for cell type specific transcriptomes
## Output: pi_0, histogram of pvals, reconstructed cell type tree, phylotree
## ============================================================================

## dependency: function to get p-value for tree structure
treeness_pval <- function(delta){
  
  treeness.pval <- 3/pi*(atan( (2*delta - 1 )/ sqrt(3)) + pi/6)
  
  return(treeness.pval)
}

## dependency: function to get detla value for a tetrad
# tetrad: genes * samples
treeness_delta <- function(tetrad){
  
  tetrad.dist <- as.matrix( dist(t(tetrad)) )
  H1 <- tetrad.dist[1,2] + tetrad.dist[3,4]
  H2 <- tetrad.dist[1,3] + tetrad.dist[2,4]
  H3 <- tetrad.dist[1,4] + tetrad.dist[2,3]
  H.vec <- sort(c( H1, H2, H3) )
  
  delta <- (H.vec[3] - H.vec[2]) / (H.vec[3] - H.vec[1])
  
  return(delta)
}

## Main function
treeness_test <- function(my.df){
  
  n.sample <- dim(my.df)[2]
  v.combn <- combn(n.sample, 4)
  
  pvals <- c()
  for( i in 1:dim(v.combn)[2] ){
    delta <- treeness_delta( my.df[ , v.combn[, i] ] )
    pvals <- c( pvals, treeness_pval( delta ) )
  }
  
  res <- qvalue(pvals)
  
  return(res)
}
