##########################################################################################################
## Correlated evolution test for two tissue types
## Main function name: corrEvo; perGene_corrEvo
## Dependency: ape
## Function dependency: get_contrast
## Input: data frame for gene expression levels in two tissue types of multiple species, known species phylogeny
## Output: level of correlated evolution for respective tissue type
##########################################################################################################

## Note: Remove mitochondrial genes before using continuous model
## Use the same function for correlated evolution for discritized, Brownian, as well as OU model

## get contrast vectors
get_contrast <- function(tissue.1, tissue.2, my.phylo, ...){
  if(  sum(colnames(tissue.1) != colnames(tissue.2)) + sum(!(colnames(tissue.1) %in% my.phylo$tip.label) ) ){
    cat('Species in data frames must match each other and species in phylogeny! \n')
    cat('Please check input...\n')
    stop(-1)
  }
  
  n.gene <- dim(tissue.1)[1]
  Contrast1 <- c()
  Contrast2 <- c()
  for(i in 1:n.gene){
    trait1 = tissue.1[i, ]
    names(trait1) = colnames(tissue.1)
    trait2 = tissue.2[i, ]
    names(trait2) = colnames(tissue.2)
    
    # Note, get contrast for each gene.
    Contrast1 <- rbind(Contrast1, pic(trait1, my.phylo, ...)) # scaled = T
    Contrast2 <- rbind(Contrast2, pic(trait2, my.phylo, ...))
  }
  
  return(list(Contrast1 = Contrast1, 
              Contrast2 = Contrast2))
}

## Get level of correlated evolution (LCE)
eval_corrEvo <- function(tissue.1, tissue.2, my.phylo, ...){
  
  if(  sum(colnames(tissue.1) != colnames(tissue.2)) + sum(!(colnames(tissue.1) %in% my.phylo$tip.label) ) ){
    cat('Species in data frames must match each other and species in phylogeny! \n')
    cat('Please check input...\n')
    stop(-1)
  }
  
  species <- colnames(tissue.1)
  n.gene <- dim(tissue.1)[1]
 
  ## The most recent common ancestor
  my.mrca <- apply(mrca(my.phylo), 1, as.character)
  ## Branching time
  my.bctime <- branching.times(my.phylo)

  ## Independent contrast for all genes in two tissues
  contrasts <- get_contrast(tissue.1, tissue.2, my.phylo, ...)
  corrEvo <- diag( cor( contrasts$Contrast1, contrasts$Contrast2 ) )
  
  return(list(dist = my.bctime,
              corrEvo = corrEvo))
}

## per-gene correlated evolution (per-gene LCE)
eval_perGene_corrEvo <- function(tissue.1, tissue.2, my.phylo, ...){
  if(  sum(colnames(tissue.1) != colnames(tissue.2)) + sum(!(colnames(tissue.1) %in% my.phylo$tip.label) ) ){
    cat('Species in data frames must match each other and species in phylogeny! \n')
    cat('Please check input...\n')
    stop(-1)
  }
  
  if(length(my.phylo$tip.label) < 6){
    cat('Number of species too low to estimate per-gene correlated evolution!\n')
    cat('Exiting...\n')
    stop(-1)
  }
  
  species <- colnames(tissue.1)
  n.gene <- dim(tissue.1)[1]
  
  ## The most recent common ancestor
  my.mrca <- apply(mrca(my.phylo), 1, as.character)
  ## Branching time
  my.bctime <- branching.times(my.phylo)
  
  ## Independent contrast for all genes in two tissues
  contrasts <- get_contrast(tissue.1, tissue.2, my.phylo, ...)
  
  LCE <- c()
  for(i in 1:n.gene){
    LCE <- c( LCE, cor(contrasts$Contrast1[i, ], contrasts$Contrast2[i,]) )
  }
  
  names(LCE) <- rownames( tissue.1 )
  
  return(LCE)
}
