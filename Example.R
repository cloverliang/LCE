##########################################################################################################
# Example part I: Normalize TPMs from multiple species
# Input: data frames of TPM values for multiple species and tissue types, ortholog gene list
# Output: A single data frame with normalized TPMs
##########################################################################################################

## Alligator & chicken digit data
load('../data/TPM_am.RData')
load('../data/TPM_gg.RData')

## alligator and chicken ortholog
load('../data/gg_am_ortholog.RData')

## cut ortholog genes
tpm.am.otlg <- tpm.am[ match( ortholog$am_gene_id, rownames(tpm.am) ), ]
tpm.gg.otlg <- tpm.gg[ match( ortholog$gg_gene_id, rownames(tpm.gg) ), ]

## normalize alligator TPMs to chicken TPMs
tpm.am.otlg <- sweep( x= tpm.am.otlg, STATS = colSums(tpm.am.otlg), MARGIN = 2, FUN = '/') * 10^6

## merge data frame
tpm.gg.am <- cbind( ortholog, tpm.gg.otlg, tpm.am.otlg)

##########################################################################################################
## Example part II: treeness test
## Input: Data frame with normalized TPM for alligator digits
## Output: p-values for all tetrads
##########################################################################################################

## prep
library(qvalue)
source('src/TreenessTest.R')

## alligator digit data
load('../data/TPM_am.RData')

## treeness test and distribution of p-values
res <- treeness_test(tpm.am)
res$pi0
pdf( file = '../results/p-value_distribution.pdf', height = 4, width = 4)
hist( res$pvalues, breaks = 30, xlab = 'p-values', 
      freq = F, main = 'Histogram of p-values' )
abline(h = res$pi0, col = 'red', lty = 2 )
dev.off()

##########################################################################################################
## Example part III: Estimation of level of correlated evolution (LCE)
## Input: Gene expression TPMs from Brawand 2011
## Output: estimated LCE
## Dataset reference: Brawand D, et al. 2011. The evolution of gene expression levels in mam- malian organs. Nature 478(7369):343–348.
##########################################################################################################

## prep
library(ape)
library(ggplot2)
source('src/CorrelatedEvolution.R')

## Brawand dataset
load('../data/Brawand_data.RData')
load('../data/Brawand_phylo.RData')

amniote.names <- amniote.data[, 1:9]
amniote.data <- amniote.data[, 10:141]

amniote.tissue = factor(sapply(strsplit(names(amniote.data), '[.]'), '[[', 2))
amniote.species = factor(sapply(strsplit(names(amniote.data), '[.]'), '[[', 1))

## correlated evolution between brain and heart
brain = amniote.data[, amniote.tissue == 'br']
heart = amniote.data[, amniote.tissue == 'ht']
brain.agg <- aggregate( t(brain), list(species = sapply(strsplit(colnames(brain), split = '[.]'), '[[', 1)) , FUN = mean)
heart.agg <- aggregate( t(heart), list(species = sapply(strsplit(colnames(heart), split = '[.]'), '[[', 1)) , FUN = mean)
brain.df <- data.frame( t( brain.agg[, 2:dim(brain.agg)[2]] ) )
heart.df <- data.frame( t(heart.agg[, 2:dim(heart.agg)[2]]) )
colnames(brain.df) <- brain.agg[,1]
colnames(heart.df) <- heart.agg[,1]

res <- eval_corrEvo(sqrt(brain.df), sqrt(heart.df), brawand.phylo)
LCE.df <- data.frame(SplitTime = res$dist, LCE = res$corrEvo)
pdf('../results/LCE_brain_heart.pdf', height = 4, width = 4)
p.LCE <- ggplot(LCE.df, aes(x = SplitTime, y = LCE)) +
            geom_point() + 
            theme_classic() + 
            ggtitle('Brain vs heart')
plot(p.LCE)
dev.off()