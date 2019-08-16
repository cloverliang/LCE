##########################################################################################################
## Get the phylogeny for a set of species
## Main function name: fetch_phylogeny
## Dependency: ape
## Input: species latin names
## Output: my.phy
## Database version: TimetreeOfLife2015
## Database: http://www.timetree.org/
## Database citation:Kumar S, Stecher G, Suleski M, Hedges SB (2017) TimeTree: A Resource for Timelines, Timetrees, and Divergence Times. Mol Biol Evol doi:10.1093/molbev/msx116
##########################################################################################################

fetch_phylogeny <- function(species = c('Felis_catus', 'Bos_taurus', 'Canis_lupus_familiaris', 'Cavia_porcellus', 
                                      'Monodelphis_domestica', 'Sus_scrofa', 'Homo_sapiens', 'Oryctolagus_cuniculus', 
                                      'Rattus_norvegicus', 'Ovis_aries', 'Mus_musculus', 'Equus_caballus'),
                           tree.file = 'data/TimetreeOfLife2015.nwk'){

  mammalia_dosReis <- read.tree(file = tree.file)
  
  cat( paste( setdiff(species, mammalia_dosReis$tip.label), 
              'are not found in the database. \n Please double check the species latin name...\n', 
              sep = ' ') ) 
  
  my.phy <- drop.tip(mammalia_dosReis, tip = which( !(mammalia_dosReis$tip.label %in% species) ))
  my.phy$node.label = (my.phy$Nnode + 2):(2*my.phy$Nnode + 1)
  
  return(my.phy)
}
