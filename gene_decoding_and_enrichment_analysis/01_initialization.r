# This analysis allows to test significant overalpping between two list of genes.
# You are required to have the first list form previous analysis, e.g. PPI interactome
# This code allows you to create the second list of genes. These genes are those that are 
# spatially correlated with imaging maps (e.g. between group t-contrast rfMRI maps), after
# FDR correction. This part of the analysis is done with a code that uses Neurovault, and 
# test the beta values of the gene-imaging correlations of each donor against 0. 
# To tun this GLM analysis make sure you uploaded your MRI map in NeuroVault 
# (https://neurovault.org/my_collections/?q=), "image" takes the ID value og NeuroVault
# Statistical testing of overlapping between the two list genes is then performed with 
# hypergeometric function.  

# before running the analysis set these paths:

source("/media/DATA1/IIT_17Sep_Tsc2/05_gene_decoding/05b_initial_analysis_with_neurovault_PAPER/gene_decode.r") # edit this, /path/02_gene_decode_R.r

dataDir <- "/media/DATA1/IIT_17Sep_Tsc2/05_gene_decoding/05b_initial_analysis_with_neurovault_PAPER/" # edit this, your top output directory
measure <- "results"
image <- "136272" # the image identifier on NeuroVault (id)

res = gene_decode(dataDir, image, measure, dbs, term, maxterms, prefix, dbs)


