# S2Bapp

apid_main.RData contains an igraph object of the protein-protein interaction network extracted from APID in May 2019. The network contains only the main component and is simplified (no self edges or multiple edges).

genbet.RData contains the general betweenness values for all the nodes in the network.

symbnodes.txt contains the gene symbols of the network nodes (which are identified by uniprot id in apid_main.RData).

ALS_list4.txt and SMA_list4.txt are example input files.

app.R contains the code for the shiny app.

