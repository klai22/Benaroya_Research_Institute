#Summer Intern | Benaroya Research Institute (Bioinformatics Core)                                                    06/2023-08/2023
#Researched antigen specificity and RNA expression differences between regulatory and conventional T cells in type 1 diabetes (T1D) patients. 
#Developed custom R pipelines to analyze scRNAseq data from 10 healthy donors and 12 T1D donors, integrating metrics from clinical assessments and public databases.
#Identified unique TCR motifs by examining network-based clusters and correlated T1D with genes exhibiting low expression. 
#Delivered presentation on islet antigen-specific Treg transcriptomic profiles and their clinical applications in T1D treatment at BRI’s summer symposium. 

#CDR3vsVDJDB.R : Processing and Visualizing CDR3 (TCR) sequences and their associated metadata; filtered by the sequences that match public VDJDB(database) to infer potential epitope cross-reactivity
#TCRdistheatmap.R : processing and visualizing TCRDIST distance matrix results along with their associated experimental/metadata (scRNA data) 
#DEA.R : Processing and visualizing differential expression analysis results from scRNA data along with assocaited metadata 
#DEA(geneset).R : repeating analysis from DEA.R, but using a supervised approach to differential expression analysis so that I could create a heatmap only displaying expression in a specific set of genes unique to certain cell types.  
