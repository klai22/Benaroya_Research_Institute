#setwd
setwd("C:/Users/klai/Box/2023SummerInternKenneth/data")

#load packages 
library("ggpubr")
library("gplots")
library("ggplot2")
library("dplyr")
library("tidyr")
library("lubridate")
library("tibble")
library("igraph")
library("ComplexHeatmap")
library("data.table")
library("RColorBrewer")
library("pheatmap")
library("gridExtra")
library("grid")
library("scales")
library("circlize") # has colorRamp2 

#IMPORTDATA 
  
  #Unzipping Expanded Differential Expression data 
    #unzip("T1DvHC_Expanded_DifferentialExpression.zip", exdir = "/Users/klai/Box/2023SummerInternKenneth/data")

  #loading the 3 sub-data files post-zipping 
    #t1dExpanded_expression.RDS: TMM-normalized, filtered, gene expression matrix on pseudobulked expanded cells
        #rows correspond to genes, and columns correspond to individual samples (sub-samples separated by cellType & studyGroup). Each element of the matrix contains the gene expression level for a specific gene in a specific sample.
      EXPexpressionmatrixdata<- readRDS("t1dExpanded_expression.RDS")
    #t1dExpanded_metadata.RDS: metadata data frame for the above
      EXPmetadata<- readRDS("t1dExpanded_metadata.RDS")
    #t1dExpandedDE.RDS: a list of dataframes showing results from limma differential expression analysis between T1D and HC 
      EXPexpressionresultsdata<- readRDS("t1dExpandedDE.RDS")
      
    #turning expression results into separate dataframes for cell types 
      EXPTregexpressionresults = EXPexpressionresultsdata[[1]]
      EXPTconvsexpressionresults = EXPexpressionresultsdata[[2]]

#RENAMING of matrices so that gene ensemble IDs match the HGNC.symbol (actual gene name) from metadata 
      #making a df of all the ensemble IDs & their matching HGNC symbols 
      Treg_symbols = as.data.frame(EXPTregexpressionresults$HGNC.symbol, rownames(EXPTregexpressionresults))
      Tconvs_symbols = as.data.frame(EXPTconvsexpressionresults$HGNC.symbol, rownames(EXPTconvsexpressionresults))
      
      Treg_symbols$symbol <- Treg_symbols$`EXPTregexpressionresults$HGNC.symbol`
      Treg_symbols$`EXPTregexpressionresults$HGNC.symbol` = NULL
      Tconvs_symbols$symbol <- Tconvs_symbols$`EXPTconvsexpressionresults$HGNC.symbol`
      Tconvs_symbols$`EXPTconvsexpressionresults$HGNC.symbol` = NULL
      
      #making a df for both of them combined (for the master plot) (found no duplicates)
      genesall_symbols <- rbind(Treg_symbols, Tconvs_symbols) 
      
      
      #Renaming matrix + results 
      matching_indices_matrix <- match(rownames(EXPexpressionmatrixdata), rownames(genesall_symbols))
      rownames(EXPexpressionmatrixdata) <- genesall_symbols$symbol[matching_indices_matrix]
      
      #matching_indices_matrix_Treg <- match(rownames(SIGEXPexpressionmatrix_reordered_Treg), rownames(SIGgenesall_symbols))
      #rownames(SIGEXPexpressionmatrix_reordered_Treg) <- SIGgenesall_symbols$symbol[matching_indices_matrix_Treg]
      
      #matching_indices_matrix_Tconvs <- match(rownames(SIGEXPexpressionmatrix_reordered_Tconvs), rownames(SIGgenesall_symbols))
      #rownames(SIGEXPexpressionmatrix_reordered_Tconvs) <- SIGgenesall_symbols$symbol[matching_indices_matrix_Tconvs]
      
      
#Isolating only rows that match Treg/Tconvs SIGNATURE GENES 
      #Creating list of gene symbols that belong to the signature genes ("gene set") for each unique cell type, these were sent to me by Alex 
        Tregsignaturegeneset <- c("TNFRSF1B","HPGD","ZNF532","STAM","LRRC32","IL1R2","ZBTB38","TRIB1","ICA1","SGMS1","IKZF2","TNFRSF9","CTLA4","FOXP3","VAV3","CSF2RB","METTL7A","IL1R1","ZC2HC1A")
      
        Thsignaturegeneset <- c("DACT1","ABCB1","PTPRK","ANK3","NELL2","HDGFL3","IL7R","ID2","LPIN2")
      
      #Creating master list for now
        #allgenesets <- c (Tregsignaturegeneset, Thsignaturegeneset)
      
      #Isolating rows of matrix that match this list (only ones that match the signature genes) [rows went from 12283-->26, 0.21% of genes in matrix were found to match all gene sets combined]
      
      #genesetmatched_expressionmatrixdata<- EXPexpressionmatrixdata[rownames(EXPexpressionmatrixdata) %in% allgenesets, ]
        #Treg isolation: [rows went from 12283--> 19 ] 
        Treggenesetmatched_expressionmatrixdata<- EXPexpressionmatrixdata[rownames(EXPexpressionmatrixdata) %in% Tregsignaturegeneset, ]
        #Th isolation: [rows went from 12283--> 7 ] 
        Thgenesetmatched_expressionmatrixdata<- EXPexpressionmatrixdata[rownames(EXPexpressionmatrixdata) %in% Thsignaturegeneset, ]
        
#FILTER DATA: PROBLEM, THERE WERE NO MATCHES BETWEEN (1) the genesets provided & (2) the list of genes which met the FDR requirement. 
  # I want to filter out only the sig. genes (expression values in matrix) with FDR-differential-expression-values [adj. p-value] that were =>5% (0.05) [FDR values in expressionresults]
      #Make a list of the genes that meet his req. 
        #SIGgenesTreg = subset(EXPTregexpressionresults,adj.P.Val <= 0.05 )
        #SIGgenesTreglist = as.list(rownames(SIGgenesTreg))
      
        #SIGgenesTconvs = subset(EXPTconvsexpressionresults,adj.P.Val <= 0.05 )
        #SIGgenesTconvslist = as.list(rownames(SIGgenesTconvs))
      
      #Make a master list of these genes (eliminating duplicates) [there were no matches between these 2 lists] 
        #allSIGgeneslist = union(SIGgenesTreglist, SIGgenesTconvslist)
      
      #subset the matrixes so that it only includes the genes that meet the FDR-differential-expression-threshold (0.05) (isolate rows w/ rownames match the list of sig genes only)
        #SIGEXPexpressionmatrix_Treggeneset = Treggenesetmatched_expressionmatrixdata[rownames(Treggenesetmatched_expressionmatrixdata) %in% allSIGgeneslist, ]
        #SIGEXPexpressionmatrix_Thgeneset = Thgenesetmatched_expressionmatrixdata[rownames(Thgenesetmatched_expressionmatrixdata) %in% allSIGgeneslist, ]
      
      
#ORGANIZING DATA / PREPPING FOR HEATMAP ANNOTATIONS (ordering matrix and pre-annotations data properly)
  #Re-ordering / grouping matrix columns based on control group (beginning of column name string) and celltype (ending of column name string)
      #converting (SIGEXPexpressionmatrix) from matrix --> df because otherwise the code below won't work 
        Treggenesetmatched_expressionmatrixdata = as.data.frame(Treggenesetmatched_expressionmatrixdata)
        Thgenesetmatched_expressionmatrixdata = as.data.frame(Thgenesetmatched_expressionmatrixdata)
      
      # Columns that start with "T1D" and end with "Treg"
      t1dtreg_cols_treggeneset <- grep("^T1D.*Treg$", names(Treggenesetmatched_expressionmatrixdata))
      t1dtreg_cols_thgeneset <- grep("^T1D.*Treg$", names(Thgenesetmatched_expressionmatrixdata))
      
      # Columns that start with "T1D" and end with "Tconventional"
      t1dtconvs_cols_treggeneset <- grep("^T1D.*Tconventional$", names(Treggenesetmatched_expressionmatrixdata))
      t1dtconvs_cols_thgeneset <- grep("^T1D.*Tconventional$", names(Thgenesetmatched_expressionmatrixdata))
      
      # Columns that start with "Control" and end with "Treg"
      ctrltreg_cols_treggeneset <- grep("^Control.*Treg$", names(Treggenesetmatched_expressionmatrixdata))
      ctrltreg_cols_thgeneset <- grep("^Control.*Treg$", names(Thgenesetmatched_expressionmatrixdata))
      
      # Columns that start with "Control" and end with "Tconventional"
      ctrltconvs_cols_treggeneset <- grep("^Control.*Tconventional$", names(Treggenesetmatched_expressionmatrixdata))
      ctrltconvs_cols_thgeneset <- grep("^Control.*Tconventional$", names(Thgenesetmatched_expressionmatrixdata))
      
      # Order the columns based on the grouping
      samplegroups_treggeneset <- c(t1dtreg_cols_treggeneset,t1dtconvs_cols_treggeneset,ctrltreg_cols_treggeneset,ctrltconvs_cols_treggeneset)
      samplegroups_thgeneset <- c(t1dtreg_cols_thgeneset,t1dtconvs_cols_thgeneset,ctrltreg_cols_thgeneset,ctrltconvs_cols_thgeneset)
      
      # Subset the dataframe with the ordered columns
      Treggeneset_expressionmatrix_reordered <- Treggenesetmatched_expressionmatrixdata[, samplegroups_treggeneset]
      Thgeneset_expressionmatrix_reordered <- Thgenesetmatched_expressionmatrixdata[, samplegroups_thgeneset]
      
      #NOTE: right now, SIGEXPexpressionmatrix_reordered is a df, if encounter problems need to turn df-->matrix 
      
  #Re-ordering metadata to match the matrix (the primary data represented in heatmap )
      
      #Isolate only the rows of metadata that are represented in the matrix 
      sampleslist_treggeneset = as.list(colnames(Treggeneset_expressionmatrix_reordered))
      sampleslist_thgeneset = as.list(colnames(Thgeneset_expressionmatrix_reordered))
      #subset the metadata so that it only includes the samples that match matrix (it seemed that EXP metadata already only had the significant sample names?, not sure if this was a coincidence or it already had been subsetted for signifigant metadata to begin with)
      metadata_treggeneset = EXPmetadata[rownames(EXPmetadata) %in% sampleslist_treggeneset, ]
      metadata_thgeneset = EXPmetadata[rownames(EXPmetadata) %in% sampleslist_thgeneset, ]
      
      #Reordering the metadata rows so that it matches the same order as the matrix 
      metadata_treggeneset_reordered <- metadata_treggeneset[match(sampleslist_treggeneset, rownames(metadata_treggeneset)), ]
      metadata_thgeneset_reordered <- metadata_thgeneset[match(sampleslist_thgeneset, rownames(metadata_thgeneset)), ]
      
#SUBSETTING MATRICES + METADATA FOR CELLTYPE IND. PLOTS 
      #Creating separate matrices for cell-type plots 
        #list of Treg & Tconvs samples 
        Tregsamples_treggeneset = as.list(rownames(metadata_treggeneset_reordered[metadata_treggeneset_reordered$type == "Treg", ]))
        Tconvssamples_treggeneset = as.list(rownames(metadata_treggeneset_reordered[metadata_treggeneset_reordered$type == "Tconventional", ]))
        
        Tregsamples_thgeneset = as.list(rownames(metadata_thgeneset_reordered[metadata_thgeneset_reordered$type == "Treg", ]))
        Tconvssamples_thgeneset = as.list(rownames(metadata_thgeneset_reordered[metadata_thgeneset_reordered$type == "Tconventional", ]))
      
        #create subset matrices for each cell type 
       Treggenesetmatrix_tregs <- Treggeneset_expressionmatrix_reordered[, names(Treggeneset_expressionmatrix_reordered) %in% Tregsamples_treggeneset]
       Treggenesetmatrix_tconvs <- Treggeneset_expressionmatrix_reordered[, names(Treggeneset_expressionmatrix_reordered) %in% Tconvssamples_treggeneset]
       
       Thgenesetmatrix_tregs <- Thgeneset_expressionmatrix_reordered[, names(Thgeneset_expressionmatrix_reordered) %in% Tregsamples_thgeneset]
       Thgenesetmatrix_tconvs <- Thgeneset_expressionmatrix_reordered[, names(Thgeneset_expressionmatrix_reordered) %in% Tconvssamples_thgeneset]
       
       #Create subset metadata for each celltype 
       treggeneset_metadata_Treg = subset(metadata_treggeneset_reordered, type == "Treg")
       treggeneset_metadata_Tconvs = subset(metadata_treggeneset_reordered, type == "Tconventional")
       
       thgeneset_metadata_Treg = subset(metadata_thgeneset_reordered, type == "Treg")
       thgeneset_metadata_Tconvs = subset(metadata_thgeneset_reordered, type == "Tconventional")
        

#CREATING HEATMAP ANNOTATIONS 
      
       #Setting color schemes 
      celltypecolors <- c("Treg" = "blue", "Tconventional" = "yellow") # defining the colors for the annotation
      studygroupcolors <- c("Control" = "aquamarine2", "T1D" = "palevioletred2") # defining the colors for the annotation
      #creating annotation obj.s 
      masterAnnot_treggeneset = HeatmapAnnotation( df=metadata_treggeneset_reordered[,c("studyGroup","type")], col=list(studyGroup=studygroupcolors, type=celltypecolors))
      TregAnnot_treggeneset = HeatmapAnnotation( df=treggeneset_metadata_Treg[,c("studyGroup","type")], col=list(studyGroup=studygroupcolors, type=celltypecolors))
      TconvsAnnot_treggeneset = HeatmapAnnotation( df=treggeneset_metadata_Tconvs[,c("studyGroup","type")], col=list(studyGroup=studygroupcolors, type=celltypecolors))
      
      masterAnnot_thgeneset = HeatmapAnnotation( df=metadata_thgeneset_reordered[,c("studyGroup","type")], col=list(studyGroup=studygroupcolors, type=celltypecolors))
      TregAnnot_thgeneset = HeatmapAnnotation( df=thgeneset_metadata_Treg[,c("studyGroup","type")], col=list(studyGroup=studygroupcolors, type=celltypecolors))
      TconvsAnnot_thgeneset = HeatmapAnnotation( df=thgeneset_metadata_Tconvs[,c("studyGroup","type")], col=list(studyGroup=studygroupcolors, type=celltypecolors))
  
      
#PLOTTING HEATMAP(S) 
      #altering column names of matrix from full sample name --> sample ID (repalcing matrix colnames w/ metadata_reordered rownames)
      sridlist_treggeneset = as.list(metadata_treggeneset_reordered$srid)
      colnames(Treggeneset_expressionmatrix_reordered) <- sridlist_treggeneset
      
      sridlistTreg_treggeneset = as.list(treggeneset_metadata_Treg$srid)
      colnames(Treggenesetmatrix_tregs) <- sridlistTreg_treggeneset
      
      sridlistTconvs_treggeneset = as.list(treggeneset_metadata_Tconvs$srid)
      colnames(Treggenesetmatrix_tconvs) <- sridlistTconvs_treggeneset
      
      
      sridlist_thgeneset = as.list(metadata_thgeneset_reordered$srid)
      colnames(Thgeneset_expressionmatrix_reordered) <- sridlist_thgeneset
      
      sridlistTreg_thgeneset = as.list(thgeneset_metadata_Treg$srid)
      colnames(Thgenesetmatrix_tregs) <- sridlistTreg_thgeneset
      
      sridlistTconvs_thgeneset = as.list(thgeneset_metadata_Tconvs$srid)
      colnames(Thgenesetmatrix_tconvs) <- sridlistTconvs_thgeneset
      
      #Normalizing rows of matrix (making colors more apparent in heatmap)
          #1. Filter out genes that don't have enough variation (greater than or =3, basically we want to only keep rows (genes) that have at least 2 samples with unique non-zero expression values)   to be properly z-scored:
            Treggeneset_expressionmatrix_reordered <- Treggeneset_expressionmatrix_reordered[ apply(Treggeneset_expressionmatrix_reordered,1, function(r) length(unique(r))>=3), ]
            Treggenesetmatrix_tregs <- Treggenesetmatrix_tregs[ apply(Treggenesetmatrix_tregs,1, function(r) length(unique(r))>=3), ]
            Treggenesetmatrix_tconvs <- Treggenesetmatrix_tconvs[ apply(Treggenesetmatrix_tconvs,1, function(r) length(unique(r))>=3), ]
            
            Thgeneset_expressionmatrix_reordered <- Thgeneset_expressionmatrix_reordered[ apply(Thgeneset_expressionmatrix_reordered,1, function(r) length(unique(r))>=3), ]
            Thgenesetmatrix_tregs <- Thgenesetmatrix_tregs[ apply(Thgenesetmatrix_tregs,1, function(r) length(unique(r))>=3), ]
            Thgenesetmatrix_tconvs <- Thgenesetmatrix_tconvs[ apply(Thgenesetmatrix_tconvs,1, function(r) length(unique(r))>=3), ]
          #2.log normalize and z-score:
            Treggeneset_expressionmatrix_reordered <- t(scale(t(log(1+Treggeneset_expressionmatrix_reordered))))
            Treggenesetmatrix_tregs <- t(scale(t(log(1+Treggenesetmatrix_tregs))))
            Treggenesetmatrix_tconvs <- t(scale(t(log(1+Treggenesetmatrix_tconvs))))
            
            Thgeneset_expressionmatrix_reordered <- t(scale(t(log(1+Thgeneset_expressionmatrix_reordered))))
            Thgenesetmatrix_tregs <- t(scale(t(log(1+Thgenesetmatrix_tregs))))
            Thgenesetmatrix_tconvs <- t(scale(t(log(1+Thgenesetmatrix_tconvs))))
            
      #converting matrix [SIGEXPexpressionmatrix_reordered]from df--> back to matrix so it is compatible with complexheatmap()
      Treggeneset_expressionmatrix_reordered = as.matrix(Treggeneset_expressionmatrix_reordered)
      Treggenesetmatrix_tregs = as.matrix(Treggenesetmatrix_tregs)
      Treggenesetmatrix_tconvs = as.matrix(Treggenesetmatrix_tconvs)
      
      Thgeneset_expressionmatrix_reordered = as.matrix(Thgeneset_expressionmatrix_reordered)
      Thgenesetmatrix_tregs = as.matrix(Thgenesetmatrix_tregs)
      Thgenesetmatrix_tconvs = as.matrix(Thgenesetmatrix_tconvs)
      
      #Creating Plots 
      MasterHeatmap_treggeneset =
        Heatmap(
          Treggeneset_expressionmatrix_reordered,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Populations)",
          show_row_names = TRUE,
          cluster_columns = FALSE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Treg Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=masterAnnot_treggeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TregHeatmap_treggeneset =
        Heatmap(
          Treggenesetmatrix_tregs,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Treg Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = FALSE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Treg Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TregAnnot_treggeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TconvsHeatmap_treggeneset =
        Heatmap(
          Treggenesetmatrix_tconvs,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Tconvs Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = FALSE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Treg Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TconvsAnnot_treggeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      MasterHeatmap_thgeneset =
        Heatmap(
          Thgeneset_expressionmatrix_reordered,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Populations)",
          show_row_names = TRUE,
          cluster_columns = FALSE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Th Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=masterAnnot_thgeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TregHeatmap_thgeneset =
        Heatmap(
          Thgenesetmatrix_tregs,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Treg Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = FALSE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Th Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TregAnnot_thgeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TconvsHeatmap_thgeneset =
        Heatmap(
          Thgenesetmatrix_tconvs,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Tconvs Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = FALSE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Th Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TconvsAnnot_thgeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      
      ##Saving it into data folder in box
      pdf("GENESET_DifferentialExpressionHeatmaps(a).pdf", width=8, height=7)
      MasterHeatmap_treggeneset = draw(MasterHeatmap_treggeneset)
      TregHeatmap_treggeneset = draw(TregHeatmap_treggeneset)
      TconvsHeatmap_treggeneset = draw(TconvsHeatmap_treggeneset)
      MasterHeatmap_thgeneset = draw(MasterHeatmap_thgeneset)
      TregHeatmap_thgeneset = draw(TregHeatmap_thgeneset)
      TconvsHeatmap_thgeneset = draw(TconvsHeatmap_thgeneset)
      ggarrange(MasterHeatmap_treggeneset,TregHeatmap_treggeneset,TconvsHeatmap_treggeneset,MasterHeatmap_thgeneset,TregHeatmap_thgeneset,TconvsHeatmap_thgeneset) 
      dev.off()
      
      
      
      
  ##Creating a Clustered Version: 
      #Creating Plots 
      MasterHeatmap_treggeneset_clustered =
        Heatmap(
          Treggeneset_expressionmatrix_reordered,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Populations)",
          show_row_names = TRUE,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Treg Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=masterAnnot_treggeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TregHeatmap_treggeneset_clustered =
        Heatmap(
          Treggenesetmatrix_tregs,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Treg Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Treg Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TregAnnot_treggeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TconvsHeatmap_treggeneset_clustered =
        Heatmap(
          Treggenesetmatrix_tconvs,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Tconvs Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Treg Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TconvsAnnot_treggeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      MasterHeatmap_thgeneset_clustered =
        Heatmap(
          Thgeneset_expressionmatrix_reordered,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Populations)",
          show_row_names = TRUE,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Th Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=masterAnnot_thgeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TregHeatmap_thgeneset_clustered =
        Heatmap(
          Thgenesetmatrix_tregs,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Treg Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Th Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TregAnnot_thgeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TconvsHeatmap_thgeneset_clustered =
        Heatmap(
          Thgenesetmatrix_tconvs,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Tconvs Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Signature Th Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("pink2", "white","deepskyblue3")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TconvsAnnot_thgeneset,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      
      ##Saving it into data folder in box
      pdf("GENESET_DifferentialExpressionHeatmaps(b).pdf", width=8, height=7)
      MasterHeatmap_treggeneset_clustered = draw(MasterHeatmap_treggeneset_clustered)
      TregHeatmap_treggeneset_clustered = draw(TregHeatmap_treggeneset_clustered)
      TconvsHeatmap_treggeneset_clustered = draw(TconvsHeatmap_treggeneset_clustered)
      MasterHeatmap_thgeneset_clustered = draw(MasterHeatmap_thgeneset_clustered)
      TregHeatmap_thgeneset_clustered = draw(TregHeatmap_thgeneset_clustered)
      TconvsHeatmap_thgeneset_clustered = draw(TconvsHeatmap_thgeneset_clustered)
      ggarrange(MasterHeatmap_treggeneset_clustered,TregHeatmap_treggeneset_clustered,TconvsHeatmap_treggeneset_clustered,MasterHeatmap_thgeneset_clustered,TregHeatmap_thgeneset_clustered,TconvsHeatmap_thgeneset_clustered) 
      dev.off()

#RESAVING ALL FIGURES AS JPEGs 
      #output_path <- "C:/Users/klai/Box/2023SummerInternKenneth/data/official_figures/"
      #jpeg(file.path(output_path, "MasterDEA(expanded).jpg"), width = 1000, height = 1000)
      #print(MasterHeatmap)
      #dev.off()
      
      #jpeg(file.path(output_path, "TregDEA(expanded).jpg"), width = 500, height = 500)
      #print(TregHeatmap)
      #dev.off()
      
      #jpeg(file.path(output_path, "TconvsDEA(expanded).jpg"), width = 500, height = 500)
      #print(TconvsHeatmap)
      #dev.off()
      
      #jpeg(file.path(output_path, "MasterDEA_clustered(expanded).jpg"), width = 1000, height = 1000)
      #print(MasterHeatmap_clustered)
      #dev.off()
      
      #jpeg(file.path(output_path, "TregDEA_clustered(expanded).jpg"), width = 500, height = 500)
      #print(TregHeatmap_clustered)
      #dev.off()
      
      #jpeg(file.path(output_path, "TconvsDEA_clustered(expanded).jpg"), width = 500, height = 500)
      #print(TconvsHeatmap_clustered)
      #dev.off()