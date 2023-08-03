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
      
#FILTER DATA 
  # I want to filter out only the sig. genes (expression values in matrix) with FDR-differential-expression-values [adj. p-value] that were =>5% (0.05) [FDR values in expressionresults]
      #Make a list of the genes that meet his req. 
      SIGgenesTreg = subset(EXPTregexpressionresults,adj.P.Val <= 0.05 )
      SIGgenesTreglist = as.list(rownames(SIGgenesTreg))
      
      SIGgenesTconvs = subset(EXPTconvsexpressionresults,adj.P.Val <= 0.05 )
      SIGgenesTconvslist = as.list(rownames(SIGgenesTconvs))
      
      #Make a master list of these genes (eliminating duplicates) [there were no matches between these 2 lists] 
      allSIGgeneslist = union(SIGgenesTreglist, SIGgenesTconvslist)
      
      #subset the matrix so that it only includes the genes that meet the FDR-differential-expression-threshold (0.05) (isolate rows w/ rownames match the list of sig genes only)
      SIGEXPexpressionmatrix = EXPexpressionmatrixdata[rownames(EXPexpressionmatrixdata) %in% allSIGgeneslist, ]
      
      
#ORGANIZING DATA / PREPPING FOR HEATMAP ANNOTATIONS (ordering matrix and pre-annotations data properly)
  #Re-ordering / grouping matrix columns based on control group (beginning of column name string) and celltype (ending of column name string)
      
      #converting (SIGEXPexpressionmatrix) from matrix --> df because otherwise the code below won't work 
      SIGEXPexpressionmatrix = as.data.frame(SIGEXPexpressionmatrix)
      
      # Columns that start with "T1D" and end with "Treg"
      t1dtreg_cols <- grep("^T1D.*Treg$", names(SIGEXPexpressionmatrix))
      
      # Columns that start with "T1D" and end with "Tconventional"
      t1dtconvs_cols <- grep("^T1D.*Tconventional$", names(SIGEXPexpressionmatrix))
      
      # Columns that start with "Control" and end with "Treg"
      ctrltreg_cols <- grep("^Control.*Treg$", names(SIGEXPexpressionmatrix))
      
      # Columns that start with "Control" and end with "Tconventional"
      ctrltconvs_cols <- grep("^Control.*Tconventional$", names(SIGEXPexpressionmatrix))
      
      # Order the columns based on the grouping
      samplegroups <- c(t1dtreg_cols,t1dtconvs_cols,ctrltreg_cols,ctrltconvs_cols)
      
      # Subset the dataframe with the ordered columns
      SIGEXPexpressionmatrix_reordered <- SIGEXPexpressionmatrix[, samplegroups]
      
      #NOTE: right now, SIGEXPexpressionmatrix_reordered is a df, if encounter problems need to turn df-->matrix 
      
  #Re-ordering metadata to match the matrix (the primary data represented in heatmap )
      
      #Isolate only the rows of metadata that are represented in the matrix (only the samples that had significant enough expression in given gene(s))
      SIGsampleslist = as.list(colnames(SIGEXPexpressionmatrix_reordered))
      
      #subset the metadata so that it only includes the samples that match matrix (it seemed that EXP metadata already only had the significant sample names?, not sure if this was a coincidence or it already had been subsetted for signifigant metadata to begin with)
      SIGmetadata = EXPmetadata[rownames(EXPmetadata) %in% SIGsampleslist, ]
      
      #Reordering the metadata rows so that it matches the same order as the matrix 
      SIGmetadata_reordered <- SIGmetadata[match(SIGsampleslist, rownames(SIGmetadata)), ]
      
#SUBSETTING MATRICES + METADATA FOR CELLTYPE IND. PLOTS 
      #Creating separate matrices for cell-type plots 
        #list of Treg & Tconvs samples 
        SIGTregsamples = as.list(rownames(SIGmetadata_reordered[SIGmetadata_reordered$type == "Treg", ]))
        SIGTconvssamples = as.list(rownames(SIGmetadata_reordered[SIGmetadata_reordered$type == "Tconventional", ]))
      
        #create subset matrices for each cell type 
       SIGEXPexpressionmatrix_reordered_Treg <- SIGEXPexpressionmatrix_reordered[, names(SIGEXPexpressionmatrix_reordered) %in% SIGTregsamples]
       SIGEXPexpressionmatrix_reordered_Tconvs <- SIGEXPexpressionmatrix_reordered[, names(SIGEXPexpressionmatrix_reordered) %in% SIGTconvssamples]
       
       #Create subset metadata for each celltype 
       SIGmetadata_reordered_Treg = subset(SIGmetadata_reordered, type == "Treg")
       SIGmetadata_reordered_Tconvs = subset(SIGmetadata_reordered, type == "Tconventional")
       
       
#RENAMING of matrices so that gene ensemble IDs match the HGNC.symbol (actual gene name) from metadata 
       #making a df of all the SIG ensemble IDs & their matching HGNC symbols 
       SIGgenesTreg_symbols = as.data.frame(SIGgenesTreg$HGNC.symbol, rownames(SIGgenesTreg))
       SIGgenesTconvs_symbols = as.data.frame(SIGgenesTconvs$HGNC.symbol, rownames(SIGgenesTconvs))
       
       SIGgenesTreg_symbols$symbol <- SIGgenesTreg_symbols$`SIGgenesTreg$HGNC.symbol`
       SIGgenesTreg_symbols$`SIGgenesTreg$HGNC.symbol` = NULL
       SIGgenesTconvs_symbols$symbol <- SIGgenesTconvs_symbols$`SIGgenesTconvs$HGNC.symbol`
       SIGgenesTconvs_symbols$`SIGgenesTconvs$HGNC.symbol` = NULL

       #making a df for both of them combined (for the master plot) (found no duplicates)
        SIGgenesall_symbols <- rbind(SIGgenesTreg_symbols, SIGgenesTconvs_symbols) 
       
       
       #Renaming matrix + results 
        matching_indices_matrix <- match(rownames(SIGEXPexpressionmatrix_reordered), rownames(SIGgenesall_symbols))
        rownames(SIGEXPexpressionmatrix_reordered) <- SIGgenesall_symbols$symbol[matching_indices_matrix]
        
        matching_indices_matrix_Treg <- match(rownames(SIGEXPexpressionmatrix_reordered_Treg), rownames(SIGgenesall_symbols))
        rownames(SIGEXPexpressionmatrix_reordered_Treg) <- SIGgenesall_symbols$symbol[matching_indices_matrix_Treg]
        
        matching_indices_matrix_Tconvs <- match(rownames(SIGEXPexpressionmatrix_reordered_Tconvs), rownames(SIGgenesall_symbols))
        rownames(SIGEXPexpressionmatrix_reordered_Tconvs) <- SIGgenesall_symbols$symbol[matching_indices_matrix_Tconvs]
        

#CREATING HEATMAP ANNOTATIONS 
      
       #Setting color schemes 
      celltypecolors <- c("Treg" = "blue", "Tconventional" = "yellow") # defining the colors for the annotation
      studygroupcolors <- c("Control" = "aquamarine2", "T1D" = "palevioletred2") # defining the colors for the annotation
      #creating annotation obj.s 
      
      masterAnnot = HeatmapAnnotation( df=SIGmetadata_reordered[,c("studyGroup","type")], col=list(studyGroup=studygroupcolors, type=celltypecolors))
      TregAnnot = HeatmapAnnotation( df=SIGmetadata_reordered_Treg[,c("studyGroup","type")], col=list(studyGroup=studygroupcolors, type=celltypecolors))
      TconvsAnnot = HeatmapAnnotation( df=SIGmetadata_reordered_Tconvs[,c("studyGroup","type")], col=list(studyGroup=studygroupcolors, type=celltypecolors))
  
      
#PLOTTING HEATMAP(S) 
      #altering column names of matrix from full sample name --> sample ID (repalcing matrix colnames w/ SIGmetadata_reordered rownames)
      sridlist = as.list(SIGmetadata_reordered$srid)
      colnames(SIGEXPexpressionmatrix_reordered) <- sridlist
      
      sridlistTreg = as.list(SIGmetadata_reordered_Treg$srid)
      colnames(SIGEXPexpressionmatrix_reordered_Treg) <- sridlistTreg
      
      sridlistTconvs = as.list(SIGmetadata_reordered_Tconvs$srid)
      colnames(SIGEXPexpressionmatrix_reordered_Tconvs) <- sridlistTconvs
      
      
      #Normalizing rows of matrix (making colors more apparent in heatmap)
          #1. Filter out genes that don't have enough variation (greater than or =3, basically we want to only keep rows (genes) that have at least 2 samples with unique non-zero expression values)   to be properly z-scored:
            SIGEXPexpressionmatrix_reordered <- SIGEXPexpressionmatrix_reordered[ apply(SIGEXPexpressionmatrix_reordered,1, function(r) length(unique(r))>=3), ]
            SIGEXPexpressionmatrix_reordered_Treg <- SIGEXPexpressionmatrix_reordered_Treg[ apply(SIGEXPexpressionmatrix_reordered_Treg,1, function(r) length(unique(r))>=3), ]
            SIGEXPexpressionmatrix_reordered_Tconvs <- SIGEXPexpressionmatrix_reordered_Tconvs[ apply(SIGEXPexpressionmatrix_reordered_Tconvs,1, function(r) length(unique(r))>=3), ]
          #2.log normalize and z-score:
            SIGEXPexpressionmatrix_reordered <- t(scale(t(log(1+SIGEXPexpressionmatrix_reordered))))
            SIGEXPexpressionmatrix_reordered_Treg <- t(scale(t(log(1+SIGEXPexpressionmatrix_reordered_Treg))))
            SIGEXPexpressionmatrix_reordered_Tconvs <- t(scale(t(log(1+SIGEXPexpressionmatrix_reordered_Tconvs))))
            
            
      #converting matrix [SIGEXPexpressionmatrix_reordered]from df--> back to matrix so it is compatible with complexheatmap()
      SIGEXPexpressionmatrix_reordered = as.matrix(SIGEXPexpressionmatrix_reordered)
      SIGEXPexpressionmatrix_reordered_Treg = as.matrix(SIGEXPexpressionmatrix_reordered_Treg)
      SIGEXPexpressionmatrix_reordered_Tconvs = as.matrix(SIGEXPexpressionmatrix_reordered_Tconvs)
      
      #Creating Plots 
      MasterHeatmap =
        Heatmap(
          SIGEXPexpressionmatrix_reordered,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Populations)",
          show_row_names = TRUE,
          cluster_columns = FALSE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("deepskyblue3", "white","pink2")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=masterAnnot,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TregHeatmap =
        Heatmap(
          SIGEXPexpressionmatrix_reordered_Treg,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Treg Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = FALSE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("deepskyblue3", "white","pink2")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TregAnnot,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TconvsHeatmap =
        Heatmap(
          SIGEXPexpressionmatrix_reordered_Tconvs,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Tconvs Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = FALSE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("deepskyblue3", "white","pink2")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TconvsAnnot,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      
      ##Saving it into data folder in box
      pdf("DifferentialExpressionHeatmaps(5a).pdf", width=8, height=7)
      MasterHeatmap = draw(MasterHeatmap)
      TregHeatmap = draw(TregHeatmap)
      TconvsHeatmap = draw(TconvsHeatmap)
      ggarrange(MasterHeatmap,TregHeatmap,TconvsHeatmap) 
      dev.off()
      
      
  ##Creating a Clustered Version: 
      #Creating Plots 
      MasterHeatmap_clustered =
        Heatmap(
          SIGEXPexpressionmatrix_reordered,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Populations)",
          show_row_names = TRUE,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("deepskyblue3", "white","pink2")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=masterAnnot,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TregHeatmap_clustered =
        Heatmap(
          SIGEXPexpressionmatrix_reordered_Treg,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Treg Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("deepskyblue3", "white","pink2")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TregAnnot,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      TconvsHeatmap_clustered =
        Heatmap(
          SIGEXPexpressionmatrix_reordered_Tconvs,
          show_column_names = TRUE,
          column_title = "Sample (Expanded Tconvs Clonotypes)",
          show_row_names = TRUE,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          name = "Expression Level",
          row_title = "Gene (transcript) Name",
          row_names_gp = gpar(fontsize = 5, fontface = "bold"),
          column_names_gp = gpar(fontsize = 6, fontface = "bold"),
          col = colorRamp2(c(-2, 0, 2),c("deepskyblue3", "white","pink2")),
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)
          ),
          row_names_side = "left",
          top_annotation=TconvsAnnot,
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          height = 8,
          width = 8
        )
      
      ##Saving it into data folder in box
      pdf("DifferentialExpressionHeatmaps(5b).pdf", width=8, height=7)
      MasterHeatmap_clustered = draw(MasterHeatmap_clustered)
      TregHeatmap_clustered = draw(TregHeatmap_clustered)
      TconvsHeatmap_clustered = draw(TconvsHeatmap_clustered)
      ggarrange(MasterHeatmap_clustered,TregHeatmap_clustered,TconvsHeatmap_clustered) 
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
      
      
    #general plan : filter by 5% FDR genes only, then by Tregs/Tconvs (by sampleID), annotate metadata prep, make heatmap 
      
      
      