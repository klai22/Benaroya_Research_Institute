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

#IMPORTDATA 
  #Loading Distance Matrices 
    distmatA = fread("BRI_alpha_junctions_distmat.tsv")
      #distmatA <- read.table("BRI_alpha_junctions_distmat.tsv", header = TRUE, sep = "\t", as.is = TRUE)
    distmatB = fread("BRI_beta_junctions_distmat.tsv")
      #distmatB = read.table("/Users/klai/Box/2023SummerInternKenneth/data/BRI_beta_junctions_distmat.tsv",row.names=1,header=TRUE,sep="\t")
      #distmatB = read.table("/Users/klai/Box/2023SummerInternKenneth/data/BRI_beta_junctions_distmat.tsv",header=TRUE,sep="\t",na.strings = c(""),comment.char = "")
  #adding junction names to rownames of matrices 
      #rownames(distmatA) = colnames(distmatA)
      #rownames(distmatB) = colnames(distamatB)
      #dimnames(distmatA) <- list(rownames(distmatA), colnames(distmatA))
      #dimnames(distmatB) <- list(rownames(distmatB), colnames(distmatB))
  #converting them back to matrices instead of data tables 
    #distmatA <- as.matrix(distmatA)
    #distmatB <- as.matrix(distmatB)

  #Importing TCR metadata 
      BRITCRs = read.table("BRI_allTCRs.txt",row.names=1,header=TRUE,sep="\t")

#DATA ORGANIZING for DISTMATs 
  #SUBSET only ROWS/COLUMNS in distmat matrix that match junctions that come from (1)P390 project & (2)TRA or TRB, but based on row #s bc, matching by characters takes way too long
    #sub-setting only TCRs from P390 
      P390TCRs = BRITCRs %>% filter(grepl("^P390", project))
    #Sub-setting P390 TCRs for TRA + TRB ( we ignored TRG + TRDs)
      P390TCRA = subset(P390TCRs, P390TCRs$chain == "TRA")
      P390TCRB = subset(P390TCRs, P390TCRs$chain == "TRB")
    #Creating a vector of only each subsets' junctions 
      junctionsP390_a = P390TCRA$junction
      junctionsP390_b = P390TCRB$junction
    #Creating a vector/list of all the rownames in distmatA where the junction seq.s MATCHED b/w junctionsP390_a & distmatA (same for B)
      indicesA <- which( colnames(distmatA) %in% junctionsP390_a )
      indicesB <- which( colnames(distmatB) %in% junctionsP390_b )
    #Isolating only the rows where rows # (indices) of distmatA matches the rownumbers in the vector "indicesA" (same for B)
      P390distmatA <- distmatA[indicesA,..indicesA]
      P390distmatB <- distmatB[indicesB,..indicesB]
  #CONVERTING ROWNAMES from #s --> Junction Seq/Names (exact same as column-names) to enable downstream analyses/organizing 
      rownames(P390distmatA) <- colnames(P390distmatA)
      rownames(P390distmatB) <- colnames(P390distmatB)
 
  #SUBSETTING the distmats (EXPANDED) : only rows/columns that came from expanded cells ( we are only interested in this for now*)
    #Creating a vector of only the expanded junctions (appeared more than once) in each Project&chain 
      expanded_Ajunctions <- unique( P390TCRA$junction[duplicated(P390TCRA$junction)] )
      expanded_Bjunctions <- unique( P390TCRB$junction[duplicated(P390TCRB$junction)] )
    #Creating a vector/list of all the rownames in P390distmatA where the junction seq.s MATCHED b/w expanded_Ajunctions & P390distmatA (same for B)
      expindicesA <- which( colnames(P390distmatA) %in% expanded_Ajunctions )
      expindicesB <- which( colnames(P390distmatB) %in% expanded_Bjunctions )
    #Isolating only the rows where rows # (expindices) of P390distmatA matches the rownumbers in the vector "expindicesA" (same for B)
      EXPP390distmatA <- P390distmatA[expindicesA,..expindicesA]
      EXPP390distmatB <- P390distmatB[expindicesB,..expindicesB]
      #CONVERTING ROWNAMES from #s --> Junction Seq/Names (exact same as column-names) to enable downstream analyses/organizing 
      rownames(EXPP390distmatA) <- colnames(EXPP390distmatA)
      rownames(EXPP390distmatB) <- colnames(EXPP390distmatB)

      
#DATA ORGANIZING for Annotations 
      ###______________________###
#importdata (this was an issue w/ line 17659, there was a "#" which made an error when reading the df. "comment.char = """ helps accounf for this issue when doing read.table)
  ##Public Database including all TCR (CDR3) seq. and the epitopes they have been known to match [multiple epitopes/TCR]
  Database = read.table("/Users/klai/Box/2023SummerInternKenneth/data/2023-06-27/SearchTable-2023-05-16_16_23_38.755_vdjdb.tsv",header=TRUE,sep="\t",na.strings = c(""),comment.char = "")
  ##Data of all the TCRs sequenced in this project 
  toRemove <- read.table("/Users/klai/Box/2023SummerInternKenneth/data/libidsToRemove.txt")[,1]
  annoP390 <- read.table("/Users/klai/Box/2023SummerInternKenneth/data/P390_annotation.txt",header=TRUE,sep="\t")
  annoP390 <- annoP390[!annoP390$libid %in% toRemove,]
  ##Data from all the TCRs that interacted w/ an epitope [Karen's experiment: 1 epitope / TCR]
  stimByLibid=readRDS("/Users/klai/Box/2023SummerInternKenneth/data/stimByLibid.RDS")
  stimByLibid[ stimByLibid == "ZNT8/PPI"] <- "ZnT8+PP1"
  annoP390$stimPool <- stimByLibid[ annoP390$libid]
  projects <- unique(annoP390$project)
  ###Converting it into a df 
  epitopedata1= as.data.frame(stimByLibid, row.names=NULL)
  epitopedata2=rownames_to_column(epitopedata1, var = "libid")
  names(epitopedata2)[names(epitopedata2) == "data1"] <- "epitope"
  ###Checking nrow of each df made
  nrow(annoP390)
  nrow(BRITCRs)
  nrow(epitopedata2)
#First, I have to isolate + add on the junction sequences for each of the rows of annoP390 that match BRITCRs. I did this based on matching of project#. 
  annoP390BRITCRs <- merge(annoP390, BRITCRs, by = "libid")
  #annoP390BRITCRs <- merge(annoP390, BRITCRs["junction"], by = "libid")
  annoP390BRITCRs2 <- annoP390BRITCRs[, c("project.x","project.y", "libid", "studyGroup", "cellType", "junction","stimPool", "chain", "sampleRepositoryID")]
#Now, I will isolate only the sequences annoP390BRITCRs (TCR seq.s from this project) that reflect the sequences which match the public database(s). 
  DBCDR3match = annoP390BRITCRs2[annoP390BRITCRs2$junction %in% Database$CDR3, ]
  ##df1_filtered <- df1[df1$id %in% df2$id, ]
  ##figuring out how representative the DB was of our data = db represents 6.6% of the og data 
  nrow(annoP390BRITCRs2)
  nrow(DBCDR3match)
  #Now I will merge this data with the Database to include epitope information 
  ##rename junction column 
  names(Database)[3] <- "junction"
  ##merging database w/ experimental data (isolated)
  DBCDR3match2 <- merge(DBCDR3match, Database, by = "junction")
  ##decide which columns you want to keep 
  cols_to_remove <- c("project.y", "complex.id", "V", "J", "Reference", "Meta", "Method", "CDR3fix", "Score")  
  DBCDR3match3 <- DBCDR3match2[, -which(names(DBCDR3match2) %in% cols_to_remove)]
  #Adding a column that includes the unique combos of sampleIDs + CellTyps 
  DBCDR3match3$sampleIDcellType <- paste(DBCDR3match3$sampleRepositoryID, DBCDR3match3$cellType, sep = " ")
  ##converting all "Tscms" in cellType --> Tconvs 
  DBCDR3match3$cellType <- gsub("Tscm", "Tconventional", DBCDR3match3$cellType)
  
      ###________subset celltype according to annoP390BRITCRs2 (all BRITCRS) ______________###
  #SUBSETTING DISTMATs bc cell type 
    #Creating a vector of junctions unique to each cellType 
      #Subsetting dataframe of all same cell type 
        DBCDR3match3_Treg = subset(annoP390BRITCRs2, annoP390BRITCRs2$cellType == "Treg")
        DBCDR3match3_Tconvs = subset(annoP390BRITCRs2, annoP390BRITCRs2$cellType == "Tconventional")
      #Creating a list of junctions belong to said cell types 
        Junctionlist_Treg = as.list(DBCDR3match3_Treg$junction)
        Junctionlist_Tconvs = as.list(DBCDR3match3_Tconvs$junction)
      #Creating a vector/list of all the rownames in EXPP390distmatA where the junction seq.s MATCHED b/w Junctionlist_Treg (and then Tcvons) & EXPP390distmatA (same for B)
        expindicesA_Treg <- which( colnames(EXPP390distmatA) %in% Junctionlist_Treg )
        expindicesB_Treg <- which( colnames(EXPP390distmatB) %in% Junctionlist_Treg )
        expindicesA_Tconvs <- which( colnames(EXPP390distmatA) %in% Junctionlist_Tconvs )
        expindicesB_Tconvs <- which( colnames(EXPP390distmatB) %in% Junctionlist_Tconvs )
      #Isolating only the rows where rows # (expindices) of P390distmatA matches the rownumbers in the vector "expindicesA" (same for B)
        EXPP390distmatA_Treg <- EXPP390distmatA[expindicesA_Treg,..expindicesA_Treg]
        EXPP390distmatB_Treg <- EXPP390distmatB[expindicesB_Treg,..expindicesB_Treg]
        EXPP390distmatA_Tconvs <- EXPP390distmatA[expindicesA_Tconvs,..expindicesA_Tconvs]
        EXPP390distmatB_Tconvs <- EXPP390distmatB[expindicesB_Tconvs,..expindicesB_Tconvs]
      #CONVERTING ROWNAMES from #s --> Junction Seq/Names (exact same as column-names) to enable downstream analyses/organizing 
        rownames(EXPP390distmatA_Treg) <- colnames(EXPP390distmatA_Treg)
        rownames(EXPP390distmatB_Treg) <- colnames(EXPP390distmatB_Treg)
        rownames(EXPP390distmatA_Tconvs) <- colnames(EXPP390distmatA_Tconvs)
        rownames(EXPP390distmatB_Tconvs) <- colnames(EXPP390distmatB_Tconvs)
  
  
    ###________subset celltype according to DBCDR3match3 (isolated junctions that matched VDJDB ______________###
  #SUBSETTING DISTMATs bc cell type 
    #Creating a vector of junctions unique to each cellType 
      #Subsetting dataframe of all same cell type 
        #DBCDR3match3_Treg = subset(DBCDR3match3, DBCDR3match3$cellType == "Treg")
        #DBCDR3match3_Tconvs = subset(DBCDR3match3, DBCDR3match3$cellType == "Tconventional")
      #Creating a list of junctions belong to said cell types 
        #Junctionlist_Treg = as.list(DBCDR3match3_Treg$junction)
        #Junctionlist_Tconvs = as.list(DBCDR3match3_Tconvs$junction)
      #Creating a vector/list of all the rownames in EXPP390distmatA where the junction seq.s MATCHED b/w Junctionlist_Treg (and then Tcvons) & EXPP390distmatA (same for B)
        #expindicesA_Treg <- which( colnames(EXPP390distmatA) %in% Junctionlist_Treg )
        #expindicesB_Treg <- which( colnames(EXPP390distmatB) %in% Junctionlist_Treg )
        #expindicesA_Tconvs <- which( colnames(EXPP390distmatA) %in% Junctionlist_Tconvs )
        #expindicesB_Tconvs <- which( colnames(EXPP390distmatB) %in% Junctionlist_Tconvs )
      #Isolating only the rows where rows # (expindices) of P390distmatA matches the rownumbers in the vector "expindicesA" (same for B)
        #EXPP390distmatA_Treg <- EXPP390distmatA[expindicesA_Treg,..expindicesA_Treg]
        #EXPP390distmatB_Treg <- EXPP390distmatB[expindicesB_Treg,..expindicesB_Treg]
        #EXPP390distmatA_Tconvs <- EXPP390distmatA[expindicesA_Tconvs,..expindicesA_Tconvs]
        #EXPP390distmatB_Tconvs <- EXPP390distmatB[expindicesB_Tconvs,..expindicesB_Tconvs]
      #CONVERTING ROWNAMES from #s --> Junction Seq/Names (exact same as column-names) to enable downstream analyses/organizing 
        #rownames(EXPP390distmatA_Treg) <- colnames(EXPP390distmatA_Treg)
        #rownames(EXPP390distmatB_Treg) <- colnames(EXPP390distmatB_Treg)
        #rownames(EXPP390distmatA_Tconvs) <- colnames(EXPP390distmatA_Tconvs)
        #rownames(EXPP390distmatB_Tconvs) <- colnames(EXPP390distmatB_Tconvs)
  
  
       ###______________________###
  #CellType & StudyGroup Annotations 
      #Extracting list (exact order) of colnames (junctions) to create a new df
        TRAseqlist_Treg = as.list(colnames(EXPP390distmatA_Treg))
        TRBseqlist_Treg = as.list(colnames(EXPP390distmatB_Treg))
        TRAseqlist_Tconvs = as.list(colnames(EXPP390distmatA_Tconvs))
        TRBseqlist_Tconvs = as.list(colnames(EXPP390distmatB_Tconvs))
      #convert to df 
        TRAseqvsSGCT_Treg <- data.frame(junction = unlist(TRAseqlist_Treg))
        TRBseqvsSGCT_Treg <- data.frame(junction = unlist(TRBseqlist_Treg))
        TRAseqvsSGCT_Tconvs <- data.frame(junction = unlist(TRAseqlist_Tconvs))
        TRBseqvsSGCT_Tconvs <- data.frame(junction = unlist(TRBseqlist_Tconvs))
      #Add on metadata (all unmatched = "NA")....TRA had 2L of NA in cellType, TRB had 0 NAs for cellType
        TRAseqvsSGCT_Treg <- merge(TRAseqvsSGCT_Treg, annoP390BRITCRs2, by = "junction", all.x = TRUE)
        TRBseqvsSGCT_Treg <- merge(TRBseqvsSGCT_Treg, annoP390BRITCRs2, by = "junction", all.x = TRUE)
        TRAseqvsSGCT_Tconvs <- merge(TRAseqvsSGCT_Tconvs, annoP390BRITCRs2, by = "junction", all.x = TRUE)
        TRBseqvsSGCT_Tconvs <- merge(TRBseqvsSGCT_Tconvs, annoP390BRITCRs2, by = "junction", all.x = TRUE)
        #getting rid of TSCM 
        TRAseqvsSGCT_Treg$cellType <- gsub("Tscm", "Tconventional", TRAseqvsSGCT_Treg$cellType)
        TRBseqvsSGCT_Treg$cellType <- gsub("Tscm", "Tconventional", TRBseqvsSGCT_Treg$cellType)
        TRAseqvsSGCT_Tconvs$cellType <- gsub("Tscm", "Tconventional", TRAseqvsSGCT_Tconvs$cellType)
        TRBseqvsSGCT_Tconvs$cellType <- gsub("Tscm", "Tconventional", TRBseqvsSGCT_Tconvs$cellType)
      #ISSUE: After merging dfs based on junction, the dimensions increased (ruining the order) bc unique junctions appear multiple times if it was found in 2 cell Types or different metadata values 
        # Concatenate values for duplicated junctions (if junction appeared multiple times in metadata, I concatenated their values together to keep og dimensions)
          TRAseqvsSGCT_Treg <- aggregate(. ~ junction, TRAseqvsSGCT_Treg, FUN = function(x) paste(unique(x), collapse = ","), na.action = na.pass)
          TRBseqvsSGCT_Treg <- aggregate(. ~ junction, TRBseqvsSGCT_Treg, FUN = function(x) paste(unique(x), collapse = ","), na.action = na.pass)
          TRAseqvsSGCT_Tconvs <- aggregate(. ~ junction, TRAseqvsSGCT_Tconvs, FUN = function(x) paste(unique(x), collapse = ","), na.action = na.pass)
          TRBseqvsSGCT_Tconvs <- aggregate(. ~ junction, TRBseqvsSGCT_Tconvs, FUN = function(x) paste(unique(x), collapse = ","), na.action = na.pass)
        # Now reordering the dfs to match that of the distmat rows / columns (their junctions)
          TRAseqvsSGCT_Treg <- TRAseqvsSGCT_Treg[match(TRAseqlist_Treg, TRAseqvsSGCT_Treg$junction), ]
          TRBseqvsSGCT_Treg <- TRBseqvsSGCT_Treg[match(TRBseqlist_Treg, TRBseqvsSGCT_Treg$junction), ]
          TRAseqvsSGCT_Tconvs <- TRAseqvsSGCT_Tconvs[match(TRAseqlist_Tconvs, TRAseqvsSGCT_Tconvs$junction), ]
          TRBseqvsSGCT_Tconvs <- TRBseqvsSGCT_Tconvs[match(TRBseqlist_Tconvs, TRBseqvsSGCT_Tconvs$junction), ]
          #df_reordered <- df[match(desired_order, df$Fruit), ]
      #Creating Ind. Subtables 
          #TRAseqvscelltype = TRAseqvsSGCT[, c( "cellType")]
          #TRAseqvsstudyGroup = TRAseqvsSGCT[, c("studyGroup")]
          #TRBseqvscelltype = TRBseqvsSGCT[, c("cellType")]
          #TRBseqvsstudyGroup = TRBseqvsSGCT[, c("studyGroup")]
      #renaming the columns titles 
          #TRAseqvscelltype = as.data.frame(TRAseqvscelltype)
          #colnames(TRAseqvscelltype) <- c("cellType")
          #TRAseqvsstudyGroup = as.data.frame(TRAseqvsstudyGroup)
          #colnames(TRAseqvsstudyGroup) <- c("studyGroup")
          #TRBseqvscelltype = as.data.frame(TRBseqvscelltype)
          #colnames(TRBseqvscelltype) <- c("cellType")
          #TRBseqvsstudyGroup = as.data.frame(TRBseqvsstudyGroup)
          #colnames(TRBseqvsstudyGroup) <- c("studyGroup")
#PLOTHEATMAP 
      #convert table--> matrix so it is compatible w/ complexHeatmap 
      EXPP390distmatA1_Treg = as.matrix(EXPP390distmatA_Treg)
      rownames(EXPP390distmatA1_Treg) <- colnames(EXPP390distmatA_Treg)
      colnames(EXPP390distmatA1_Treg) <- colnames(EXPP390distmatA_Treg)
      
      EXPP390distmatB1_Treg = as.matrix(EXPP390distmatB_Treg)
      rownames(EXPP390distmatB1_Treg) <- colnames(EXPP390distmatB_Treg)
      colnames(EXPP390distmatB1_Treg) <- colnames(EXPP390distmatB_Treg)
      
      EXPP390distmatA1_Tconvs = as.matrix(EXPP390distmatA_Tconvs)
      rownames(EXPP390distmatA1_Tconvs) <- colnames(EXPP390distmatA_Tconvs)
      colnames(EXPP390distmatA1_Tconvs) <- colnames(EXPP390distmatA_Tconvs)
      
      EXPP390distmatB1_Tconvs = as.matrix(EXPP390distmatB_Tconvs)
      rownames(EXPP390distmatB1_Tconvs) <- colnames(EXPP390distmatB_Tconvs)
      colnames(EXPP390distmatB1_Tconvs) <- colnames(EXPP390distmatB_Tconvs)
      
      
      #Creating the Annotations 
      #TRACTAnnot = HeatmapAnnotation(bar = TRAseqvscelltype$cellType, col = list(cellType = c("Treg" = "blue", "Tconventional" = "yellow","Tconventional,Treg"="green", "Treg,Tconventional"="green","NA" = "white")))
      #TRBCTAnnot = HeatmapAnnotation(bar = TRBseqvscelltype$cellType, col = list(cellType = c("Treg" = "blue", "Tconventional" = "yellow","Tconventional,Treg"="green", "Treg,Tconventional"="green","NA" = "white")))
      #TRASGAnnot = HeatmapAnnotation(bar = TRAseqvsstudyGroup$studyGroup, col = list(cellType = c("Control" = "aquamarine2", "T1D" = "palevioletred2","T1D,Control"="purple3", "Control,T1D"="purple3","NA" = "white")))
      #TRBSGAnnot = HeatmapAnnotation(bar = TRBseqvsstudyGroup$studyGroup, col = list(cellType = c("Control" = "aquamarine2", "T1D" = "palevioletred2","T1D,Control"="purple3", "Control,T1D"="purple3","NA" = "white")))
        #Setting color schemes 
          celltypecolors <- c("Treg" = "blue", "Tconventional" = "yellow",
                          "Tconventional,Treg"="green", "Treg,Tconventional"="green","NA" = "white") # defining the colors for the annotation
          studygroupcolors <- c("Control" = "aquamarine2", "T1D" = "palevioletred2",
                            "T1D,Control"="purple3", "Control,T1D"="purple3","NA" = "white") # defining the colors for the annotation
        #creating annotation obj.s 
          
          TRAAnnot_Treg = HeatmapAnnotation( df=TRAseqvsSGCT_Treg[,c("studyGroup","cellType")], col=list(studyGroup=studygroupcolors, cellType=celltypecolors))
      
          TRBAnnot_Treg = HeatmapAnnotation( df=TRBseqvsSGCT_Treg[,c("studyGroup","cellType")], col=list(studyGroup=studygroupcolors, cellType=celltypecolors))
          
          TRAAnnot_Tconvs = HeatmapAnnotation( df=TRAseqvsSGCT_Tconvs[,c("studyGroup","cellType")], col=list(studyGroup=studygroupcolors, cellType=celltypecolors))
          
          TRBAnnot_Tconvs = HeatmapAnnotation( df=TRBseqvsSGCT_Tconvs[,c("studyGroup","cellType")], col=list(studyGroup=studygroupcolors, cellType=celltypecolors))
      
      #Plot Heat Maps  
          #TEMPLATE: Heatmap(EXPP390distmatA, show_column_names=TRUE,top_annotation=ha,cluster_columns=FALSE,name="DISTMATA",column_title="TCR Sequence")
          
      ###Attempt #2 ...
          Expanded_TRA_Distmat_Heatmap_Treg =
            Heatmap(
              EXPP390distmatA1_Treg,
              show_column_names = FALSE,
              column_title = "TREG Alpha Junction TCR Distmat",
              show_row_names = TRUE,
              cluster_columns = TRUE,
              name = "TCR Distance",
              row_title = "TRA Sequences",
              row_names_gp = gpar(fontsize = 4, fontface = "bold"),
              col = colorRampPalette(c("deepskyblue3", "white","pink2"))(100),
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10, fontface = "bold"),
                labels_gp = gpar(fontsize = 6)
              ),
              row_names_side = "right",
              top_annotation=TRAAnnot_Treg,
              show_row_dend = FALSE,
              show_column_dend = FALSE
            )
            
          Expanded_TRB_Distmat_Heatmap_Treg =
            Heatmap(
              EXPP390distmatB1_Treg,
              show_column_names = FALSE,
              column_title = "TREG Beta Junction TCR Distmat",
              show_row_names = TRUE,
              cluster_columns = TRUE,
              name = "TCR Distance",
              row_title = "TRB Sequences",
              row_names_gp = gpar(fontsize = 4, fontface = "bold"),
              col = colorRampPalette(c("deepskyblue3", "white","pink2"))(100),
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10, fontface = "bold"),
                labels_gp = gpar(fontsize = 6)
              ),
              row_names_side = "right",
              top_annotation=TRBAnnot_Treg,
              show_row_dend = FALSE,
              show_column_dend = FALSE
            )
          
          Expanded_TRA_Distmat_Heatmap_Tconvs =
            Heatmap(
              EXPP390distmatA1_Tconvs,
              show_column_names = FALSE,
              column_title = "TCONVS Alpha Junction TCR Distmat",
              show_row_names = TRUE,
              cluster_columns = TRUE,
              name = "TCR Distance",
              row_title = "TRA Sequences",
              row_names_gp = gpar(fontsize = 4, fontface = "bold"),
              col = colorRampPalette(c("deepskyblue3", "white","pink2"))(100),
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10, fontface = "bold"),
                labels_gp = gpar(fontsize = 6)
              ),
              row_names_side = "right",
              top_annotation=TRAAnnot_Tconvs,
              show_row_dend = FALSE,
              show_column_dend = FALSE
            )
          
          Expanded_TRB_Distmat_Heatmap_Tconvs =
            Heatmap(
              EXPP390distmatB1_Tconvs,
              show_column_names = FALSE,
              column_title = "TCONVS Beta Junction TCR Distmat",
              show_row_names = TRUE,
              cluster_columns = TRUE,
              name = "TCR Distance",
              row_title = "TRB Sequences",
              row_names_gp = gpar(fontsize = 4, fontface = "bold"),
              col = colorRampPalette(c("deepskyblue3", "white","pink2"))(100),
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10, fontface = "bold"),
                labels_gp = gpar(fontsize = 6)
              ),
              row_names_side = "right",
              top_annotation=TRBAnnot_Tconvs,
              show_row_dend = FALSE,
              show_column_dend = FALSE
            )
            

      ##Saving it into data folder in box
      pdf("ExpandedTCRDistmatHeatmaps(4).pdf", width=13, height=13)
      Expanded_TRA_Distmat_Heatmap_Treg = draw(Expanded_TRA_Distmat_Heatmap_Treg)
      Expanded_TRB_Distmat_Heatmap_Treg = draw(Expanded_TRB_Distmat_Heatmap_Treg)
      Expanded_TRA_Distmat_Heatmap_Tconvs = draw(Expanded_TRA_Distmat_Heatmap_Tconvs)
      Expanded_TRB_Distmat_Heatmap_Tconvs = draw(Expanded_TRB_Distmat_Heatmap_Tconvs)
      ggarrange(Expanded_TRA_Distmat_Heatmap_Treg, Expanded_TRB_Distmat_Heatmap_Treg,Expanded_TRA_Distmat_Heatmap_Tconvs, Expanded_TRB_Distmat_Heatmap_Tconvs) 
      dev.off()
      
#RESAVING AS JPEGS 
      output_path <- "C:/Users/klai/Box/2023SummerInternKenneth/data/official_figures/"
      jpeg(file.path(output_path, "TRADISTMAT(expanded)TREG.jpg"), width = 1000, height = 1000)
      print(Expanded_TRA_Distmat_Heatmap_Treg)
      dev.off()
      
      jpeg(file.path(output_path, "TRBDISTMAT(expanded)TREG.jpg"), width = 1000, height = 1000)
      print(Expanded_TRB_Distmat_Heatmap_Treg)
      dev.off()
      
      output_path <- "C:/Users/klai/Box/2023SummerInternKenneth/data/official_figures/"
      jpeg(file.path(output_path, "TRADISTMAT(expanded)TCONVS.jpg"), width = 1000, height = 1000)
      print(Expanded_TRA_Distmat_Heatmap_Tconvs)
      dev.off()
      
      jpeg(file.path(output_path, "TRBDISTMAT(expanded)TCONVS.jpg"), width = 1000, height = 1000)
      print(Expanded_TRB_Distmat_Heatmap_Tconvs)
      dev.off()
      
      
      
    


###______SCRAP______###

      ## Plot a heatmap with the genes DE in the unstiulated condition between the knockouts
      ens <- rownames(nostimKODE)[ nostimKODE$adj.P.Val <= 0.05] # nostimKODE is a data frame containing differential expression adjusted p-values. Each row is a gene.
      design_qc <- design_qc[ order(design_qc$stimulation, design_qc$studyGroup),] # design_qc is the metadata
      libs <- design_qc$libid # libs are the samples
      tp <- t(scale(t(log2(0.5+counts_pc_norm[ens, libs])))) # counts_pc_norm is the gene expression matrix. Here, I'm log-transforming the gene expression data and then z-scoring it
      rownames(tp) <- nostimKODE[ens,"mgi_symbol"] # Here I'm renaming the rows of the matrix to be the gene symbol rather than the ensembl ID. mgi_symbol is more mouse, but your data will have "HGNC.symbol" probably, for human.
      genotypeColors <- c("D1A dlckCre+"="red","Wildtype"="blue") # defining the colors for the annotation
      stimulationColors <- c("Il-6, Il-1B, anti-Il-4, anti-Il-12, anti-IFNg, TGFb"="orange","none"="grey") # defining the colors for the annotation
      
      ha <- HeatmapAnnotation( df=design_qc[,c("studyGroup","stimulation")], col=list(studyGroup=genotypeColors, stimulation=stimulationColors))
      png("../../data/2023-07-05/plots/14DEGenes_Dyrk1aKO_heatmap.png",width=600,height=300)
      Heatmap(tp, show_column_names=FALSE,top_annotation=ha, cluster_columns=FALSE,name="scaled\nlog\nexpression",column_title="14 genes DE at <=5% FDR\nbetween Dyrk1a-KO and Wiltype\nin unstimulated cells")
      dev.off() 