#setwd
setwd("C:/Users/klai/Box/2023SummerInternKenneth/data/2023-06-27")

#load packages 
library("ggpubr")
library("ggplot2")
library("dplyr")
library("tidyr")
library("lubridate")
library("tibble")
library("igraph")
library("ComplexHeatmap")

#importdata (this was an issue w/ line 17659, there was a "#" which made an error when reading the df. "comment.char = """ helps accounf for this issue when doing read.table)
  ##Public Database including all TCR (CDR3) seq. and the epitopes they have been known to match [multiple epitopes/TCR]
  Database = read.table("SearchTable-2023-05-16_16_23_38.755_vdjdb.tsv",header=TRUE,sep="\t",na.strings = c(""),comment.char = "")
  ##all the TCR (CDR3) seq.s from BRI ##NOTE: "junction" = "TCR seq." = "CDR3 Seq."
  BRITCRs = read.table("/Users/klai/Box/2023SummerInternKenneth/data/BRI_allTCRs.txt",row.names=1,header=TRUE,sep="\t")
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
  DBCDR3match3$cellType <- gsub("Tscm", "Tconventional", DBCDR3match3$cellType)
  DBCDR3match3$sampleIDcellType <- paste(DBCDR3match3$sampleRepositoryID, DBCDR3match3$cellType, sep = " ")
  ##converting all "Tscms" in cellType --> Tconvs 
  
#Calculating Proportion: total # of libids (from a given sub-project) [which matched seq. in VDJDB] that match a given epitope / total # of libids in given project (including those that did not have matches in VDJDB)
####MAKING THIS A FUNCTION SO THAT I CAN GENERATE PLOTS FOR EACH EPITOPE LVL: epitope, epitope gene, epitope spp.########
##note: "Epitope" is just a default, you can just delete the "= "Epitope"" and it would still work, function knows to take whatever word you give it in "" when you run it. 
generateCDR3vsVDJDBplot <- function( epitopeColumnName = "Epitope"){
  #Numerator 
  ## Initialize empty vectors to store the results
  numerator <- vector()
  epitope <- vector()
  sampleIDcellType <- vector()
  studyGroup <- vector()
  cellType <- vector()
  sampleRepositoryID <- vector()
  ## Loop through unique srID-epitope combinations
  for (s in unique(DBCDR3match3$sampleIDcellType)) {
    for (e in unique(DBCDR3match3[DBCDR3match3$sampleIDcellType == s, epitopeColumnName])) {
      # Subset the dataframe for the current srID-epitope combination
      subset_df <- DBCDR3match3[DBCDR3match3$sampleIDcellType == s & DBCDR3match3[, epitopeColumnName] == e, ]
      # Calculate the numerator (count of unique libids)
      num <- length(unique(subset_df$libid))
      # Store the results in the vectors
      numerator <- c(numerator, num)
      epitope <- c(epitope, e)
      sampleIDcellType <- c(sampleIDcellType, s)
      studyGroup <- c(studyGroup, unique(subset_df$studyGroup))
      cellType <- c(cellType, unique(subset_df$cellType))
      sampleRepositoryID <- c(sampleRepositoryID, unique(subset_df$sampleRepositoryID))
    }
  }
  ## Create the new dataframe
  new_df <- data.frame(numerator, epitope, studyGroup, cellType, sampleIDcellType, sampleRepositoryID)
  ## Print the new dataframe
  print(new_df)
  
  
  #Denominator = total libids / project (non-VDJDB isolates)
  ##Rename project sampleRepositoryID in annoP390BRITCR2
  #annoP390BRITCRs3 <- annoP390BRITCRs2 %>%
    #rename(project = project.x)
  #Adding a column that includes the unique combos of sampleIDs + CellTyps 
  annoP390BRITCRs2$sampleIDcellType <- paste(annoP390BRITCRs2$sampleRepositoryID, annoP390BRITCRs2$cellType, sep = " ")
  ## Calculate the total number of libids for each project in df2 (annoP390BRITCRs3)
  libid_counts <- annoP390BRITCRs2 %>%
    group_by(sampleIDcellType) %>%
    summarise(denominator = n())
  ## Merge libid_counts with new_df based on the project column
  new_df <- new_df %>%
    left_join(libid_counts, by = "sampleIDcellType")
  ## Print the updated new_df dataframe
  print(new_df)
  
  
  #Calc proportions 
  new_df$proportion = new_df$numerator / new_df$denominator
  
  
  #PLOT PROPTIONS
  ##converting all "Tscms" in cellType --> Tconvs 
  new_df$cellType <- gsub("Tscm", "Tconventional", new_df$cellType)
  ##eliminating any epitope plots that do not represent both cellTypes ( we only want plots w/ proportion data for BOTH Tconvs + Tregs) 
  ### Step 1: Count the number of unique cellTypes for each epitope
  celltype_count <- aggregate(cellType ~ epitope, new_df, function(x) length(unique(x)))
  ### Step 2: Filter the dataframe to keep only individuals who ate both apples and pears
  filtered_epitopes <- celltype_count$epitope[celltype_count$cellType == 2]
  df_filtered <- new_df[new_df$epitope %in% filtered_epitopes, ]
  ##plotting 
  VDJDBepitopeplot=ggplot(df_filtered[!is.na(df_filtered$epitope),], aes(x=cellType,y=proportion, color=studyGroup, group=sampleRepositoryID)) + geom_boxplot(aes(x=cellType, group=cellType),color="black") + geom_line() + geom_point() + facet_wrap(~epitope)
  ##ending the function and indicating the plot object as the official output of the function you made 
  
  return(VDJDBepitopeplot)
}


#Applying the function you made to actually generate the 3 plots of interest 
  CDR3vsVDJDB_epitope <- generateCDR3vsVDJDBplot("Epitope")
  CDR3vsVDJDB_epitopegene <- generateCDR3vsVDJDBplot("Epitope.gene")
  CDR3vsVDJDB_epitopespp <- generateCDR3vsVDJDBplot("Epitope.species")
  #printing out graphs 
  print(CDR3vsVDJDB_epitope)
  print(CDR3vsVDJDB_epitopegene)
  print(CDR3vsVDJDB_epitopespp)
  
  
  ##Saving it into data folder in box
  pdf("VDJdbEpitopeProportions_TconvVTreg[Karen].pdf", width=19, height=8)
  ggarrange(CDR3vsVDJDB_epitope, CDR3vsVDJDB_epitopegene, CDR3vsVDJDB_epitopespp, nrow = 1, ncol = 3) 
  dev.off()

###SKIP THI STEP 7/28/2023, KAREN MEETING = ISOLATING THE EXPANDED CELL TYPES in DBCDR3match3 ### 
  #ISOLATING JUNCTIONS THAT WERE EXPANDED# 
    #METHOD1: Find which junctions were expanded
      #tcrsP390 <- BRITCRs[ grepl("P390",BRITCRs$project),]
      #expandedJunctions <- unique( tcrsP390$junction[ duplicated( tcrsP390$full_nt_sequence)] )
  
    #METHOD2:  tcrsP390 are the tcrs that are just from P390
      #tcrsP390$sampleRepositoryID <- annoP390[ BRITCRs$libid, "sampleRepositoryID"]
      #tcrsP390$donorNt <- paste( tcrsP390$sampleRepositoryID, tcrsP390$full_nt_sequence)
      #expandedDonorNts <- tcrsP390$donorNt[ duplicated(tcrsP390$donorNt) ]
      #expandedJunctions2 <- unique( tcrsP390$junction[tcrsP390$donorNt %in% expandedDonorNts ] )
  
  #isolating libids that match this expanded junction list (eliminating rows that do not have a junction that match this expansion list)-->67% of data was expanded 
      #DBCDR3match3 <- DBCDR3match3[DBCDR3match3$junction %in% expandedJunctions, ]
  
      
#######- THIS IS WHERE I LEFT OFF: below is the for loop to calc the t.test/pvalues, you now need to make a df of df_filtered where you have rows that have columns listing the proportinos for each cell type, but they have to be in the same row so that they represent the same donor. Once you have this df, you can plug it into the for loop below and it will do a paired t-test. rmbr what makes it a PAIRED t-test is the fact that the proportions came from the same donor. Furthermore, the dimensinos of the groups needs to be the same in order to do a PAIREd t-test which is why this is important so we can exclude any of the data points for each plot that did not have a correlating dot in each cell type that came from the same patient. -#######
  
#Creating a new function to generate the dfs instead of making plots# 
  generateCDR3vsVDJDBdf <- function( epitopeColumnName = "Epitope"){
    #Numerator 
    ## Initialize empty vectors to store the results
    numeratorA <- vector()
    numeratorB <- vector()
    epitope <- vector()
    sampleIDcellType <- vector()
    studyGroup <- vector()
    cellType <- vector()
    sampleRepositoryID <- vector()
    ## Loop through unique srID-epitope combinations
    for (s in unique(DBCDR3match3$sampleIDcellType)) {
      for (e in unique(DBCDR3match3[DBCDR3match3$sampleIDcellType == s, epitopeColumnName])) {
        # Subset the dataframe for the current srID-epitope combination
        subset_df <- DBCDR3match3[DBCDR3match3$sampleIDcellType == s & DBCDR3match3[, epitopeColumnName] == e, ]
        # Calculate the numerator (count of unique libids)
        numA <- length(unique(subset(subset_df, Gene == "TRA")$libid))
        numB <- length(unique(subset(subset_df, Gene == "TRB")$libid))
        #num <- length(unique(subset_df$libid))
        # Store the results in the vectors
        numeratorA <- c(numeratorA, numA)
        numeratorB <- c(numeratorB, numB)
        epitope <- c(epitope, e)
        sampleIDcellType <- c(sampleIDcellType, s)
        studyGroup <- c(studyGroup, unique(subset_df$studyGroup))
        cellType <- c(cellType, unique(subset_df$cellType))
        sampleRepositoryID <- c(sampleRepositoryID, unique(subset_df$sampleRepositoryID))
      }
    }
    ## Create the new dataframe
    new_df <- data.frame(numeratorA, numeratorB, epitope, studyGroup, cellType, sampleIDcellType, sampleRepositoryID)
    ## Print the new dataframe
    print(new_df)
    
    
    #Denominator = total libids / project (non-VDJDB isolates)
    ##Rename project sampleRepositoryID in annoP390BRITCR2
    #annoP390BRITCRs3 <- annoP390BRITCRs2 %>%
    #rename(project = project.x)
    #Adding a column that includes the unique combos of sampleIDs + CellTyps 
    annoP390BRITCRs2$sampleIDcellType <- paste(annoP390BRITCRs2$sampleRepositoryID, annoP390BRITCRs2$cellType, sep = " ")
    ## Calculate the total number of libids for each project in df2 (annoP390BRITCRs3)
    libid_countsA <- annoP390BRITCRs2 %>%
      filter(chain == "TRA") %>%
      group_by(sampleIDcellType) %>%
      summarise(denominatorA = n())
    ## Merge libid_counts with new_df based on the project column
    new_df <- new_df %>%
      left_join(libid_countsA, by = "sampleIDcellType")
    ## Calculate the total number of libids for each project in df2 (annoP390BRITCRs3)
    libid_countsB <- annoP390BRITCRs2 %>%
      filter(chain == "TRB") %>%
      group_by(sampleIDcellType) %>%
      summarise(denominatorB = n())
    ## Merge libid_counts with new_df based on the project column
    new_df <- new_df %>%
      left_join(libid_countsB, by = "sampleIDcellType")
    ## Print the updated new_df dataframe
    print(new_df)
    
    
    #Calc proportions 
    new_df$proportionA = new_df$numeratorA / new_df$denominatorA
    new_df$proportionB = new_df$numeratorB / new_df$denominatorB
    
    
    #Further Filtering Data 
    ##converting all "Tscms" in cellType --> Tconvs 
    new_df$cellType <- gsub("Tscm", "Tconventional", new_df$cellType)
    ##converting all Tscms --> Tconvs for $sampleIDcellType 
    new_df$sampleIDcellType <- gsub("Tscm$", "Tconventional", new_df$sampleIDcellType)
    ##eliminating any epitope plots that do not represent both cellTypes ( we only want plots w/ proportion data for BOTH Tconvs + Tregs) 
    ### Step 1: Count the number of unique cellTypes for each epitope
    celltype_count <- aggregate(cellType ~ epitope, new_df, function(x) length(unique(x)))
    ### Step 2: Filter the dataframe to keep only data that had the same patient, but both cell types 
    filtered_epitopes <- celltype_count$epitope[celltype_count$cellType == 2]
    df_filtered <- new_df[new_df$epitope %in% filtered_epitopes, ]
    
    return(df_filtered)
  }
  

  #Applying the function you made to actually generate the 3 dfs of interest 
  CDR3vsVDJDB_epitopedf <- generateCDR3vsVDJDBdf("Epitope")
  CDR3vsVDJDB_epitopegenedf <- generateCDR3vsVDJDBdf("Epitope.gene")
  CDR3vsVDJDB_epitopesppdf <- generateCDR3vsVDJDBdf("Epitope.species")
  #printing out dfs 
  print(CDR3vsVDJDB_epitopedf)
  print(CDR3vsVDJDB_epitopegenedf)
  print(CDR3vsVDJDB_epitopesppdf)
  

  
  
#Pre-STATISTICS# - isolating only the unique epitope/sampleID combinations that have rows for BOTH Tregs & Tconvs. This way we can isolate only the proportions that come from the same donor, making the df now eligible to be run through paired t-test

# For a particular epitope, return only cells from donors that have both Tregs and Tconventional Cells
getTregTconvDonorsPerEpitope <- function(df, epitope){
  df <- df[ df$epitope == epitope,] # Subset df to contain just entries associated with epitope
  df$cellTypeSRID <- paste( df$cellType, df$sampleRepositoryID) 
  nodups <- df[ !duplicated(df$cellTypeSRID), ]
  srids <- nodups[ duplicated(nodups$sampleRepositoryID), "sampleRepositoryID"] # These SRIDs have both Tregs and Tconventional cells
  return( df[ df$sampleRepositoryID %in% srids, ] ) 
}
getTregTconvDonors <- function(df){
  epitopes <- unique(df$epitope)
  df <- do.call( rbind, lapply(epitopes, function(epitope) getTregTconvDonorsPerEpitope(df,epitope)))
  df <- df[order(df$epitope, df$sampleRepositoryID),]
  return(df)
}
  
  
  #Applying the function you made to actually generate the 3 dfs of interest (make new subset dfs for paired t test)
  CDR3vsVDJDB_epitopedfpaired <- getTregTconvDonors(CDR3vsVDJDB_epitopedf)
  CDR3vsVDJDB_epitopegenedfpaired <- getTregTconvDonors(CDR3vsVDJDB_epitopegenedf)
  CDR3vsVDJDB_epitopesppdfpaired <- getTregTconvDonors(CDR3vsVDJDB_epitopesppdf)
 ##Double Checking to convert all rows that Tscm 
  CDR3vsVDJDB_epitopedfpaired$cellType <- gsub("Tscm", "Tconventional", CDR3vsVDJDB_epitopedfpaired$cellType)
  CDR3vsVDJDB_epitopegenedfpaired$cellType <- gsub("Tscm", "Tconventional", CDR3vsVDJDB_epitopegenedfpaired$cellType)
  CDR3vsVDJDB_epitopesppdfpaired$cellType <- gsub("Tscm", "Tconventional", CDR3vsVDJDB_epitopesppdfpaired$cellType)
  ##Now, double ehcking/ fixing the sampleIDcellType column instead. 
  CDR3vsVDJDB_epitopedfpaired$sampleIDcellType <- gsub("Tscm$", "Tconventional", CDR3vsVDJDB_epitopedfpaired$sampleIDcellType)
  CDR3vsVDJDB_epitopegenedfpaired$sampleIDcellType <- gsub("Tscm$", "Tconventional", CDR3vsVDJDB_epitopegenedfpaired$sampleIDcellType)
  CDR3vsVDJDB_epitopesppdfpaired$sampleIDcellType <- gsub("Tscm$", "Tconventional", CDR3vsVDJDB_epitopesppdfpaired$sampleIDcellType)
  #checking: "notebook <- sum(grepl("Tscm$", CDR3vsVDJDB_epitopedfpaired$sampleIDcellType))"
  
#SUBSETTING DATA BY CHAIN TYPES 
  #Creating new dfs that delete the proportion column that is not related to said chain 
    CDR3vsVDJDB_epitopedfpairedA = CDR3vsVDJDB_epitopedfpaired[, !colnames(CDR3vsVDJDB_epitopedfpaired) %in% "proportionB"]
    CDR3vsVDJDB_epitopegenedfpairedA = CDR3vsVDJDB_epitopegenedfpaired[, !colnames(CDR3vsVDJDB_epitopegenedfpaired) %in% "proportionB"]
    CDR3vsVDJDB_epitopesppdfpairedA= CDR3vsVDJDB_epitopesppdfpaired[, !colnames(CDR3vsVDJDB_epitopesppdfpaired) %in% "proportionB"]
    CDR3vsVDJDB_epitopedfpairedB = CDR3vsVDJDB_epitopedfpaired[, !colnames(CDR3vsVDJDB_epitopedfpaired) %in% "proportionA"]
    CDR3vsVDJDB_epitopegenedfpairedB  = CDR3vsVDJDB_epitopegenedfpaired[, !colnames(CDR3vsVDJDB_epitopegenedfpaired) %in% "proportionA"]
    CDR3vsVDJDB_epitopesppdfpairedB = CDR3vsVDJDB_epitopesppdfpaired[, !colnames(CDR3vsVDJDB_epitopesppdfpaired) %in% "proportionA"]
  #Renaming proportionA or proportionB with just "proportion" so it will work with the fxn below for stats testing 
    colnames(CDR3vsVDJDB_epitopedfpairedA)[colnames(CDR3vsVDJDB_epitopedfpairedA) == "proportionA"] <- "proportion"
    colnames(CDR3vsVDJDB_epitopegenedfpairedA)[colnames(CDR3vsVDJDB_epitopegenedfpairedA) == "proportionA"] <- "proportion"
    colnames(CDR3vsVDJDB_epitopesppdfpairedA)[colnames(CDR3vsVDJDB_epitopesppdfpairedA) == "proportionA"] <- "proportion"
    colnames(CDR3vsVDJDB_epitopedfpairedB)[colnames(CDR3vsVDJDB_epitopedfpairedB) == "proportionB"] <- "proportion"
    colnames(CDR3vsVDJDB_epitopegenedfpairedB)[colnames(CDR3vsVDJDB_epitopegenedfpairedB) == "proportionB"] <- "proportion"
    colnames(CDR3vsVDJDB_epitopesppdfpairedB)[colnames(CDR3vsVDJDB_epitopesppdfpairedB) == "proportionB"] <- "proportion"
    
#STATISTICS# - calculate paired t-test on each of the epitope data sets
  pairedttest <- function(df){
  
  ###Have to turn out proportions --> normalized distribution form ["normalizing our data"] 
  logitTransform <- function(p){
    p[ p== 0] <- 0.001
    return( log(p/(1-p)))
  }
  
    # Get unique epitope types
  unique_epitopes <- unique(df$epitope)
  
  # Create an empty list to store the t-test results
  t_test_results <- list()
  
  #df <- df[ order(df$sampleRepositoryID)]
  df <- df %>% arrange(sampleRepositoryID)
  
  # Loop over each unique epitope
  for (epitope in unique_epitopes) {
    print(epitope)
    # Subset the data for the current epitope and cellTypes A and B
    stats_df <- df[df$epitope == epitope & (df$cellType == "Treg" | df$cellType == "Tconventional"), ]
    
    # Extract the proportions for cellType A and B
    proportions_A <- stats_df[stats_df$cellType == "Treg", "proportion"]
    proportions_B <- stats_df[stats_df$cellType == "Tconventional", "proportion"]
    ##only include subsets that have 3+ data points for proportions 
    if( length(proportions_A) < 3 | length(proportions_B) < 3 ){
      print("Not enough data")
      t_test_results[[epitope]] <- 1
      next
    }
    
    # Perform t-test
    t_test_result <- t.test( logitTransform(proportions_A), logitTransform(proportions_B), paired=TRUE)
    
    # Store the t-test result in the list
    t_test_results[[epitope]] <- t_test_result$p.value
  }
  

  return(t_test_results)
  }
  
  ##Applying stats function to get actual t test info. 
  
  #Applying the function you made to actually generate the 3 dfs of interest (make new subset dfs for paired t test)
  CDR3vsVDJDB_epitopettestA <- pairedttest(CDR3vsVDJDB_epitopedfpairedA)
  CDR3vsVDJDB_epitopegenettestA <- pairedttest(CDR3vsVDJDB_epitopegenedfpairedA)
  CDR3vsVDJDB_epitopesppttestA <- pairedttest(CDR3vsVDJDB_epitopesppdfpairedA)
  CDR3vsVDJDB_epitopettestB <- pairedttest(CDR3vsVDJDB_epitopedfpairedB)
  CDR3vsVDJDB_epitopegenettestB <- pairedttest(CDR3vsVDJDB_epitopegenedfpairedB)
  CDR3vsVDJDB_epitopesppttestB <- pairedttest(CDR3vsVDJDB_epitopesppdfpairedB)
  
  
#REPLOTTING only the graphs for PAIRED T-TEST + attaching their correlating p-values. 
  generatepairedttestplots <- function( df_epitopes = CDR3vsVDJDB_epitopedfpaired, epitope_pvalues = CDR3vsVDJDB_epitopettest){
    # Iterate over each epitope name and create the corresponding plot
    for (epitope in unique(df_epitopes$epitope)) {
      # Subset the data for the current epitope
      epitope_data <- subset(df_epitopes, epitope == epitope)
      
      epitopes <- names(epitope_pvalues)
      epitope_pvalues <- as.numeric(epitope_pvalues); names(epitope_pvalues) <- epitopes
      
      toplot <- df_epitopes[!is.na(df_epitopes$epitope),]
      toplot$epitopePval <- paste( toplot$epitope, "\np-value: ", round(epitope_pvalues[toplot$epitope], 4) )
      
      # Create the plot for the current epitope
      testplot <- ggplot(toplot, aes(x=cellType,y=proportion, color=studyGroup, group=sampleRepositoryID)) + geom_boxplot(aes(x=cellType, group=cellType),color="black") + geom_line() + geom_point() + facet_wrap(~epitopePval)
      
      # Get the corresponding p-value from the list
      #p_value <- epitope_pvalues[epitope]
      
      # Add the p-value annotation to the plot
      #testplot <- testplot +
      #annotate("text", x = x_position, y = y_position, label = paste("p-value =", p_value))
      
      # Print or save the plot
      return(testplot)
      # or save the plot using ggsave()
    }
  }
  
   #running the function generate the 3 plots. 
  CDR3vsVDJDB_epitopestatA = generatepairedttestplots (CDR3vsVDJDB_epitopedfpairedA, CDR3vsVDJDB_epitopettestA)
  CDR3vsVDJDB_epitopegenestatA = generatepairedttestplots (CDR3vsVDJDB_epitopegenedfpairedA, CDR3vsVDJDB_epitopegenettestA)
  CDR3vsVDJDB_epitopesppstatA = generatepairedttestplots (CDR3vsVDJDB_epitopesppdfpairedA, CDR3vsVDJDB_epitopesppttestA)
  CDR3vsVDJDB_epitopestatB = generatepairedttestplots (CDR3vsVDJDB_epitopedfpairedB, CDR3vsVDJDB_epitopettestB)
  CDR3vsVDJDB_epitopegenestatB = generatepairedttestplots (CDR3vsVDJDB_epitopegenedfpairedB, CDR3vsVDJDB_epitopegenettestB)
  CDR3vsVDJDB_epitopesppstatB = generatepairedttestplots (CDR3vsVDJDB_epitopesppdfpairedB, CDR3vsVDJDB_epitopesppttestB)
  
  CDR3vsVDJDB_epitopestatAt <- CDR3vsVDJDB_epitopestatA + labs(title = "TRA Proportions/Epitope Seq.")
  CDR3vsVDJDB_epitopegenestatAt <- CDR3vsVDJDB_epitopegenestatA + labs(title = "TRA Proportions/Epitope Gene")
  CDR3vsVDJDB_epitopesppstatAt <- CDR3vsVDJDB_epitopesppstatA + labs(title = "TRA Proportions/Epitope Spp.")
  CDR3vsVDJDB_epitopestatBt <- CDR3vsVDJDB_epitopestatB + labs(title = "TRB Proportions/Epitope Seq.")
  CDR3vsVDJDB_epitopegenestatBt <- CDR3vsVDJDB_epitopegenestatB + labs(title = "TRB Proportions/Epitope Gene")
  CDR3vsVDJDB_epitopesppstatBt <- CDR3vsVDJDB_epitopesppstatB + labs(title = "TRB Proportions/Epitope Spp.")
  
  print(CDR3vsVDJDB_epitopestatA)
  print(CDR3vsVDJDB_epitopegenestatA)
  print(CDR3vsVDJDB_epitopesppstatA)
  print(CDR3vsVDJDB_epitopestatB)
  print(CDR3vsVDJDB_epitopegenestatB)
  print(CDR3vsVDJDB_epitopesppstatB)
  
  ##Saving it into data folder in box
  pdf("VDJdbEpitopeProportions_TconvVTreg[pairedt-test](KAREN)(2).pdf", width=19, height=8)
  ggarrange(CDR3vsVDJDB_epitopestatAt, CDR3vsVDJDB_epitopegenestatAt, CDR3vsVDJDB_epitopesppstatAt,CDR3vsVDJDB_epitopestatBt, CDR3vsVDJDB_epitopegenestatBt, CDR3vsVDJDB_epitopesppstatBt, nrow = 2, ncol = 3) 
  dev.off()
  
  pdf("VDJdbEpitopeProportions_TconvVTreg[pairedt-test](KAREN)no-title(2).pdf", width=19, height=8)
  ggarrange(CDR3vsVDJDB_epitopestatA, CDR3vsVDJDB_epitopegenestatA, CDR3vsVDJDB_epitopesppstatA,CDR3vsVDJDB_epitopestatB, CDR3vsVDJDB_epitopegenestatB, CDR3vsVDJDB_epitopesppstatB, nrow = 2, ncol = 3) 
  dev.off()
  

  #RESAVING ALL FIGURES AS JPEGs 
  ggsave("/Users/klai/Box/2023SummerInternKenneth/data/official_figures/epitopeseqproportions[Karen](2).jpg", plot = CDR3vsVDJDB_epitope, width = 7.5, height = 7.5)
  ggsave("/Users/klai/Box/2023SummerInternKenneth/data/official_figures/epitopegeneproportions[Karen](2).jpg", plot = CDR3vsVDJDB_epitopegene, width = 7.5, height = 7.5)
  ggsave("/Users/klai/Box/2023SummerInternKenneth/data/official_figures/epitopesppproportions[Karen](2).jpg", plot = CDR3vsVDJDB_epitopespp, width = 7.5, height = 7.5)
  ggsave("/Users/klai/Box/2023SummerInternKenneth/data/official_figures/epitopeseqproportionsstatA[Karen](2).jpg", plot = CDR3vsVDJDB_epitopestatAt, width = 7.5, height = 7.5)
  ggsave("/Users/klai/Box/2023SummerInternKenneth/data/official_figures/epitopegeneproportionsstatA[Karen](2).jpg", plot = CDR3vsVDJDB_epitopegenestatAt, width = 7.5, height = 7.5)
  ggsave("/Users/klai/Box/2023SummerInternKenneth/data/official_figures/epitopesppproportionsstatA[Karen](2).jpg", plot = CDR3vsVDJDB_epitopesppstatAt, width = 7.5, height = 7.5)
  ggsave("/Users/klai/Box/2023SummerInternKenneth/data/official_figures/epitopeseqproportionsstatB[Karen](2).jpg", plot = CDR3vsVDJDB_epitopestatBt, width = 7.5, height = 7.5)
  ggsave("/Users/klai/Box/2023SummerInternKenneth/data/official_figures/epitopegeneproportionsstatB[Karen](2).jpg", plot = CDR3vsVDJDB_epitopegenestatBt, width = 7.5, height = 7.5)
  ggsave("/Users/klai/Box/2023SummerInternKenneth/data/official_figures/epitopesppproportionsstatB[Karen](2).jpg", plot = CDR3vsVDJDB_epitopesppstatBt, width = 7.5, height = 7.5)

  
  
  
  
  
###__________________SCRAP__________________###

  
