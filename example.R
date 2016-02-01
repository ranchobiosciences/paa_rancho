# Rancho PAA example script
# John Obenauer, Ph.D., and Ivan Grishagin, Ph.D.
# January 2016
# 
# This R script shows an example of how to analyze Invitrogen 
# ProtoArrays using the modifications to the PAA package 
# made by Rancho BioSciences.  The example data set used is 
# Mark Gerstein's data, consisting of positive control samples 
# and a normal sample.  In a typical R environment, this script 
# can be run using:
#   R --no-save < gerstein.R

# Install R and Bioconductor packages, if needed. 
# These are commented out by default, but if the library 
# commands below give an error, then uncomment these lines 
# to make sure the dependencies are installed.
#install.packages("gplots")
#install.packages("gtools")
#install.packages("openxlsx")
#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("MASS")

# Load required R and Bioconductor packages
library(gplots)
library(gtools)
library(limma)
library(MASS)
library(openxlsx)

# Load Rancho-modified PAA package
source("PAA_Rancho_v21.R")

# Set path for input *.gpr files and targets file (sample information file)
# Change these if your GPR files and targets file are in another location
gpr.path = "gpr"
targets.path = "gpr/targets.txt"

# Read in the targets file (sample information file).  
targets = readTargets(targets.path)

# Load GPR files and make expression list (elist) object
elist.raw = loadGPR(gpr.path=gpr.path, targets.path=targets.path)

# Select control spots on the array to use in normalization.  
# By default, Invitrogen ProtoArrays use these two sets of 
# control spots, HumanIgG and Anti-HumanIgG.  There is an 
# option to add a V5 spike-in control, which improves 
# normalization.  The Gerstein data uses the spike-in control, 
# so this is indicated in the command below.
controls.names<-c("HumanIgG","Anti-HumanIgG", "V5_")
# Many experiments will need this command instead
#controls.names<-c("HumanIgG","Anti-HumanIgG")

# Make a list to store the control spots that will be used
controls.elist = elist.raw

# Create an output directory for results
output.path<-paste0(getwd(),"/output")
dir.create(output.path, showWarnings = FALSE)

# If the ProtoArrays were run in several batches, a two-step batch 
# correction process is recommended.  A "Batch" column in the targets 
# file is used to distinguish different batches.  Its values can be 
# anything (like "1" and "2", "A" and "B", "2016-05-01" and "2016-05-15").  
# Although two batches will be shown in this example, any number of 
# batches can be present in the data.  

# The Gerstein data has all the samples run separately with Hi gain and 
# Lo gain.  We will treat these Hi and Lo samples as separate batches.  

# Check whether batches are present in the data by looking for multiple 
# unique values in the Batch column of the targets file
if (length(unique(elist.raw$targets$Batch)) > 1) {

    # Get unique group names
    targ.groups = unique(elist.raw$targets$Group)

    # Initialize empty list of discarded rows (spots on the arrays)
    all.disc.rows<-NULL
        
    # For each group found in the data (such as "Disease" and "Healthy", 
    # identify values that differ too much according to batch (not group) 
    # and remove these from the analyses.  This prevents false positives.
    for (gr in targ.groups){
        gr.IDs = elist.raw$targets[elist.raw$targets$Group == gr,"ArrayID"]
        elist = elist.raw[,gr.IDs]
        
        # Find rows that need to be discarded for the current group
        disc.row.num = batchFilter.anova(elist=elist,  
            p.thresh=0.0001, fold.thresh=5, output.path=output.path)
                                                
        # Keep track of discarded rows for all groups
        all.disc.rows = c(all.disc.rows,disc.row.num)
            
    }
    
    # If any rows were filtered out, create a new elist with these rows removed
    if (!is.null(all.disc.rows)){
            # Make a subset of the original elist with the discarded rows removed
            elist.filter = elist.raw[-unique(all.disc.rows),]
    } else {
            message("batchFilter: nothing to filter!")
            elist.filter = elist.raw
    }
}

# Write out the total number of unique spots filtered out
message(paste0("Number of unique features discarded: ", length(unique(all.disc.rows))))

# If batch correction is not needed for a data set, "elist.raw" can be substituted for 
# "elist.filter" in subsequent commands.  Alternatively, the line "elist.filter = elist.raw" 
# can be added to avoid altering the subsequent lines.

# Perform background correction
elist.filter = backgroundCorrect(elist.filter, method="minimum")
controls.elist = backgroundCorrect(controls.elist, method="minimum")

# Fill in values of control spots to be used for normalization, and 
# remove these spots from the other elists
controls.elist = make.controls.elist(elist=controls.elist, controls.names=controls.names)
elist.filter = remove.control.spots(elist=elist.filter,controls.names=controls.names)

# RLM normalization
elist.norm = normalizeArrays(method="rlm", elist=elist.filter, controls.elist=controls.elist,
    gpr.path=gpr.path, targets.path=targets.path, output.path=output.path, contr.names=NULL)

# Batch correction, second step
# In our data sets with multiple batches, we found the best results were obtained when both 
# batch filtering (done earlier) and batch correction are used.  The filtering step only removes 
# outlier values that are batch associated; this correction step modifies the values up or down 
# to account for the batch effect.  
if (length(unique(elist.raw$targets$Batch)) > 1) {
    elist.norm = batchAdjust(elist=elist.norm, is.logged=TRUE)
}

# Write out normalized and batch corrected data, in logged and unlogged forms
elist.norm.output = cbind.data.frame(Description = elist.norm$genes[,"Name"], elist.norm$E)
write.table(elist.norm.output, paste0(output.path, "/gerstein_normalized_log.txt"), 
    sep="\t", row.names=FALSE, col.names=TRUE)
elist.norm.output.unlog = elist.norm.output
elist.norm.output.unlog[,2:ncol(elist.norm.output.unlog)] = 2^elist.norm$E
write.table(elist.norm.output.unlog, paste0(output.path, "/gerstein_normalized_unlog.txt"), 
    sep="\t", row.names=FALSE, col.names=TRUE)

# Set desired thresholds for fold change and adjusted p-value
foldChange=1.5
pvalue=0.05

# Create table to map Block/Row/Column values to protein names for output files 
brc2name = new.env()
for (i in 1:nrow(elist.raw)) {
    brc2name[[ as.character(paste(elist.raw$genes$Block[i], elist.raw$genes$Row[i], 
        elist.raw$genes$Column[i], sep=",")) ]] = elist.raw$genes$Name[i]
}

# Compare the experimental samples to controls
experiment = "PCS100"
control = "PCS000"
output.path.expt = paste0(output.path, "/", experiment);
dir.create(output.path.expt, showWarnings = FALSE)
expt.ids = targets$ArrayID[targets$Group == experiment]
cont.ids = targets$ArrayID[targets$Group == control]
expt.idx = match(expt.ids, colnames(elist.norm))
cont.idx = match(cont.ids, colnames(elist.norm))

# Unlog values for volcano plots
elist.norm.unlog = elist.norm
elist.norm.unlog$E = 2^elist.norm.unlog$E
volcanoPlot(elist=elist.norm.unlog, group1=expt.idx, group2=cont.idx, method = "mMs",
    p.thresh = pvalue, fold.thresh = foldChange, output.path = output.path.expt,
    tag = "_mM_", above=500, between=200, title="Positive Controls vs. Normal (M-statistic)")
volcanoPlot(elist=elist.norm.unlog, group1=expt.idx, group2=cont.idx, method = "tTest",
    p.thresh = pvalue, fold.thresh = foldChange, output.path = output.path.expt,
    tag = "_tT_", above=500, between=200, title="Positive Controls vs. Normal (t-Test)")

# Define row names (gene info) and column names (sample info)
rownames(elist.norm.unlog) = paste0(elist.norm.unlog$genes$Block, ",", elist.norm.unlog$genes$Row, ",",  
    elist.norm.unlog$genes$Column)
colnames(elist.norm.unlog) = elist.norm.unlog$targets$SampleName
# Alternatively, this line labels each colum with both Group name and SampleName
#colnames(elist.norm.unlog) = paste0(elist.norm.unlog$targets$Group, " ", elist.norm.unlog$targets$SampleName)

# Perform differential analysis
diffAnalysis.result = diffAnalysis(input=elist.norm.unlog$E, label1=expt.idx, label2=cont.idx,
    class1=experiment, class2=control, sortBY="fold_change", above=500, between=200,
    output.path=output.path.expt, filename = paste0("diffAnalysis_mMs_", experiment, ".txt"))

# Add protein names to differential expression results
Protein = as.character(unlist(mget(diffAnalysis.result$BRC, brc2name, ifnotfound=NA)))
diffAnalysis.result = cbind(Protein, diffAnalysis.result)

# Select features based on t-test
selected.tTest.list = diffAnalysis.result[diffAnalysis.result$tTest < pvalue & 
    (diffAnalysis.result$fold_change > foldChange | diffAnalysis.result$fold_change < 1 / foldChange),]
selected.tTest = diffAnalysis.result$BRC[diffAnalysis.result$tTest<pvalue & 
    (diffAnalysis.result$fold_change > foldChange | diffAnalysis.result$fold_change < 1 / foldChange)]
if(length(selected.tTest)){

    # Write selected results to file
    openxlsx:::write.xlsx(selected.tTest.list, file = paste0(output.path.expt, "/results-tTest-FC_",
        foldChange, "-pval_", pvalue, "_", experiment, ".xlsx"))
    
    # Draw heatmap
    selected.tTest.idx = match(selected.tTest, rownames(elist.norm.unlog))
    elist.heatmap = elist.norm.unlog[selected.tTest.idx, c(expt.idx, cont.idx)]
    #title=paste0(experiment, " vs. ", control) # this would typically be used
    title="Positive Controls vs. Normal (t-Test)"    # overwrite default title for Gerstein data
    if(length(selected.tTest) > 100) {
            cexRowSize = 0.1
    } else if(length(selected.tTest) > 40) {
            cexRowSize = 0.4
    } else {
            cexRowSize = 0.8
    }
    ColSideColors<-c(rep("firebrick",length(expt.idx)), rep("green4",length(cont.idx)))
    heatmap.2.mod(elist.heatmap$E,             
        Rowv = TRUE,                             #row and col clustering
        Colv = TRUE,
        ColSideColors=ColSideColors,
        col=blueyellow(300),                   #color scheme
        scale="row",                           #scale by column
        main=title,                            #title of the heatmap
        margins = c(15,15),                    #margins around heatmap    
        key=TRUE,                              #legend is present
        key.title="",                          #no title for legend
        key.xlab="",                           #no label for legend axis
        key.ylab="",
        keysize=1.1,                           #size of legend strip
        symkey=FALSE,                          #do not make the colors symmetrical around 0
        density.info="density",                   #no density information
        trace="none",                          #
        cexRow=cexRowSize,                     #column labels' size
        labRow=row.names(elist.heatmap),       #row labels
        labCol=colnames(elist.heatmap)         #column labels: limit to 30 chars
    )
        
}

# Select features based on M-statistic
selected.mMs.list = diffAnalysis.result[diffAnalysis.result$mMs < pvalue & 
    (diffAnalysis.result$fold_change > foldChange | diffAnalysis.result$fold_change < 1 / foldChange),]
selected.mMs = diffAnalysis.result$BRC[diffAnalysis.result$mMs < pvalue & 
    (diffAnalysis.result$fold_change > foldChange | diffAnalysis.result$fold_change < 1 / foldChange)]
if(length(selected.mMs)){

    # Write selected results to file
    openxlsx:::write.xlsx(selected.mMs.list, file = paste0(output.path.expt, "/results-mMs-FC_",
        foldChange, "-pval_", pvalue, "_", experiment, ".xlsx"))
    
    # Draw heatmap
    selected.mMs.idx = match(selected.mMs, rownames(elist.norm.unlog))
    elist.heatmap = elist.norm.unlog[selected.mMs.idx, c(expt.idx, cont.idx)]
    #title=paste0(experiment, " vs. ", control, " (M-statistic)")
    title="Positive Controls vs. Normal (M-statistic)"
    if(length(selected.mMs) > 100) {
        cexRowSize = 0.1
    } else if(length(selected.mMs) > 40) {
        cexRowSize = 0.4
    } else {
        cexRowSize = 0.8
    }
    ColSideColors<-c(rep("firebrick",length(expt.idx)), rep("green4",length(cont.idx)))
    heatmap.2.mod(elist.heatmap$E,             
        Rowv = TRUE,                             #row and col clustering
        Colv = TRUE,
        ColSideColors=ColSideColors,
        col=blueyellow(300),                   #color scheme
        scale="row",                           #scale by column
        main=title,                            #title of the heatmap
        margins = c(15,15),                    #margins around heatmap    
        key=TRUE,                              #legend is present
        key.title="",                          #no title for legend
        key.xlab="",                           #no label for legend axis
        key.ylab="",
        keysize=1.1,                           #size of legend strip
        symkey=FALSE,                          #do not make the colors symmetrical around 0
        density.info="density",                   #no density information
        trace="none",                          #
        cexRow=cexRowSize,                     #column labels' size
        labRow=row.names(elist.heatmap),       #row labels
        labCol=colnames(elist.heatmap)         #column labels: limit to 30 chars
    )
        
}
                
	
