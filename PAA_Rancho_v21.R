##########################################################################
#           ProteinArrayAnalyzer (PAA) - Rancho                          #
##########################################################################

##########################################################################
#                            FUNCTIONS                                   #
##########################################################################

#+++++++++++++++++++++++++++ loadGPR +++++++++++++++++++++++++++++++++++++++
loadGPR <- function(gpr.path=NULL, 
                    targets.path=NULL, 
                    array.type="ProtoArray",
                    protoarray.aggregation="min", 
                    array.columns=list(E="F635 Median",
                                       Eb="B635 Median"),
                    array.annotation=c("Block", 
                                       "Column",
                                       "Row", 
                                       "Description", 
                                       "Name",
                                       "ID")){
        
        if(is.null(gpr.path) || 
           is.null(targets.path)) {
                stop("ERROR: Not all mandatory arguments have been defined!")
        }
        targets <- readTargets(targets.path)
        elist <- read.maimages(files=targets,
                               path=gpr.path,
                               source="genepix.median", 
                               columns=array.columns,
                               annotation=array.annotation)
        
        # If no match, grep will return 'integer(0)' resulting in an empty elist$E
        # in the following lines. Hence the following if-statements are necessary:
        if(any(grep("Empty", elist$genes$Description))){  
                elist <- elist[-grep("Empty", elist$genes$Description),]
        }
        
        colnames(elist$E) <- targets$ArrayID
        colnames(elist$Eb) <- colnames(elist$E)
        
        if(array.type=="ProtoArray" && protoarray.aggregation=="min"){
                row.len <- nrow(elist$E)
                col.len <- ncol(elist$E)
                tmp.col.len <- (row.len*col.len)/2
                elist$E[row(elist$E)[,1]%%2==1] <-
                        matrix(apply(matrix(elist$E,2,tmp.col.len),2,min),
                               row.len/2,col.len)
                elist <- elist[-row(elist)[,1]%%2==1,]
                
        }else if(array.type=="ProtoArray" && protoarray.aggregation=="mean"){
                elist$E[row(elist$E)[,1]%%2==1,] <-
                        (elist$E[row(elist$E)[,1]%%2==1,] +
                                 elist$E[row(elist$E)[,1]%%2==0,])/2
                
                elist <- elist[-row(elist)[,1]%%2==1,]
        }
        return(elist)
}
#+++++++++++++++++++++++++++ loadGPR +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ remove.control.spots +++++++++++++++++++++++++++++++++++++++
remove.control.spots<-function(elist=NULL,
                               controls.names=NULL){
        #if no elist - abort
        if(is.null(elist)) {
                stop("ERROR: No elist supplied!")
        }
        
        if(is.null(controls.names)){
                #remove control spots if based on description 
                if(any(grep("Control", elist$genes$Description))){
                        elist <- elist[-grep("Control", elist$genes$Description),]
                        cat("Control spots removed from main elist.\n")
                } else {
                        warning("No rows with description \"Control\" found and no controls' names were supplied!")
                }
        } else {
                #look up control spots based on the supplied names
                #empty vector to stor row indices of control spots
                controlsRowsVector<-NULL
                #for each control name look it up in elist
                for (contrName in controls.names){
                        controlsRowsVector<-c(controlsRowsVector,grep(contrName, elist$genes$Name))
                }
                #remove controls elist, if found
                if(length(controlsRowsVector)==0){
                        warning("No controls from the supplied list were found")
                }
                elist<-elist[-controlsRowsVector,]
                
                cat(c(length(controlsRowsVector),
                      "control spots based on supplied controls' names were removed.\n"))
        } 
        return(elist)
}
#+++++++++++++++++++++++++++ remove.control.spots +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ make.controls.elist +++++++++++++++++++++++++++++++++++++++
make.controls.elist<-function(elist=NULL,
                              controls.names=NULL){
        #if no elist or controls.names - abort
        if(is.null(elist) || is.null(controls.names)) {
                stop("ERROR: No elist and/or controls names supplied!")
        }
        
        #empty vector to stor row indices of control spots
        controlsRowsVector<-NULL
        #for each control name look it up in elist
        #and add to vector of row indices
        for (contrName in controls.names){
                controlsRowsVector<-c(controlsRowsVector,
                                      grep(contrName, elist$genes$Name))
        }
        #if not hits, abort
        if (length(controlsRowsVector)==0) stop("List of controls is empty!\nCheck names!")
        #otherwise, prepare controls elist
        #also, only consider unique controls' indices, otherwise they may be considered more than once
        controls.elist<-elist[unique(controlsRowsVector),]
        
        return(controls.elist)
}
#+++++++++++++++++++++++++++ make.controls.elist +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ batchFilter +++++++++++++++++++++++++++++++++++++++
batchFilter <- function(elist=NULL, 
                        lot1=NULL, 
                        lot2=NULL, 
                        p.thresh=0.05,
                        fold.thresh=1.5, 
                        output.path=NULL){
    
    if(is.null(elist) || 
       is.null(lot1) || 
       is.null(lot2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    row.number <- nrow(elist$E)
    col.number <- ncol(elist$E)
    discard <- c()
    output <- matrix(nrow=row.number, ncol=4)
    colnames(output) <- c("BRC", 
                          "P-value", 
                          "Fold change", 
                          "Description")

    for(zeile in 1:row.number) {
        p <- try(t.test(x=elist$E[zeile,lot1],
                        y=elist$E[zeile,lot2])$p.value, 
                 FALSE)
        if(!is.numeric(p)){p <- 1}
        
        fold <- mean(elist$E[zeile, lot1])/mean(elist$E[zeile, lot2])
        
        output[zeile,] <- c(paste(elist$genes[zeile,1], " ",
                                  elist$genes[zeile,3], " ", 
                                  elist$genes[zeile,2], sep=""), 
                            p, 
                            fold,
                            elist$genes[zeile,4])
        
        if(p < p.thresh && (fold > fold.thresh || fold < 1/fold.thresh)) {
                #---> ISSUE: there  may be no exact representation for floats!!!
                discard <- c(discard,zeile)
        }
    }
    message(paste("batchFilter - number of features to discard: ",
                  length(discard), sep=""), "\n")
    
    if(is.null(output.path)) {
        plot(x=log2(as.numeric(output[,3])), 
             y=-log10(as.numeric(output[,2])),
             xlab="log2(fold change)", 
             ylab="-log10(p-value)",
             main="batch filter volcano")
            
        if(length(discard)>0){
                points(x=log2(as.numeric(output[discard,3])),
                       y=-log10(as.numeric(output[discard,2])), 
                       col="red")
        }
            
        abline(h=-log10(p.thresh), lty=2, col="red")
        abline(v=log2(fold.thresh), lty=2, col="red")
        abline(v=log2(1/fold.thresh), lty=2, col="red")
        if(length(discard)>0){
            elist <- elist[-discard,]
        }
    }else{
        write.table(x=output, 
                    file=paste(output.path,
          "/batch_filter.txt",sep=""), 
          sep="\t", eol="\n",
          row.names=FALSE,
          quote=FALSE)
        write.table(x=output[discard,], 
                    file=paste(output.path,
          "/batch_filter_discarded.txt",sep=""), sep="\t", eol="\n",
          row.names=FALSE, quote=FALSE)
        
        tiff(paste(output.path,"/batch_filter.tiff",sep=""), 
             width=2000,
             height=2000,
             pointsize=15, 
             compression="lzw", 
             res=300)
        
            plot(x=log2(as.numeric(output[,3])),
                 y=-log10(as.numeric(output[,2])),
                 xlab="log2(fold change)",
                 ylab="-log10(p-value)", 
                 main="batch filter volcano")
            
            if(length(discard)>0){
                points(x=log2(as.numeric(output[discard,3])),
                       y=-log10(as.numeric(output[discard,2])), 
                       col="red")
            }
            abline(h=-log10(p.thresh), lty=2, col="red")
            abline(v=log2(fold.thresh), lty=2, col="red")
            abline(v=log2(1/fold.thresh), lty=2, col="red")
        dev.off()
        
        if(length(discard)>0){
            tiff(paste(output.path,"/batch_filter_discarded.tiff",sep=""),
                 width=2000,
                 height=2000, 
                 pointsize=15,
                 compression="lzw",
                 res=300)
                
                plot(x=log2(as.numeric(output[-discard,3])),
                     y=-log10(as.numeric(output[-discard,2])),
                     xlab="log2(fold change)",
                     ylab="-log10(p-value)",
                     main="batch filter volcano")
            dev.off()
            
            #elist <- elist[-discard,]
        }
    }
    return(discard)
}
#+++++++++++++++++++++++++++ batchFilter +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ batchAdjust +++++++++++++++++++++++++++++++++++++++
batchAdjust <- function(elist=NULL, 
                        is.logged=NULL){
        if(is.null(elist) || 
           is.null(is.logged)) {
                stop("ERROR: Not all mandatory arguments have been defined!")
        }
        if(is.logged){
                elist$E <-
                        ComBat(dat=elist$E,
                               mod=factor(elist$targets$Group),
                               batch=factor(elist$targets$Batch))
        }else if(!is.logged){
                elist$E <-
                        2^ComBat(dat=log2(elist$E),
                                 mod=factor(elist$targets$Group),
                                 batch=factor(elist$targets$Batch))
        }else{
                stop("ERROR in batchAdjust: is.logged not 'TRUE' or 'FALSE'.")
        }
        return(elist)
}
#+++++++++++++++++++++++++++ batchAdjust +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++ normalize.elist +++++++++++++++++++++++++++++++++++++
# Function added by Rancho.  Takes elist and rlm results as input and returns 
# a normalized elist.  Logarithmic intensity of each spot is corrected according 
# to y' = y-coef1-coef2
#where  y - raw logarithmic intensity 
#       coef1 - interslide coefficient (between arrays)
#       coef2 - intraslide coefficient (between blocks)
normalize.elist<-function(elist=elist,
                          n.arrays=n.arrays,
                          n.blocks=n.blocks,
                          rlm.result=rlm.result){
        # store logged readout values in a new matrix
        normalizedNew <- log2(elist$E)
        
        # merge all coefficients pertaining to the effect
        # of arrays into one vector
        coef1vect<-c(sapply(2:(1+n.arrays),
                            function(arrNum){
                                    rep(rlm.result$coefficients[arrNum],
                                        nrow(elist$E))
                            }))
        # vector -> matrix
        coef1M<-matrix(coef1vect,ncol=n.arrays)
        
        # record how many spots there are in each block 
        blockTable<-as.numeric(table(elist$genes$Block))
        
        # merge all coefficients pertaining to the effect
        # of blocks into one vector
        coef2vect<-unlist(sapply(1:n.blocks,
                                 function(blockNum){
                                         numReps<-blockTable[blockNum]
                                         vect<-rep(rlm.result$coefficients[(1+n.arrays+blockNum)],
                                                   numReps)
                                         return(vect)
                                 }))
        # vector -> matrix
        coef2M<-matrix(coef2vect,
                       nrow(elist$E),
                       ncol=n.arrays)
        # subtract coefficients from the original logged values
        normalizedNew<-log2(elist$E)-coef1M-coef2M
        return(normalizedNew)
}
#+++++++++++++++++++++++ normalize.elist +++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++ normalizeRLM +++++++++++++++++++++++++++++++++++++++
normalizeRLM <- function(elist=NULL, controls.elist=NULL, gpr.path=NULL,
                         targets.path=NULL, contr.names=NULL, output.path=NULL){    
        cat("Running Rancho-modified normalizeRLM...\n")
        # Check if controls.elist was supplied,
        # if not - create
        if(is.null(controls.elist)){
                controls.elist <-
                        loadGPR(gpr.path=gpr.path, targets.path=targets.path)
                        
                # Only use the "Control" value for chips that have them (more recent ones)
                # (IvanG) wrong elist object: need to use controls.elist instead of elist. Fixed.
                
                if (any(grep("Control", controls.elist$genes$Description))) {
                        controls.elist <- 
                                controls.elist[grep("Control", controls.elist$genes$Description),]
                }
                # (IvanG) include only HumanIgG and Anti-HumanIgG spots (by default).
                controls.elist <-
                        controls.elist[grep("HumanIg", controls.elist$genes$Name),]
        }
        
        # (IvanG) form a list of control spot names if they were not supplied 
        # each control name is a unique identifier formed as follows:
        # [protein spot name]-r[row of the spot in a block]c[column of the spot in a block]
        # this way identifier is completely unique and even if the protein is present
        # more than once per block, it will be accounted for correctly
        if(is.null(contr.names)){
                #remove "~N/A" from controls' names (if at all present)
                controls.elist$genes$Name<-unlist(strsplit(controls.elist$genes$Name,split = "~N/A"))
                #add a unique identifier: modify original controls.elist
                controls.elist$genes$Name<-paste0(controls.elist$genes$Name,
                                                  "-c",controls.elist$genes$Column,
                                                  "r",controls.elist$genes$Row)
                #get a list of unique controls' names
                contr.names<-unique(controls.elist$genes$Name)
        }
       
        n.controls <- length(contr.names)  
        n.arrays <- ncol(controls.elist)        
        # (IvanG) normally there are 48 blocks, but just in case, let's make it flexible
        # denfensive programming is defensive, yeah
        n.blocks <- length(unique(controls.elist$gene$Block))
        
        # (IvanG) log transform and 
        # merge all control readout data columnwise
        # i.e. array after array
        # As such, length(x)=nrow(E)*ncol(E)
        x <- c(log2(controls.elist$E))                
        
        # (IvanG) Create matrix of explanatory variables for rlm
        # rows: correspond to the rows of x
        dummies<-matrix(0, 
                        ncol=n.arrays+n.blocks+n.controls, 
                        nrow=length(x))
        
        # fill the matrix with 1's in such manner that for each row
        # the combination of 1's and 0's is unique
        current.row <- 1
        for(arrN in 1:n.arrays){
                for(blockN in 1:n.blocks){
                        #names of controls in current block
                        conrols.in.current.block<-
                                controls.elist$genes$Name[controls.elist$genes$Block==blockN]
                        #loop through ALL controls ...
                        for(controlN in 1:n.controls){
                                #but if the current control is NOT in current block, 
                                #move to the next one
                                if (!contr.names[controlN] %in% conrols.in.current.block) {next}
                                        #fill 1's in appropriate positions: array, block, and control
                                        dummies[current.row, arrN] <- 1
                                        dummies[current.row, n.arrays+blockN] <- 1
                                        dummies[current.row, n.arrays+n.blocks+controlN] <- 1
                                        current.row <- current.row + 1
                        }
                }
        }
        
        # add 3 rows under the dummies matrix 
        # (likely to denote which columns pertain to arrays, which ones to blocks, 
        # and which ones to controls)
        # with  -1's under first n.arrays columns, then under n.blocks.contr columns,
        # then under n.controls columns
        dummies <- rbind(dummies, 
                         c(rep(-1, n.arrays), rep(0,n.blocks), rep(0,n.controls)),
                         c(rep(0, n.arrays), rep(-1,n.blocks), rep(0,n.controls)), 
                         c(rep(0, n.arrays), rep(0, n.blocks), rep(-1,n.controls)))
        # add intercept (?)
        dummies <- cbind(rep(1, {length(x)+3}), dummies)
        
        
        rlm.result <- rlm(y=c(x,0,0,0), x=dummies, method="M", psi=psi.bisquare)
        
        # Apply the fit coefficients to all of original values to normalize them
        result <- elist
        result$E <- normalize.elist(elist=elist,
                                    n.arrays=n.arrays,
                                    n.blocks=n.blocks,
                                    rlm.result=rlm.result)  
        
        # Normalize the control values in the same way
        contr.normalized <- normalize.elist(elist=controls.elist,
                                            n.arrays=n.arrays,
                                            n.blocks=n.blocks,
                                            rlm.result=rlm.result) 
        
        # create boxplots for raw and normalized values
        if(!is.null(output.path)){
                output.path<-paste0(output.path
                                    ,"/rlm_plots")
                dir.create(output.path, showWarnings = FALSE)
                
                png(paste0(output.path,"/rlm_fitted_values.png"),
                    width=4000, height=2000, pointsize=10, res=300)
                boxplot(matrix(rlm.result$fitted.values[1:length(x)], 
                               ncol=n.arrays), 
                        col="gray", cex=0.6, cex.axis=0.6)
                dev.off()
                
                png(paste0(output.path,"/rlm_all_raw.png"),
                    width=4000, height=2000, pointsize=10, res=300)
                boxplot(log2(elist$E), col="gray", cex=0.6, cex.axis=0.6)
                dev.off()
                
                png(paste0(output.path,"/rlm_all_norm.png"),
                    width=4000, height=2000, pointsize=10, res=300)
                boxplot(result$E, col="gray", cex=0.6, cex.axis=0.6)
                dev.off()
                
                png(paste0(output.path,"/rlm_control_raw.png"),
                    width=4000, height=2000, pointsize=10, res=300)
                boxplot(log2(controls.elist$E), col="gray", cex=0.6, cex.axis=0.6)
                dev.off()
                
                png(paste0(output.path,"/rlm_control_norm.png"),
                    width=4000, height=2000, pointsize=10, res=300)
                boxplot(contr.normalized, col="gray", cex=0.6, cex.axis=0.6)
                dev.off()
                
        }
        return(result)
}
#++++++++++++++++++++++++ normalizeRLM +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ normalizeArrays +++++++++++++++++++++++++++++++++++
# (IvanG) re-formatted. For the pleasure of the eye and the solace of the soul.
normalizeArrays <- function(elist=NULL, 
                            method="quantile",
                            cyclicloess.method="pairs", 
                            group1=NULL, 
                            group2=NULL, 
                            controls.elist=NULL,
                            gpr.path=NULL, 
                            targets.path=NULL, 
                            contr.names=NULL, 
                            output.path=NULL){
        
        if(is.null(elist)) {
                stop("ERROR: Not all mandatory arguments have been defined!")
        }
        
        
        if(method=="none"){
                elist.normalized<-elist
                elist.normalized$E<-log2(elist$E)
                return(elist.normalized)
        }
        
        if(!is.null(group1) && 
           !is.null(group2)) {
                if(method=="quantile"){
                        tmp1<-normalizeBetweenArrays(elist[,group1],
                                                     method=method, 
                                                     ties=FALSE)
                        tmp2<-normalizeBetweenArrays(elist[,group2], 
                                                     method=method, 
                                                     ties=FALSE)
                        elist.normalized <- cbind(tmp1, tmp2)
                }else if(method=="cyclicloess"){
                        tmp1 <- normalizeBetweenArrays(elist[,group1], 
                                                       method=method,
                                                       cyclicloess.method)  
                        #cyclicloess.method: ("pairs"|"fast"|"affy")
                        tmp2 <- normalizeBetweenArrays(elist[,group2], 
                                                       method=method,
                                                       cyclicloess.method)
                        elist.normalized <- cbind(tmp1, tmp2)
                }else if(method=="vsn"){
                        tmp1 <- normalizeVSN(elist[,group1])
                        tmp2 <- normalizeVSN(elist[,group2])
                        elist.normalized <- cbind(tmp1, tmp2)
                }else if(method=="rlm"){
                        tmp1 <- normalizeRLM(elist=elist[,group1],
                                             controls.elist=controls.elist, 
                                             gpr.path=gpr.path,
                                             targets.path=targets.path,
                                             contr.names=NULL,
                                             output.path=output.path)
                        # (IvanG) was "group1" below, needs to be "group2". Fixed.
                        tmp2 <- normalizeRLM(elist=elist[,group2],
                                             controls.elist=controls.elist, 
                                             gpr.path=gpr.path,
                                             targets.path=targets.path,
                                             contr.names=NULL,
                                             output.path=output.path)
                        elist.normalized <- cbind(tmp1, tmp2)
                }
        }else{
                if(method=="quantile"){
                        elist.normalized <-
                                normalizeBetweenArrays(elist, 
                                                       method=method,
                                                       ties=FALSE)
                }else if(method=="cyclicloess"){
                        elist.normalized <-
                                normalizeBetweenArrays(elist,
                                                       method=method, 
                                                       cyclicloess.method)
                        #cyclicloess.method: ("pairs"|"fast"|"affy")
                }else if(method=="vsn"){
                        elist.normalized <- normalizeVSN(elist)
                }else if(method=="rlm"){
                        elist.normalized <-
                                normalizeRLM(elist=elist, 
                                             controls.elist=controls.elist,
                                             gpr.path=gpr.path, 
                                             targets.path=targets.path, 
                                             contr.names=NULL,
                                             output.path=output.path)
                }
        }
        return(elist.normalized)
}
#+++++++++++++++++++++++ normalizeArrays +++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ selectFeatures.plotHeatmap ++++++++++++++++++++++++++++++++++++
selectFeatures.plotHeatmap <- function(features=NULL,
                                       diffAnalysis.results=NULL,
                                       elist=NULL, 
                                       columns1=NULL, 
                                       columns2=NULL,
                                       output.path=NULL,
                                       filename=NULL,
                                       NM_index=TRUE,
                                       title="Biomarker Heatmap" ){
        require(gtools)
    
        if(is.null(features) || 
           is.null(elist) || 
           is.null(columns1) || 
           is.null(columns2)) {
                stop("ERROR: Not all mandatory arguments have been defined!")
        }
        
        #(IvanG) use columns1 and 2 (column names)
        #to subset the arrays of interest
        datamatrix <- elist$E[,c(columns1,columns2)]
        brc <- paste(elist$genes[,1],
                     elist$genes[,3], 
                     elist$genes[,2], 
                     sep=" ")
        
        rows <- match(features,brc)
        
        datamatrix <- datamatrix[rows,,drop=FALSE]
        
        if(NM_index){
                NM.vector<-elist$genes[rows,"Name"]
                NM.vector<-sapply(NM.vector,
                                  function(entry){
                                          unlist(strsplit(
                                                  unlist(strsplit(entry,split=":"))[2],
                                                  split = "~"))[1]  
                                  })
                        
                rownames(datamatrix) <- paste(brc[rows],
                                              NM.vector,
                                              sep=" ")
        }else{
                rownames(datamatrix) <- brc[rows]
        }
        my.dist <- function(x){
                as.dist(1 - abs(cor(t(x), method="pearson")))
                #methods: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
                
        }
        my.hclust <- function(d){
                hclust(d, method="complete")
                #methods:"ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
        }
        
        if(nrow(datamatrix)==1){
                return(datamatrix)
        }
        
        #setup two colors to quickly denote samples belonging to different groups
        RowSideColors<-c(rep("firebrick",length(columns1)), rep("green4",length(columns2)))
        
        if(!is.null(output.path)){
                if(is.null(filename)){
                        filename<-"biomarker_heatmap.pdf"
                }
                
                pdf(paste(output.path,
                           filename,
                           sep="/"),
                     pointsize=15) 
                
                par(cex.main=1)                                      #set main title font size
                heatmap.2.mod(t(datamatrix),             
                              Rowv=TRUE,                             #row and col clustering
                              Colv = TRUE,
                              #hclustfun = my.hclust,
                              #distfun = my.dist,
                              #dendrogram="both",                    #dendrogram  for both
                              RowSideColors=RowSideColors,
                              col=blueyellow(300),                   #color scheme
                              scale="none",                          #scale by column
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
                              #cexRow=,                              #column labels' size
                              labRow=row.names(t(datamatrix)),       #row labels
                              labCol=colnames(t(datamatrix))         #column labels: limit to 30 chars
                )
                dev.off()
        } else {
                par(cex.main=1)                                      #set main title font size
                heatmap.2.mod(t(datamatrix),             
                              Rowv=TRUE,                             #row and col clustering
                              Colv = TRUE,
                              #hclustfun = my.hclust,
                              #distfun = my.dist,
                              #dendrogram="both",                    #dendrogram  for both
                              RowSideColors=RowSideColors,
                              col=blueyellow(300),                  #color scheme
                              scale="none",                          #scale by column
                              main=title,                            #title of the heatmap
                              margins = c(15,15),                    #margins around heatmap    
                              key=TRUE,                              #legend is present
                              key.title="",                      #no title for legend
                              key.xlab="",                           #no label for legend axis
                              key.ylab="",
                              keysize=1.1,                           #size of legend strip
                              symkey=FALSE,                          #do not make the colors symmetrical around 0
                              density.info="density",                   #no density information
                              trace="none",                          #
                              #cexRow=,                              #column labels' size
                              labRow=row.names(t(datamatrix)),       #row labels
                              labCol=colnames(t(datamatrix))         #column labels: limit to 30 chars
                )
        }
        #declare matrix to return
        to.return<-NULL
        #if diffAnalysis data were supplied, include
        if(!is.null(diffAnalysis.results)){
                #select rows corresponding to necessary features
                diffA_rows<-match(features,diffAnalysis.results$BRC)
                to.return<-diffAnalysis.results[diffA_rows,]
                #check if the rows selected from elist
                #and the rows selected from diff analysis df are identical
                if (!identical(brc[rows],to.return$BRC)) {
                        message("selectFeatures.plotHeatmap: Diffanalysis <-> elist mismatch problem!\n")
                }
        }
        #add annotation,
        #datamatrix,
        to.return<-cbind(to.return,
                         elist$genes[rows,],
                         datamatrix)
        #...and proper row names
        rownames(to.return)<-rownames(datamatrix)
        return(to.return)
}
#++++++++++++++++++++++ selectFeatures.plotHeatmap ++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ pvaluePlot +++++++++++++++++++++++++++++++++++++++
pvaluePlot <- function(elist=NULL, 
                       group1=NULL,
                       group2=NULL, 
                       method="tTest",
                       output.path=NULL, 
                       tag="",
                       above=500,
                       between=200,
                       adjust=FALSE){
    
    if(is.null(elist) || 
       is.null(group1) || 
       is.null(group2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    row.number <- nrow(elist$E)
    col.number <- ncol(elist$E)
    pvalues <- vector(mode="numeric", length=row.number)

    if(method=="mMs"){
            mMs.matrix1<-mMsMatrixRancho(length(group1),length(group2))
            mMs.matrix2<-mMsMatrixRancho(length(group2),length(group1))
    }

    for(zeile in 1:row.number) {
            vector1<-elist$E[zeile, group1]
            vector2<-elist$E[zeile, group2]
            
        fold <- mean(vector1)/mean(vector2)
        if(method=="tTest"){
                p <- try(t.test(x=vector1,
                                y=vector2)$p.value, TRUE)
                if(!is.numeric(p)){
                        p <- 1
                }
        }else if(method=="mMs"){
                min.p1 <- calc.mMs.pval(vector1,vector2,mMs.matrix1,above,between)
                min.p2 <- calc.mMs.pval(vector2,vector1,mMs.matrix2,above,between)
                p <- min(min.p1,min.p2)
        }
        pvalues[zeile] <- p
    }
    
    if(adjust){
        pvalues <- p.adjust(p=pvalues, method="fdr")
    }

    number1 <- length(pvalues[pvalues < 0.05])
    number2 <- length(pvalues[pvalues < 0.01])
    number3 <- length(pvalues[pvalues < 0.05/row.number])
    ylimit <- min(log2(0.04/row.number), 
                  log2(min(pvalues))-1)

    if(is.null(output.path)) {
        if(adjust){
            plot(sort(pvalues), 
                 ylab="fdr", 
                 main=paste("FDRs\n(",
                            number1,
                            " < 0.05, ", 
                            number2, 
                            " < 0.01)", 
                            sep=""))
            abline(h=0.01, 
                   lty=2, 
                   col="red")
            abline(h=0.05, 
                   lty=2, 
                   col="red")
        }else{
            plot(log2(sort(pvalues)),
                 ylim=c(ylimit,0),
                 ylab="log2(p-value)",
                 main=paste("p-values\n(", 
                            number1, 
                            " < 0.05, ",
                            number2,
                            " < 0.01, ",
                            number3,
                            " < 0.05 after Bonferroni)", 
                            sep=""))
            abline(h=log2(0.05), 
                   lty=2,
                   col="red")
            abline(h=log2(0.01), 
                   lty=2,
                   col="green")
            abline(h=log2(0.05/row.number), 
                   lty=2, 
                   col="blue")
            legend("bottomright", 
                   legend=c("p-value = 0.05",
                            "p-value = 0.01",
                            "bonferroni"),
                   col=c("red", "green", "blue"), 
                   lty=c(2,2,2))
        }
    } else {
        if(adjust){
            tiff(paste(output.path,
                       "/pvalues_fdr",
                       tag, 
                       ".tiff",
                       sep=""),
              width=2000, 
              height=2000, 
              pointsize=15, 
              compression="lzw",
              res=300)
                
                plot(sort(pvalues),
                     ylab="fdr", 
                     main=paste("FDRs\n(", 
                                number1,
                                " < 0.05, ",
                                number2, 
                                " < 0.01)", 
                                sep=""))
                abline(h=0.01,
                       lty=2, 
                       col="red")
                abline(h=0.05, 
                       lty=2, 
                       col="red")
            dev.off()
            
        }else{
            tiff(paste(output.path,"/pvalues",
                       tag, 
                       ".tiff",
                       sep=""),
                 width=2000,
                 height=2000,
                 pointsize=15,
                 compression="lzw",
                 res=300)
                
                plot(log2(sort(pvalues)), 
                     ylim=c(ylimit,0),
                     ylab="log2(p-value)",
                     main=paste("p-values\n(", 
                                number1,
                                " < 0.05, ", 
                                number2, 
                                " < 0.01, ", 
                                number3,
                                " < 0.05 after Bonferroni)", 
                                sep=""))
                abline(h=log2(0.05),
                       lty=2, 
                       col="red")
                abline(h=log2(0.01), 
                       lty=2,
                       col="green")
                abline(h=log2(0.05/row.number),
                       lty=2, 
                       col="blue")
                legend("bottomright", 
                       legend=c("p-value = 0.05",
                                "p-value = 0.01",
                                "bonferroni"), 
                       col=c("red", 
                             "green",
                             "blue"), 
                       lty=c(2,2,2))
            dev.off()
        }
    }
}
#+++++++++++++++++++++++++++ pvaluePlot +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ volcanoPlot +++++++++++++++++++++++++++++++++++++++
volcanoPlot <- function(elist=NULL, 
                        group1=NULL, 
                        group2=NULL, 
                        method="tTest",
                        p.thresh=NULL, 
                        fold.thresh=NULL, 
                        output.path=NULL, 
                        tag="", 
                        above=500, 
                        between=200,
                        title="Volcano Plot"){
    
        if(is.null(elist) || 
           is.null(group1) || 
           is.null(group2)) {
                stop("ERROR: Not all mandatory arguments have been defined!")
        }
        row.number <- nrow(elist$E)
        col.number <- ncol(elist$E)
        mark <- c()
        output <- matrix(nrow=row.number, ncol={4+ncol(elist$E)})
        colnames(output) <- c("BRC", 
                              "P-value", 
                              "Fold change", 
                              "Description",
                              colnames(elist$E))
        
        if(method=="mMs"){
                mMs.matrix1<-mMsMatrixRancho(length(group1),length(group2))
                mMs.matrix2<-mMsMatrixRancho(length(group2),length(group1))
        }
        for(zeile in 1:row.number) {
                vector1<-elist$E[zeile, group1]
                vector2<-elist$E[zeile, group2]
                
                fold <- mean(vector1)/mean(vector2)
                
                if(method=="tTest"){
                        p <- try(t.test(x=vector1,
                                        y=vector2)$p.value, TRUE)
                        if(!is.numeric(p)){
                                p <- 1
                        }
                }else if(method=="mMs"){
                        
                        min.p1 <- calc.mMs.pval(vector1,vector2,mMs.matrix1,above,between)
                        min.p2 <- calc.mMs.pval(vector2,vector1,mMs.matrix2,above,between)
                        
                        p <- min(min.p1,min.p2)
                }
                if(!is.null(p.thresh) && 
                   !is.null(fold.thresh)){   # --> A & B
                        if(p < p.thresh && 
                           (fold > fold.thresh || 
                           fold < 1/fold.thresh)) {
                                mark <- c(mark,zeile)
                        }
                }else if(!is.null(p.thresh) && 
                         is.null(fold.thresh)) {    # --> A & -B
                        if(p < p.thresh) {
                                mark <- c(mark,zeile)
                        }
                }else if(is.null(p.thresh) && 
                         !is.null(fold.thresh)) {    # --> -A & B
                        if(fold > fold.thresh || 
                           fold < 1/fold.thresh) {
                                mark <- c(mark,zeile)
                        }
                }
                output[zeile,] <- c(paste(elist$genes[zeile,1],
                                          elist$genes[zeile,3], 
                                          elist$genes[zeile,2], 
                                          sep=" "), 
                                    p, 
                                    fold,
                                    elist$genes[zeile,4], 
                                    elist$E[zeile,])
        }
        #write.table(x=output, file=paste(output.path, "/volcano.txt", sep=""),
        #  sep="\t", eol="\n", row.names=FALSE, quote=FALSE)
        if(is.null(output.path) && 
           is.null(p.thresh) && 
           is.null(fold.thresh)){
                # --> -A & -B & -C
                plot(x=log2(as.numeric(output[,3])), 
                     y=-log10(as.numeric(output[,2])),
                     xlab="log2(fold change)", 
                     ylab="-log10(p-value)", 
                     main=title)
                
        }else if(!is.null(output.path) && 
                 is.null(p.thresh) &&
                 is.null(fold.thresh)) {    # --> A & -B & -C
                
                tiff(paste(output.path,"/volcano", 
                           tag, 
                           ".tiff",
                           sep=""), 
                     width=2000,
                     height=2000, 
                     pointsize=15, 
                     compression="lzw", 
                     res=300)
                
                plot(x=log2(as.numeric(output[,3])),
                     y=-log10(as.numeric(output[,2])),
                     xlab="log2(fold change)",
                     ylab="-log10(p-value)",
                     main=title)
                dev.off()
                
        }else if(is.null(output.path) && 
                 !is.null(p.thresh) &&
                 is.null(fold.thresh)) { # --> -A & B & -C
                
                plot(x=log2(as.numeric(output[,3])),
                     y=-log10(as.numeric(output[,2])),
                     xlab="log2(fold change)", 
                     ylab="-log10(p-value)",
                     main=title)
                
                points(x=log2(as.numeric(output[mark,3])),
                       y=-log10(as.numeric(output[mark,2])),
                       col="red")
                abline(h=-log10(p.thresh),
                       lty=2, col="red")
                
        }else if(!is.null(output.path) && 
                 !is.null(p.thresh) &&
                 is.null(fold.thresh)) {   # --> A & B & -C
                
                tiff(paste(output.path,"/volcano",
                           tag, 
                           ".tiff",
                           sep=""),
                     width=2000,
                     height=2000, 
                     pointsize=15, 
                     compression="lzw",
                     res=300)
                
                plot(x=log2(as.numeric(output[,3])), 
                     y=-log10(as.numeric(output[,2])),
                     xlab="log2(fold change)",
                     ylab="-log10(p-value)",
                     main=title)
                points(x=log2(as.numeric(output[mark,3])),
                       y=-log10(as.numeric(output[mark,2])), 
                       col="red")
                abline(h=-log10(p.thresh), 
                       lty=2, 
                       col="red")
                dev.off()  
                
        }else if(is.null(output.path) &&
                 is.null(p.thresh) &&
                 !is.null(fold.thresh)) {   # --> -A & -B & C
                
                plot(x=log2(as.numeric(output[,3])),
                     y=-log10(as.numeric(output[,2])),
                     xlab="log2(fold change)", 
                     ylab="-log10(p-value)", 
                     main=title)
                
                points(x=log2(as.numeric(output[mark,3])),
                       y=-log10(as.numeric(output[mark,2])),
                       col="red")
                
                abline(v=log2(fold.thresh), 
                       lty=2, 
                       col="red")
                abline(v=log2(1/fold.thresh), 
                       lty=2, 
                       col="red")
                
        }else if(!is.null(output.path) && 
                 is.null(p.thresh) &&
                 !is.null(fold.thresh)) {    # --> A & -B & C
                
                tiff(paste(output.path,"/volcano", 
                           tag, 
                           ".tiff",
                           sep=""),
                     width=2000,
                     height=2000, 
                     pointsize=15,
                     compression="lzw",
                     res=300)
                
                plot(x=log2(as.numeric(output[,3])),
                     y=-log10(as.numeric(output[,2])),
                     xlab="log2(fold change)", 
                     ylab="-log10(p-value)", 
                     main=title)
                
                points(x=log2(as.numeric(output[mark,3])),
                       y=-log10(as.numeric(output[mark,2])), 
                       col="red")
                
                abline(v=log2(fold.thresh),
                       lty=2, 
                       col="red")
                
                abline(v=log2(1/fold.thresh), 
                       lty=2, 
                       col="red")
                dev.off()
                
        }else if(is.null(output.path) &&
                 !is.null(p.thresh) &&
                 !is.null(fold.thresh)) {    # --> -A & B & C
                
                plot(x=log2(as.numeric(output[,3])),
                     y=-log10(as.numeric(output[,2])),
                     xlab="log2(fold change)",
                     ylab="-log10(p-value)", 
                     main=title)
                
                points(x=log2(as.numeric(output[mark,3])),
                       y=-log10(as.numeric(output[mark,2])),
                       col="red")
                abline(h=-log10(p.thresh),
                       lty=2,
                       col="red")
                abline(v=log2(fold.thresh), 
                       lty=2,
                       col="red")
                abline(v=log2(1/fold.thresh), 
                       lty=2, 
                       col="red")
                
        }else if(!is.null(output.path) && 
                 !is.null(p.thresh) &&
                 !is.null(fold.thresh)) {    # --> A & B & C
                
                tiff(paste(output.path,"/volcano", 
                           tag, 
                           ".tiff",
                           sep=""), 
                     width=2000,
                     height=2000, 
                     pointsize=15, 
                     compression="lzw", 
                     res=300)
                
                plot(x=log2(as.numeric(output[,3])),
                     y=-log10(as.numeric(output[,2])),
                     xlab="log2(fold change)",
                     ylab="-log10(p-value)",
                     main=title)
                
                points(x=log2(as.numeric(output[mark,3])),
                       y=-log10(as.numeric(output[mark,2])),
                       col="red")
                
                abline(h=-log10(p.thresh),
                       lty=2, 
                       col="red")
                
                abline(v=log2(fold.thresh),
                       lty=2, 
                       col="red")
                
                abline(v=log2(1/fold.thresh),
                       lty=2,
                       col="red")
                
                dev.off()
        }
}
#+++++++++++++++++++++++++++ volcanoPlot +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ diffAnalysis ++++++++++++++++++++++++++++++++++++++
diffAnalysis <- function(input=NULL, 
                         label1=NULL, 
                         label2=NULL, 
                         class1=NULL,
                         class2=NULL,
                         sortBY="fold_change",
                         features=NULL, 
                         feature.names=NULL,
                         output.path=NULL, 
                         filename="diffAnalysis.txt",
                         above=500,
                         between=200){
        
    cat("Running Rancho-modified diffAnalysis...\n")
        
    if(is.null(input) || 
       is.null(label1) || 
       is.null(label2) || 
       is.null(class1) ||
       is.null(class2)) {
        stop("ERROR: Not all mandatory arguments have been defined!")
    }
    row.number <- nrow(input)
    col.number <- ncol(input)
    output <- matrix(nrow=row.number, ncol=12)
    colnames(output) <- c("BRC", 
                          "tTest", 
                          "FDR_t", 
                          "mMs",
                          "FDR_mMs", 
                          "fold_change", 
                          paste("mean", class1, sep="_"),
                          paste("mean", class2, sep="_"), 
                          paste("median", class1, sep="_"),
                          paste("median", class2, sep="_"), 
                          paste("sd", class1, sep="_"),
                          paste("sd", class2, sep="_"))
    
    #calculate the minimal M statistic p-value matrices
    #for comparisons between groups 1 vs. 2 and 2 vs. 1
    mMs.matrix1<-mMsMatrixRancho(length(label1),length(label2))
    mMs.matrix2<-mMsMatrixRancho(length(label2),length(label1))
    
    for(zeile in 1:row.number) {
        #row names
        output[zeile,1] <- rownames(input)[zeile]
        
        #define vector1 and vector2 
        vector1<-as.numeric(input[zeile,label1])
        vector2<-as.numeric(input[zeile,label2])

        p <- try(t.test(x=vector1,
                        y=vector2)$p.value, 
                 TRUE)
        if(!is.numeric(p)){
                output[zeile,2] <- 1
        }else{
                output[zeile,2] <- p
        }
        
        min.p1 <- calc.mMs.pval(vector1,vector2,mMs.matrix1,above,between)
        min.p2 <- calc.mMs.pval(vector2,vector1,mMs.matrix2,above,between)
        mean1<-mean(vector1)
        mean2<-mean(vector2)
        
        output[zeile,4] <- min(min.p1,min.p2) 
        
        output[zeile,6] <- mean1/mean2
        output[zeile,7] <- mean1
        output[zeile,8] <- mean2
        output[zeile,9] <- median(vector1)
        output[zeile,10] <- median(vector2)
        output[zeile,11] <- sd(vector1)       
        output[zeile,12] <- sd(vector2)
    }
    output[,3] <- p.adjust(p=as.numeric(output[,2]), method="fdr")
    output[,5] <- p.adjust(p=as.numeric(output[,4]), method="fdr")
    
    #(IvanG) added stringsAsFactors
    outputDF<-data.frame(output,
                         stringsAsFactors = FALSE)
    #(IvanG) convert DF columns to numeric
    for (colNum in 2:ncol(output)){
            outputDF[,colNum]<-as.numeric(outputDF[,colNum])
    }
    #sort by specified column
    if (!is.null(sortBY) && 
        (sortBY %in% colnames(outputDF))){
            if (sortBY=="fold_change"){
                    outputDF<-outputDF[order(outputDF[,sortBY],decreasing = TRUE),]
            } else {
                    outputDF<-outputDF[order(outputDF[,sortBY],decreasing = FALSE),]
            }
    }
    
    if(!is.null(output.path)){

            if(is.null(features)){
                    #export full DF sorted by specified column
                    write.table(x=outputDF, file=paste(output.path, 
                                                     filename,
                                                     sep="/"), 
                                sep="\t", eol="\n", row.names=FALSE, quote=FALSE)
            }else{
                    #if specific features have been requested, export them unsorted
                    if(is.null(feature.names)){
                            write.table(x=output[features,], file=paste(output.path,
                                                                        filename,
                                                                        sep="/"), 
                                        sep="\t", eol="\n",
                                        row.names=FALSE, quote=FALSE)
                    }else{
                            write.table(x=cbind(output[features,],feature.names),
                                        file=paste(output.path, 
                                                   filename, 
                                                   sep="/"),
                                        sep="\t", eol="\n", row.names=FALSE, quote=FALSE)
                    }
            }
    }
  
    return(outputDF)
}
#++++++++++++++++++++++++++ diffAnalysis +++++++++++++++++++++++++++++++++++++++

#####################################################################################################
#################################### added by Rancho ################################################
#####################################################################################################

#+++++++++++++++++++++++++++ diffAnalysis.tTest ++++++++++++++++++++++++++++++++++++++
diffAnalysis.tTest <- function(input=NULL, 
                               label1=NULL, 
                               label2=NULL, 
                               class1=NULL,
                               class2=NULL, 
                               sortBY="fold_change",
                               features=NULL, 
                               feature.names=NULL,
                               output.path=NULL, 
                               filename="diffAnalysis.txt"){
        
        cat("Running Rancho-modified diffAnalysis.tTest...\n")
        
        if(is.null(input) || 
           is.null(label1) || 
           is.null(label2) || 
           is.null(class1) ||
           is.null(class2)) {
                stop("ERROR: Not all mandatory arguments have been defined!")
        }
        row.number <- nrow(input)
        col.number <- ncol(input)
        colNames <- c("BRC", 
                      "tTest", 
                      "FDR_t",
                      "fold_change", 
                      paste("mean", class1, sep="_"),
                      paste("mean", class2, sep="_"), 
                      paste("median", class1, sep="_"),
                      paste("median", class2, sep="_"), 
                      paste("sd", class1, sep="_"),
                      paste("sd", class2, sep="_"))
        
        output <- matrix(nrow=row.number, ncol=length(colNames))
        colnames(output)<-colNames
        
        for(zeile in 1:row.number) {
                output[zeile,1] <- rownames(input)[zeile]
                
                p <- try(t.test(x=as.numeric(input[zeile,label1]),
                                y=as.numeric(input[zeile,label2]))$p.value, TRUE)
                if(!is.numeric(p)){output[zeile,2] <- 1}else{output[zeile,2] <- p}
                
                #calculate means and medians for each row, 
                #and fill up the table
                mean1<-mean(as.numeric(input[zeile,label1]))
                mean2<-mean(as.numeric(input[zeile,label2]))
                
                output[zeile,4] <- mean1/mean2
                output[zeile,5] <- mean1
                output[zeile,6] <- mean2
                output[zeile,7] <- median(as.numeric(input[zeile, label1]))
                output[zeile,8] <- median(as.numeric(input[zeile, label2]))
                output[zeile,9] <- sd(as.numeric(input[zeile, label1]))       
                output[zeile,10] <- sd(as.numeric(input[zeile, label2]))
        }
        output[,3] <- p.adjust(p=as.numeric(output[,2]), method="fdr")
        
        #(IvanG) added stringsAsFactors
        outputDF<-data.frame(output,
                             stringsAsFactors = FALSE)
        #(IvanG) convert DF columns to numeric
        for (colNum in 2:ncol(output)){
                outputDF[,colNum]<-as.numeric(outputDF[,colNum])
        }
        #sort by specified column
        if (!is.null(sortBY) && 
            (sortBY %in% colnames(outputDF))){
                if (sortBY=="fold_change"){
                        outputDF<-outputDF[order(outputDF[,sortBY],decreasing = TRUE),]
                } else {
                        outputDF<-outputDF[order(outputDF[,sortBY],decreasing = FALSE),]
                }
        }
        
        if(!is.null(output.path)){
                
                if(is.null(features)){
                        #export full DF sorted by specified column
                        write.table(x=outputDF, 
                                    file=paste(output.path, 
                                               filename,
                                               sep="/"), 
                                    sep="\t", 
                                    eol="\n", 
                                    row.names=FALSE, 
                                    quote=FALSE)
                }else{
                        #if specific features have been requested, export them unsorted
                        if(is.null(feature.names)){
                                write.table(x=output[features,], 
                                            file=paste(output.path,
                                                       filename,
                                                       sep="/"), 
                                            sep="\t", 
                                            eol="\n",
                                            row.names=FALSE, 
                                            quote=FALSE)
                        }else{
                                write.table(x=cbind(output[features,],feature.names),
                                            file=paste(output.path, 
                                                       filename, 
                                                       sep="/"),
                                            sep="\t", 
                                            eol="\n", 
                                            row.names=FALSE, 
                                            quote=FALSE)
                        }
                }
        }
        return(outputDF)
}
#++++++++++++++++++++++++++ diffAnalysis.tTest +++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ yellowblue ++++++++++++++++++++++++++++++++++++
yellowblue<-function(n){
        colorpanel(n, "yellow", "black", "blue")
}
#++++++++++++++++++++++ yellowblue ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ blueyellow ++++++++++++++++++++++++++++++++++++
blueyellow<-function(n){
        colorpanel(n, "blue","black", "yellow")
}
#++++++++++++++++++++++ blueyellow ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ heatmap.2.mod ++++++++++++++++++++++++++++++++++++

#modified heatmap.2 function
#2015-03-13
#replaced cexRow default formula with 70/nr*(1/(2*log10(nr)))
#2015-03-30
#replaced cexRow default formula with ifelse(nr<15,2,70/nr*(1/(2*log10(nr))))
#to accomodate cases when very few genes are returned
#2015-06-29
#same done for columns

heatmap.2.mod<-function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                         distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                            "row", "column", "none"), reorderfun = function(d, w) reorder(d, 
                                                                                                                                          w), symm = FALSE, scale = c("none", "row", "column"), 
                         na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
                         symbreaks = any(x < 0, na.rm = TRUE) || scale != "none", 
                         col = "heat.colors", colsep, rowsep, sepcolor = "white", 
                         sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", 
                         na.color = par("bg"), trace = c("column", "row", "both", 
                                                         "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), 
                         linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors, 
                         cexRow = ifelse(nr<25,0.7,70/nr*(1/(2*log10(nr)))), cexCol = ifelse(nc<25,0.9,60/nc*(1/(2*log10(nc)))), labRow = NULL, 
                         labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, 
                                                                                 NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5, 
                         key = TRUE, keysize = 1.5, density.info = c("histogram", 
                                                                     "density", "none"), denscol = tracecol, symkey = any(x < 
                                                                                                                                  0, na.rm = TRUE) || symbreaks, densadj = 0.25, key.title = NULL, 
                         key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, key.ytickfun = NULL, 
                         key.par = list(), main = NULL, xlab = NULL, ylab = NULL, 
                         lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, ...) 
{
        scale01 <- function(x, low = min(x), high = max(x)) {
                x <- (x - low)/(high - low)
                x
        }
        retval <- list()
        scale <- if (symm && missing(scale)) 
                "none"
        else match.arg(scale)
        dendrogram <- match.arg(dendrogram)
        trace <- match.arg(trace)
        density.info <- match.arg(density.info)
        if (length(col) == 1 && is.character(col)) 
                col <- get(col, mode = "function")
        if (!missing(breaks) && (scale != "none")) 
                warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
                        "specified can produce unpredictable results.", "Please consider using only one or the other.")
        if (is.null(Rowv) || is.na(Rowv)) 
                Rowv <- FALSE
        if (is.null(Colv) || is.na(Colv)) 
                Colv <- FALSE
        else if (all(Colv == "Rowv")) 
                Colv <- Rowv
        if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
                stop("`x' must be a numeric matrix")
        nr <- di[1]
        nc <- di[2]
        if (nr <= 1 || nc <= 1) 
                stop("`x' must have at least 2 rows and 2 columns")
        if (!is.numeric(margins) || length(margins) != 2) 
                stop("`margins' must be a numeric vector of length 2")
        if (missing(cellnote)) 
                cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
        if (!inherits(Rowv, "dendrogram")) {
                if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) && 
                    (dendrogram %in% c("both", "row"))) {
                        if (is.logical(Colv) && (Colv)) 
                                dendrogram <- "column"
                        else dendrogram <- "none"
                        warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                                dendrogram, "'. Omitting row dendogram.")
                }
        }
        if (!inherits(Colv, "dendrogram")) {
                if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) && 
                    (dendrogram %in% c("both", "column"))) {
                        if (is.logical(Rowv) && (Rowv)) 
                                dendrogram <- "row"
                        else dendrogram <- "none"
                        warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                                dendrogram, "'. Omitting column dendogram.")
                }
        }
        if (inherits(Rowv, "dendrogram")) {
                ddr <- Rowv
                rowInd <- order.dendrogram(ddr)
                if (length(rowInd) > nr || any(rowInd < 1 | rowInd > 
                                               nr)) 
                        stop("Rowv dendrogram doesn't match size of x")
        }
        else if (is.integer(Rowv)) {
                browser()
                hcr <- hclustfun(distfun(x))
                ddr <- as.dendrogram(hcr)
                ddr <- reorderfun(ddr, Rowv)
                rowInd <- order.dendrogram(ddr)
                if (nr != length(rowInd)) 
                        stop("row dendrogram ordering gave index of wrong length")
        }
        else if (isTRUE(Rowv)) {
                Rowv <- rowMeans(x, na.rm = na.rm)
                hcr <- hclustfun(distfun(x))
                ddr <- as.dendrogram(hcr)
                ddr <- reorderfun(ddr, Rowv)
                rowInd <- order.dendrogram(ddr)
                if (nr != length(rowInd)) 
                        stop("row dendrogram ordering gave index of wrong length")
        }
        else {
                rowInd <- nr:1
        }
        if (inherits(Colv, "dendrogram")) {
                ddc <- Colv
                colInd <- order.dendrogram(ddc)
                if (length(colInd) > nc || any(colInd < 1 | colInd > 
                                               nc)) 
                        stop("Colv dendrogram doesn't match size of x")
        }
        else if (identical(Colv, "Rowv")) {
                if (nr != nc) 
                        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
                if (exists("ddr")) {
                        ddc <- ddr
                        colInd <- order.dendrogram(ddc)
                }
                else colInd <- rowInd
        }
        else if (is.integer(Colv)) {
                hcc <- hclustfun(distfun(if (symm) 
                        x
                        else t(x)))
                ddc <- as.dendrogram(hcc)
                ddc <- reorderfun(ddc, Colv)
                colInd <- order.dendrogram(ddc)
                if (nc != length(colInd)) 
                        stop("column dendrogram ordering gave index of wrong length")
        }
        else if (isTRUE(Colv)) {
                Colv <- colMeans(x, na.rm = na.rm)
                hcc <- hclustfun(distfun(if (symm) 
                        x
                        else t(x)))
                ddc <- as.dendrogram(hcc)
                ddc <- reorderfun(ddc, Colv)
                colInd <- order.dendrogram(ddc)
                if (nc != length(colInd)) 
                        stop("column dendrogram ordering gave index of wrong length")
        }
        else {
                colInd <- 1:nc
        }
        retval$rowInd <- rowInd
        retval$colInd <- colInd
        retval$call <- match.call()
        x <- x[rowInd, colInd]
        x.unscaled <- x
        cellnote <- cellnote[rowInd, colInd]
        if (is.null(labRow)) 
                labRow <- if (is.null(rownames(x))) 
                        (1:nr)[rowInd]
        else rownames(x)
        else labRow <- labRow[rowInd]
        if (is.null(labCol)) 
                labCol <- if (is.null(colnames(x))) 
                        (1:nc)[colInd]
        else colnames(x)
        else labCol <- labCol[colInd]
        if (scale == "row") {
                retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
                x <- sweep(x, 1, rm)
                retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
                x <- sweep(x, 1, sx, "/")
        }
        else if (scale == "column") {
                retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
                x <- sweep(x, 2, rm)
                retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
                x <- sweep(x, 2, sx, "/")
        }
        if (missing(breaks) || is.null(breaks) || length(breaks) < 
            1) {
                if (missing(col) || is.function(col)) 
                        breaks <- 16
                else breaks <- length(col) + 1
        }
        if (length(breaks) == 1) {
                if (!symbreaks) 
                        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                                      length = breaks)
                else {
                        extreme <- max(abs(x), na.rm = TRUE)
                        breaks <- seq(-extreme, extreme, length = breaks)
                }
        }
        nbr <- length(breaks)
        ncol <- length(breaks) - 1
        if (class(col) == "function") 
                col <- col(ncol)
        min.breaks <- min(breaks)
        max.breaks <- max(breaks)
        x[x < min.breaks] <- min.breaks
        x[x > max.breaks] <- max.breaks
        if (missing(lhei) || is.null(lhei)) 
                lhei <- c(keysize, 4)
        if (missing(lwid) || is.null(lwid)) 
                lwid <- c(keysize, 4)
        if (missing(lmat) || is.null(lmat)) {
                lmat <- rbind(4:3, 2:1)
                if (!missing(ColSideColors)) {
                        if (!is.character(ColSideColors) || length(ColSideColors) != 
                            nc) 
                                stop("'ColSideColors' must be a character vector of length ncol(x)")
                        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                                              1)
                        lhei <- c(lhei[1], 0.2, lhei[2])
                }
                if (!missing(RowSideColors)) {
                        if (!is.character(RowSideColors) || length(RowSideColors) != 
                            nr) 
                                stop("'RowSideColors' must be a character vector of length nrow(x)")
                        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                                                                   1), 1), lmat[, 2] + 1)
                        lwid <- c(lwid[1], 0.2, lwid[2])
                }
                lmat[is.na(lmat)] <- 0
        }
        if (length(lhei) != nrow(lmat)) 
                stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
        if (length(lwid) != ncol(lmat)) 
                stop("lwid must have length = ncol(lmat) =", ncol(lmat))
        op <- par(no.readonly = TRUE)
        on.exit(par(op))
        layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
        if (!missing(RowSideColors)) {
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        }
        if (!missing(ColSideColors)) {
                par(mar = c(0.5, 0, 0, margins[2]))
                image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        }
        par(mar = c(margins[1], 0, 0, margins[2]))
        x <- t(x)
        cellnote <- t(cellnote)
        if (revC) {
                iy <- nr:1
                if (exists("ddr")) 
                        ddr <- rev(ddr)
                x <- x[, iy]
                cellnote <- cellnote[, iy]
        }
        else iy <- 1:nr
        image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
                      c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
              breaks = breaks, ...)
        retval$carpet <- x
        if (exists("ddr")) 
                retval$rowDendrogram <- ddr
        if (exists("ddc")) 
                retval$colDendrogram <- ddc
        retval$breaks <- breaks
        retval$col <- col
        if (!invalid(na.color) & any(is.na(x))) {
                mmat <- ifelse(is.na(x), 1, NA)
                image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
                      col = na.color, add = TRUE)
        }
        if (is.null(srtCol)) 
                axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
                             offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
                     padj = adjCol[2])
        else {
                if (is.numeric(srtCol)) {
                        if (missing(adjCol) || is.null(adjCol)) 
                                adjCol = c(1, NA)
                        xpd.orig <- par("xpd")
                        par(xpd = NA)
                        xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                                     tick = 0)
                        text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
                                     strheight("M"), labels = labCol, adj = adjCol, 
                             cex = cexCol, srt = srtCol)
                        par(xpd = xpd.orig)
                }
                else warning("Invalid value for srtCol ignored.")
        }
        if (is.null(srtRow)) {
                axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
                     tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
        }
        else {
                if (is.numeric(srtRow)) {
                        xpd.orig <- par("xpd")
                        par(xpd = NA)
                        ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                                     line = -0.5, tick = 0)
                        text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
                             y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
                             srt = srtRow)
                        par(xpd = xpd.orig)
                }
                else warning("Invalid value for srtRow ignored.")
        }
        if (!is.null(xlab)) 
                mtext(xlab, side = 1, line = margins[1] - 1.25)
        if (!is.null(ylab)) 
                mtext(ylab, side = 4, line = margins[2] - 1.25)
        if (!missing(add.expr)) 
                eval(substitute(add.expr))
        if (!missing(colsep)) 
                for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
                                          xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                                                  1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
        if (!missing(rowsep)) 
                for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                                                        1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                                                                                               1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                                          col = sepcolor, border = sepcolor)
        min.scale <- min(breaks)
        max.scale <- max(breaks)
        x.scaled <- scale01(t(x), min.scale, max.scale)
        if (trace %in% c("both", "column")) {
                retval$vline <- vline
                vline.vals <- scale01(vline, min.scale, max.scale)
                for (i in colInd) {
                        if (!is.null(vline)) {
                                abline(v = i - 0.5 + vline.vals, col = linecol, 
                                       lty = 2)
                        }
                        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
                        xv <- c(xv[1], xv)
                        yv <- 1:length(xv) - 0.5
                        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
                }
        }
        if (trace %in% c("both", "row")) {
                retval$hline <- hline
                hline.vals <- scale01(hline, min.scale, max.scale)
                for (i in rowInd) {
                        if (!is.null(hline)) {
                                abline(h = i - 0.5 + hline.vals, col = linecol, 
                                       lty = 2)
                        }
                        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
                        yv <- rev(c(yv[1], yv))
                        xv <- length(yv):1 - 0.5
                        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
                }
        }
        if (!missing(cellnote)) 
                text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
                     col = notecol, cex = notecex)
        par(mar = c(margins[1], 0, 0, 0))
        if (dendrogram %in% c("both", "row")) {
                flag <- try(gplots:::plot.dendrogram(ddr, horiz = TRUE, axes = FALSE, 
                                                     yaxs = "i", leaflab = "none"))
                if ("try-error" %in% class(flag)) {
                        cond <- attr(flag, "condition")
                        if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
                                stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
                }
        }
        else plot.new()
        par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
        if (dendrogram %in% c("both", "column")) {
                flag <- try(gplots:::plot.dendrogram(ddc, axes = FALSE, xaxs = "i", 
                                                     leaflab = "none"))
                if ("try-error" %in% class(flag)) {
                        cond <- attr(flag, "condition")
                        if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
                                stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
                }
        }
        else plot.new()
        if (!is.null(main)) 
                title(main, cex.main = 1.5 * op[["cex.main"]])
        if (key) {
                mar <- c(5, 4, 2, 1)
                if (!is.null(key.xlab) && is.na(key.xlab)) 
                        mar[1] <- 2
                if (!is.null(key.ylab) && is.na(key.ylab)) 
                        mar[2] <- 2
                if (!is.null(key.title) && is.na(key.title)) 
                        mar[3] <- 1
                par(mar = mar, cex = 0.5, mgp = c(2, 1, 0))
                if (length(key.par) > 0) 
                        do.call(par, key.par)
                tmpbreaks <- breaks
                if (symkey) {
                        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
                        min.raw <- -max.raw
                        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
                        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
                }
                else {
                        min.raw <- min(x, na.rm = TRUE)
                        max.raw <- max(x, na.rm = TRUE)
                }
                z <- seq(min.raw, max.raw, length = length(col))
                image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
                      xaxt = "n", yaxt = "n")
                par(usr = c(0, 1, 0, 1))
                if (is.null(key.xtickfun)) {
                        lv <- pretty(breaks)
                        xv <- scale01(as.numeric(lv), min.raw, max.raw)
                        xargs <- list(at = xv, labels = lv)
                }
                else {
                        xargs <- key.xtickfun()
                }
                xargs$side <- 1
                do.call(axis, xargs)
                if (is.null(key.xlab)) {
                        if (scale == "row") 
                                key.xlab <- "Row Z-Score"
                        else if (scale == "column") 
                                key.xlab <- "Column Z-Score"
                        else key.xlab <- "Value"
                }
                if (!is.na(key.xlab)) {
                        mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5)
                }
                if (density.info == "density") {
                        dens <- density(x, adjust = densadj, na.rm = TRUE)
                        omit <- dens$x < min(breaks) | dens$x > max(breaks)
                        dens$x <- dens$x[-omit]
                        dens$y <- dens$y[-omit]
                        dens$x <- scale01(dens$x, min.raw, max.raw)
                        lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                              lwd = 1)
                        if (is.null(key.ytickfun)) {
                                yargs <- list(at = pretty(dens$y)/max(dens$y) * 
                                                      0.95, labels = pretty(dens$y))
                        }
                        else {
                                yargs <- key.ytickfun()
                        }
                        yargs$side <- 2
                        do.call(axis, yargs)
                        if (is.null(key.title)) 
                                key.title <- "Color Key\nand Density Plot"
                        if (!is.na(key.title)) 
                                title(key.title)
                        par(cex = 0.5)
                        if (is.null(key.ylab)) 
                                key.ylab <- "Density"
                        if (!is.na(key.ylab)) 
                                mtext(side = 2, key.ylab, line = par("mgp")[1], 
                                      padj = 0.5)
                }
                else if (density.info == "histogram") {
                        h <- hist(x, plot = FALSE, breaks = breaks)
                        hx <- scale01(breaks, min.raw, max.raw)
                        hy <- c(h$counts, h$counts[length(h$counts)])
                        lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                              col = denscol)
                        if (is.null(key.ytickfun)) {
                                yargs <- list(at = pretty(hy)/max(hy) * 0.95, 
                                              labels = pretty(hy))
                        }
                        else {
                                yargs <- key.ytickfun()
                        }
                        yargs$side <- 2
                        do.call(axis, yargs)
                        if (is.null(key.title)) 
                                key.title <- "Color Key\nand Histogram"
                        if (!is.na(key.title)) 
                                title(key.title)
                        par(cex = 0.5)
                        if (is.null(key.ylab)) 
                                key.ylab <- "Count"
                        if (!is.na(key.ylab)) 
                                mtext(side = 2, key.ylab, line = par("mgp")[1], 
                                      padj = 0.5)
                }
                else if (is.null(key.title)) 
                        title("Color Key")
                if (trace %in% c("both", "column")) {
                        vline.vals <- scale01(vline, min.raw, max.raw)
                        if (!is.null(vline)) {
                                abline(v = vline.vals, col = linecol, lty = 2)
                        }
                }
                if (trace %in% c("both", "row")) {
                        hline.vals <- scale01(hline, min.raw, max.raw)
                        if (!is.null(hline)) {
                                abline(v = hline.vals, col = linecol, lty = 2)
                        }
                }
        }
        else plot.new()
        retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                        high = retval$breaks[-1], color = retval$col)
        if (!is.null(extrafun)) 
                extrafun()
        invisible(retval)
}
#++++++++++++++++++++++ heatmap.2.mod ++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ ComBat +++++++++++++++++++++++++++++++++
ComBat<-
        function (dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE) 
        {
                cat("Running Rancho-modified ComBat...\n")
                
                if (length(dim(batch)) > 1) {
                        stop("This version of ComBat only allows one batch variable")
                }
                batch <- as.factor(batch)
                batchmod <- model.matrix(~-1 + batch)
                cat("Found", nlevels(batch), "batches.\n")
                n.batch <- nlevels(batch)
                batches <- list()
                for (i in 1:n.batch) {
                        batches[[i]] <- which(batch == levels(batch)[i])
                }
                n.batches <- sapply(batches, length)
                n.array <- sum(n.batches)
                design <- cbind(batchmod, mod)
                check <- apply(design, 2, function(x) all(x == 1))
                design <- as.matrix(design[, !check])
                cat("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)\n")
                if (qr(design)$rank < ncol(design)) {
                        if (ncol(design) == (n.batch + 1)) {
                                stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
                        }
                        if (ncol(design) > (n.batch + 1)) {
                                if ((qr(design[, -c(1:n.batch)])$rank < ncol(design[, 
                                                                                    -c(1:n.batch)]))) {
                                        stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
                                }
                                else {
                                        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
                                }
                        }
                }
                NAs = any(is.na(dat))
                if (NAs) {
                        cat(c("Found", sum(is.na(dat)), "Missing Data Values\n"), 
                            sep = " ")
                }
                cat("Standardizing Data across genes\n")
                if (!NAs) {
                        B.hat <- solve(t(design) %*% design) %*% t(design) %*% 
                                t(as.matrix(dat))
                } else {
                        B.hat = apply(dat, 1, Beta.NA, design)
                }
                
                grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]
                
                if (!NAs) {
                        var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array, 
                                                                              n.array)
                } else {
                        var.pooled <- apply(dat - t(design %*% B.hat), 1, var, 
                                            na.rm = T)
                }
                stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
                if (!is.null(design)) {
                        tmp <- design
                        tmp[, c(1:n.batch)] <- 0
                        stand.mean <- stand.mean + t(tmp %*% B.hat)
                }
                s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, 
                                                                         n.array)))
                cat("Fitting L/S model and finding priors\n")
                batch.design <- design[, 1:n.batch]
                if (!NAs) {
                        gamma.hat <- solve(t(batch.design) %*% batch.design) %*% 
                                t(batch.design) %*% t(as.matrix(s.data))
                } else {
                        gamma.hat = apply(s.data, 1, Beta.NA, batch.design)
                }
                delta.hat <- NULL
                for (i in batches) {
                        delta.hat <- rbind(delta.hat, apply(s.data[, i], 1, var, 
                                                            na.rm = T))
                }
                gamma.bar <- apply(gamma.hat, 1, mean)
                t2 <- apply(gamma.hat, 1, var)
                a.prior <- apply(delta.hat, 1, aprior)
                b.prior <- apply(delta.hat, 1, bprior)
                if (prior.plots & par.prior) {
                        par(mfrow = c(2, 2))
                        tmp <- density(gamma.hat[1, ])
                        plot(tmp, type = "l", main = "Density Plot")
                        xx <- seq(min(tmp$x), max(tmp$x), length = 100)
                        lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
                        qqnorm(gamma.hat[1, ])
                        qqline(gamma.hat[1, ], col = 2)
                        tmp <- density(delta.hat[1, ])
                        invgam <- 1/rgamma(ncol(delta.hat), a.prior[1], b.prior[1])
                        tmp1 <- density(invgam)
                        plot(tmp, typ = "l", main = "Density Plot", ylim = c(0, 
                                                                             max(tmp$y, tmp1$y)))
                        lines(tmp1, col = 2)
                        qqplot(delta.hat[1, ], invgam, xlab = "Sample Quantiles", 
                               ylab = "Theoretical Quantiles")
                        lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
                        title("Q-Q Plot")
                }
                gamma.star <- delta.star <- NULL
                if (par.prior) {
                        cat("Finding parametric adjustments\n")
                        for (i in 1:n.batch) {
                                temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i, 
                                                                                 ], delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i], 
                                               b.prior[i])
                                gamma.star <- rbind(gamma.star, temp[1, ])
                                delta.star <- rbind(delta.star, temp[2, ])
                        }
                } else {
                        cat("Finding nonparametric adjustments\n")
                        for (i in 1:n.batch) {
                                temp <- int.eprior(as.matrix(s.data[, batches[[i]]]), 
                                                   gamma.hat[i, ], delta.hat[i, ])
                                gamma.star <- rbind(gamma.star, temp[1, ])
                                delta.star <- rbind(delta.star, temp[2, ])
                        }
                }
                cat("Adjusting the Data\n")
                bayesdata <- s.data
                j <- 1
                for (i in batches) {
                        bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, 
                                                                           ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1, 
                                                                                                                               n.batches[j])))
                        j <- j + 1
                }
                bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, 
                                                                      n.array)))) + stand.mean
                return(bayesdata)
        }
#+++++++++++++++++++++++++++ ComBat +++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ aprior +++++++++++++++++++++++++++++++++
aprior<-
        function (gamma.hat) 
        {
                m = mean(gamma.hat)
                s2 = var(gamma.hat)
                (2 * s2 + m^2)/s2
        }
#+++++++++++++++++++++++++++ aprior +++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ bprior +++++++++++++++++++++++++++++++++
bprior<-
        function (gamma.hat) 
        {
                m = mean(gamma.hat)
                s2 = var(gamma.hat)
                (m * s2 + m^3)/s2
        }
#+++++++++++++++++++++++++++ bprior +++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ it.sol +++++++++++++++++++++++++++++++++
it.sol<-
        function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 1e-04) 
        {
                n <- apply(!is.na(sdat), 1, sum)
                g.old <- g.hat
                d.old <- d.hat
                change <- 1
                count <- 0
                while (change > conv) {
                        g.new <- postmean(g.hat, g.bar, n, d.old, t2)
                        sum2 <- apply((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, 
                                      1, sum, na.rm = T)
                        d.new <- postvar(sum2, n, a, b)
                        change <- max(abs(g.new - g.old)/g.old, abs(d.new - d.old)/d.old)
                        g.old <- g.new
                        d.old <- d.new
                        count <- count + 1
                }
                adjust <- rbind(g.new, d.new)
                rownames(adjust) <- c("g.star", "d.star")
                adjust
        }
#+++++++++++++++++++++++++++ it.sol +++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ postmean +++++++++++++++++++++++++++++++++
postmean<-
        function (g.hat, g.bar, n, d.star, t2) 
        {
                (t2 * n * g.hat + d.star * g.bar)/(t2 * n + d.star)
        }
#+++++++++++++++++++++++++++ postmean +++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ postvar +++++++++++++++++++++++++++++++++
postvar<-
        function (sum2, n, a, b) 
        {
                (0.5 * sum2 + b)/(n/2 + a - 1)
        }
#+++++++++++++++++++++++++++ postvar +++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ make.list.of.processed.elists +++++++++++++++++++++++++++++++++++++++
make.list.of.processed.elists<-function(elists.list=NULL,
                                        controls.elists.list=NULL,
                                        gpr.path=NULL,
                                        targets.path=NULL,
                                        output.path.list=NULL,
                                        norm.methods=NULL,
                                        return_orig=FALSE){
        
        processed.elists.list<-vector("list",length(elists.list))
        
        for (index in 1:length(elists.list)){
                #it is possible to return original
                #not normalized, but logged elist
                if (return_orig){
                        norm.methods<-
                                c("none", norm.methods)
                }  
                        
                #declare list of lists
                processed.elists.list[[index]]<-
                        vector("list",length(norm.methods))
                
                names(processed.elists.list[[index]])<-
                        c(norm.methods)
                
                #perfrom normalizations 
                for (methodN in 1:length(norm.methods)){
                        processed.elists.list[[index]][[methodN]]<-
                                normalizeArrays(method=norm.methods[methodN],
                                                elist=elists.list[[index]], 
                                                controls.elist=controls.elists.list[[index]],
                                                gpr.path=gpr.path,
                                                targets.path=targets.path,
                                                output.path=output.path.list[[index]],
                                                contr.names=NULL)
                }
        } 
        
        return(processed.elists.list)
}
#+++++++++++++++++++++++++++ make.list.of.processed.elists +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ make.batchAdjust +++++++++++++++++++++++++++++++++++++++
make.batchAdjust<-function(list.of.elists=NULL,
                           is.logged=TRUE){
        for(elistNum in 1:length(list.of.elists)){
                list.of.elists[[elistNum]]<-
                        batchAdjust(elist=list.of.elists[[elistNum]],
                                    is.logged=is.logged)
                
        }
        return(list.of.elists)
}
#+++++++++++++++++++++++++++ make.batchAdjust +++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++ make.norm.tables +++++++++++++++++++++++++++++++++++++++
make.norm.tables<-function(all.elists.list=NULL,
                           description.colname="Description",
                           output.path=output.path,
                           tag=NULL){
        
        
        for (index in 1:length(all.elists.list)){
                elist<-all.elists.list[[index]]
                
                #name tag w/o extension
                filename.noext=paste0(output.path,"/",
                                     Sys.Date(),
                                     "_PAA_",
                                     tag, "_",
                                     names(all.elists.list)[index])
                
                output.table<-cbind.data.frame(Description=elist$genes[,description.colname],
                                                  elist$E)
                
                #write to excel file
                openxlsx:::write.xlsx(output.table,
                                      paste0(filename.noext,
                                             "-log.xlsx"))
                #unlog and write to file
                output.table.unlog<-output.table
                output.table.unlog[,2:ncol(output.table.unlog)]<-2^elist$E
                openxlsx:::write.xlsx(output.table.unlog,
                                      paste0(filename.noext,
                                             "-unlog.xlsx"))
        }
}
#+++++++++++++++++++++++++++ make.norm.tables +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ diff.analysis +++++++++++++++++++++++++++++++++++++++
diff.analysis<-function(all.elists.list=NULL,
                        groups=NULL,
                        is.logged=TRUE,
                        p_value_plots=FALSE,
                        volcano_plots=TRUE,
                        foldChange=foldChange,
                        pvalue=pvalue,
                        sortBY="fold_change",
                        output.path=output.path,
                        include.mMs=TRUE){
        
        #declare a list to store diff analysis features
        diffAnalysis.list<-vector("list", length(all.elists.list))
        names(diffAnalysis.list)<-names(all.elists.list)
        
        #UPD 2015-07-15: mMs fixed
        if (include.mMs){
                plot.method.options<-c("tTest","mMs")
        } else {
                plot.method.options<-"tTest"
        }
       
        #make t-Test and mMs volcano plots
        for (index in 1:length(all.elists.list)){
                
                #unlog the elist if needed
                elist<-all.elists.list[[index]]
                if (is.logged){
                        elist$E<-2^elist$E
                }
                
                #make the volcano plots: one for each method
                #and pvalue plots: two for each method (no adj and adj)
                for (method in plot.method.options){
                        #make the file tag
                        tag<-paste0("-", method, "-",
                                    names(all.elists.list)[index])
                        #volcano plots
                        volcanoPlot(elist=elist,
                                    group1=groups[[1]],
                                    group2=groups[[2]],
                                    method = method,
                                    p.thresh = pvalue,
                                    fold.thresh = foldChange,
                                    output.path = output.path,
                                    tag = tag,
                                    above=500,
                                    between=200)
                        #pvalue plots
                        if(p_value_plots){
                                for (adj in c(TRUE,FALSE)){
                                        pvaluePlot(elist=elist,
                                                   group1=groups[[1]],
                                                   group2=groups[[2]],
                                                   method = method,
                                                   adjust=adj,
                                                   output.path = output.path,
                                                   tag = tag,
                                                   above=500,
                                                   between=200)
                                }
                        }
                }
                #prepare data for univariate diff analysis
                #data is currently UN-logged!
                elist.E<-elist$E
                #add BRC identifiers to the elist.E copy
                rownames(elist.E)<-paste(elist$genes$Block, 
                                         elist$genes$Row,
                                         elist$genes$Column)
                
                #perform univariate analysis
                #and store the data in the list
                if (include.mMs){
                        diffAnalysis.list[[index]] <- 
                                diffAnalysis(input=elist.E,
                                             label1=groups[[1]], 
                                             label2=groups[[2]],
                                             class1=names(groups)[1], 
                                             class2=names(groups)[2],
                                             sortBY=sortBY,
                                             above=500,
                                             between=200,
                                             output.path=output.path,
                                             filename = paste0("diffAnalysis-tTest-mMs_",
                                                               names(all.elists.list)[index],
                                                               ".txt"))
                } else {
                        diffAnalysis.list[[index]] <- 
                                diffAnalysis.tTest(input=elist.E,
                                                   label1=groups[[1]], 
                                                   label2=groups[[2]],
                                                   class1=names(groups)[1], 
                                                   class2=names(groups)[2],
                                                   sortBY=sortBY,
                                                   output.path=output.path,
                                                   filename = paste0("diffAnalysis-tTest-",
                                                                     names(all.elists.list)[index],
                                                                     ".txt"))
                }
                
                
        }
        return(diffAnalysis.list)
}
#+++++++++++++++++++++++++++ diff.analysis +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ selectResults +++++++++++++++++++++++++++++++++++++++
selectResults<-function(all.elists.list = NULL,
                        diffAnalysis.list = NULL,
                        foldChange = 3,
                        pvalue = 0.001,
                        groups=NULL,
                        output.path=NULL,
                        sel.to.return="tTest"){
        require(openxlsx)
        
        selected.list<-vector("list", 
                              length(all.elists.list))
        names(selected.list)<-names(all.elists.list)
        
        for (index in 1:length(all.elists.list)){
                
                elist<-all.elists.list[[index]]
                results.table<-diffAnalysis.list[[index]]
                
                selectedFeatures.tTest<-
                        results.table$BRC[results.table$tTest<pvalue &
                                                  (results.table$fold_change>foldChange |
                                                           results.table$fold_change<1/foldChange)]
                if(length(selectedFeatures.tTest)==0){
                        cat(c("tTest: no features could be selected for",
                              names(all.elists.list)[index],
                              "case, ouput path:\n",
                              output.path,
                              "\n"))
                        selected.tTest<-NULL
                } else {
                        selected.tTest<-
                                selectFeatures.plotHeatmap(features=selectedFeatures.tTest, 
                                                           elist=elist, 
                                                           diffAnalysis.results=results.table,
                                                           columns1=groups[[1]], 
                                                           columns2=groups[[2]],
                                                           output.path=output.path, 
                                                           filename=paste0("biomarker_heatmap-tTest-",
                                                                           names(all.elists.list)[index],
                                                                           ".pdf"),
                                                           NM_index=TRUE)
                        
                        openxlsx:::write.xlsx(selected.tTest,
                                              file = paste0(output.path,
                                                            "/results-tTest-FC_",
                                                            foldChange,
                                                            "-pval_",
                                                            pvalue,"-",
                                                            names(all.elists.list)[index],
                                                            ".xlsx"))
                }
                
                selectedFeatures.mMs<-
                        results.table$BRC[results.table$mMs<pvalue &
                                                  (results.table$fold_change>foldChange |
                                                           results.table$fold_change<1/foldChange)]
                if(length(selectedFeatures.mMs)==0){
                        cat(c("mMs: no features could be selected for",
                              names(all.elists.list)[index],
                              "case, ouput path:\n",
                              output.path,
                              "\n"))
                        selected.mMs<-NULL
                } else {
                        
                        selected.mMs<-
                                selectFeatures.plotHeatmap(features=selectedFeatures.mMs, 
                                                           elist=elist, 
                                                           diffAnalysis.results=results.table,
                                                           columns1=groups[[1]], 
                                                           columns2=groups[[2]],
                                                           output.path=output.path, 
                                                           filename=paste0("biomarker_heatmap-mMs-",
                                                                           names(all.elists.list)[index],
                                                                           ".pdf"),
                                                           NM_index=TRUE)
                        
                        openxlsx:::write.xlsx(selected.mMs,
                                              file = paste0(output.path,
                                                            "/results-mMs-FC_",
                                                            foldChange,
                                                            "-pval_",
                                                            pvalue,"-",
                                                            names(all.elists.list)[index],
                                                            ".xlsx"))
                }
                if(sel.to.return=="tTest"){
                        selected.list[[index]]<-selected.tTest
                        
                } else {
                        selected.list[[index]]<-selected.mMs
                        
                }
        }
        return(selected.list)
}
#+++++++++++++++++++++++++++ selectResults +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ batchFilter.anova +++++++++++++++++++++++++++++++++++++++
batchFilter.anova <- function(elist=NULL, 
                              p.thresh=0.05,
                              fold.thresh=1.5, 
                              output.path=NULL){
        
        if(is.null(elist)) {
                stop("ERROR: Not all mandatory arguments have been defined!")
        }
        row.number <- nrow(elist$E)
        col.number <- ncol(elist$E)
        discard <- c()
        output <- matrix(nrow=row.number, ncol=4)
        colnames(output) <- c("BRC", "ANOVA p-value", "Fold Change", "Description")
        #how many batches and batch combinations are there?
        batches<-unique(elist$targets$Batch)
        Nbatches<-length(batches)
        cat(c("Found batches:", Nbatches, "\n"))
        #if less than two batches, return null
        if (Nbatches<2) {return(NULL)}
        
        batch.combs<-combn(Nbatches,2)
        
        
        for(zeile in 1:row.number) {
                #make df for anova analysis
                data<-data.frame(Value=elist$E[zeile,],
                                 Batch=elist$targets$Batch)
                #perform anova on readings vs batch
                
                p<-oneway.test(Value~Batch,data)$p.value
                if(is.na(p)) {p<-1}
          
                #get the means for each batch
                means<-rep(NA,Nbatches)
                for (batchInd in 1:Nbatches){
                        means[batchInd]<-mean(data[data$Batch==batches[batchInd],"Value"])
                }
                
                #cycle through all combinations, find fold changes
                #and store them in a separate vector
                FCvector<-rep(NA,ncol(batch.combs))
                for (batchInd in 1:ncol(batch.combs)){
                        mean1<-means[batch.combs[1,batchInd]]
                        mean2<-means[batch.combs[2,batchInd]]
                        #determine FC
                        FCvector<-mean1/mean2
                        
                }
                #if max FC value > then 1/minFC,
                #take it as FC, else - otherwise
                FC<-ifelse(max(FCvector)>1/min(FCvector),
                           max(FCvector),
                           min(FCvector))
                
                #prep output table
                output[zeile,] <- c(paste0(elist$genes[zeile,1], ",",
                                           elist$genes[zeile,3], ",", 
                                           elist$genes[zeile,2]), 
                                    p, 
                                    FC,
                                    elist$genes[zeile,4])
                #determine if the feature (row) is to be rejected
                if(p < p.thresh && (FC > fold.thresh || FC < 1/fold.thresh)) {
                        discard <- c(discard,zeile)
                }
        }
        message(paste0("batchFilter - number of features to discard: ",
                       length(discard)), "\n")
        
        if(is.null(output.path)) {
                plot(x=log2(as.numeric(output[,3])), 
                     y=-log10(as.numeric(output[,2])),
                     xlab="log2(fold change)", 
                     ylab="-log10(p-value)",
                     main="batch filter volcano")
                
                if(length(discard)>0){
                        points(x=log2(as.numeric(output[discard,3])),
                               y=-log10(as.numeric(output[discard,2])), 
                               col="red")
                }
                
                abline(h=-log10(p.thresh), lty=2, col="red")
                abline(v=log2(fold.thresh), lty=2, col="red")
                abline(v=log2(1/fold.thresh), lty=2, col="red")
                if(length(discard)>0){
                        elist <- elist[-discard,]
                }
        }else{
                write.table(x=output, 
                            file=paste(output.path,
                                       "/batch_filter.txt",sep=""), 
                            sep="\t", eol="\n",
                            row.names=FALSE,
                            quote=FALSE)
                write.table(x=output[discard,], 
                            file=paste(output.path,
                                       "/batch_filter_discarded.txt",sep=""), sep="\t", eol="\n",
                            row.names=FALSE, quote=FALSE)
                
                tiff(paste(output.path,"/batch_filter_before.tiff",sep=""), 
                     width=2000,
                     height=2000,
                     pointsize=15, 
                     compression="lzw", 
                     res=300)
                
                plot(x=log2(as.numeric(output[,3])),
                     y=-log10(as.numeric(output[,2])),
                     xlab="log2(fold change)",
                     ylab="-log10(p-value)", 
                     main="batch filter volcano")
                
                if(length(discard)>0){
                        points(x=log2(as.numeric(output[discard,3])),
                               y=-log10(as.numeric(output[discard,2])), 
                               col="red")
                }
                abline(h=-log10(p.thresh), lty=2, col="red")
                abline(v=log2(fold.thresh), lty=2, col="red")
                abline(v=log2(1/fold.thresh), lty=2, col="red")
                dev.off()
                
                if(length(discard)>0){
                        tiff(paste(output.path,"/batch_filter_after.tiff",sep=""),
                             width=2000,
                             height=2000, 
                             pointsize=15,
                             compression="lzw",
                             res=300)
                        
                        plot(x=log2(as.numeric(output[-discard,3])),
                             y=-log10(as.numeric(output[-discard,2])),
                             xlab="log2(fold change)",
                             ylab="-log10(p-value)",
                             main="batch filter volcano")
                        dev.off()
                        
                }
        }
        return(discard)
}
#+++++++++++++++++++++++++++ batchFilter.anova +++++++++++++++++++++++++++++++++++++++

#------------------------- calc.mMs.pval --------------------------------------------------------
calc.mMs.pval<-function(vector1=NULL,
                        vector2=NULL,
                        mMs.matrix=NULL,
                        above=1500,
                        between=400){
        #length of vector 2
        vect2_len<-length(vector2)
        
        #declare vector to stor p-values
        pval_vect<-rep(NA,vect2_len)
        
        #sort vector2 (smallest to largest)
        vect2.sorted<-sort(vector2)
        
        #length of vector2 is max M stat "order"
        #length of vector1 is max M stat M* value
        for (order in 1:vect2_len){
                
                #calculate M* value
                Mstar<-mCountRancho(vector1=vector1,
                                    vector2=vect2.sorted[1:(vect2_len-order+1)],
                                    above=above,
                                    between=between)
                #locate appropriate p-value for a given M* and order
                pval_vect[order]<-mMs.matrix[order,Mstar+1]
        }
        #return minimum p-value
        return(min(pval_vect))
}
#------------------------- calc.mMs.pval --------------------------------------------------------

#------------------------- mCountRancho --------------------------------------------------------
mCountRancho <- function(vector1=NULL, 
                         vector2=NULL, 
                         above=1500,
                         between=400) {
        if(is.null(vector1) || 
           is.null(vector2)) {
                stop("ERROR: Not all mandatory arguments have been defined!")
        }
        #add the "between" value to all values of vector2
        vector2.mod <- vector2+between
        
        #subtract "above" values (arbitrary noise threshold)
        #from both vector1 and modified vector2
        #finally, find ratio of vector 1 to max of vector2
        #(if max vector 2 value is less than 1 - return 1)
        ratio <- (vector1-above)/max(vector2.mod-above,1)
        
        #final count is a number of members of vector1
        #which are larger than max value in vector2
        count <- length(ratio[ratio >= 1])
        
        return(count)
}
#------------------------- mCountRancho --------------------------------------------------------

#------------------------- mMsMatrixRancho --------------------------------------------------------
mMsMatrixRancho <- function(n1=NULL, 
                            n2=NULL){
        if(is.null(n1) || 
           is.null(n2)) {
                stop("ERROR: Not all mandatory arguments have been defined!")
        }
        #there are 
        #m=n1 white balls in the urn
        #and n=n2 black balls
        #we can draw k balls (more on that later)
        #amongst those we draw q=n1-M* white balls
        #because at M*=0 we can *expect* there are still n1 white balls
        #and at M*=n1 there are definitely no white balls
        #now it is only reasonable that we get q white balls 
        #in q+1 or more draws, with max draws capping at q+n2
        #thus total draws k=(q+1):(q+n2) because at order 1 we can draw n2 times
        #and at order n2 we can draw only 1 time (i.e. members of group 2 to compare)
        
        #make the matrix
        output.matrix<-mapply(function(q,m,n){
                k<-(q+n2):(q+1)
                phyper(q,m,n,k)
                },
                q=n1:0,
                m=n1,
                n=n2)
        return(output.matrix)
}
#------------------------- mMsMatrixRancho --------------------------------------------------------

#####################################################################################################
#################################### added by Rancho ################################################
#####################################################################################################
