library(limma)
library(DESeq2)
library(edgeR)

##_______________________Limma Voom___________________________________ 
#process count data 

myCPM <- cpm(countDF)
CPM <- cpm(countDF)
thresh <- myCPM > 0.5
keep <- rowSums(thresh) >= 2 #try different filters 
keep.exprs <- rowSums(CPM>1)>=3 #filter2
counts.keep <- countDF[keep,]
counts_keep <- countDF[keep.exprs,]
plot(myCPM[,1],countDF[,1])
y <- DGEList(counts.keep)
y2 <- DGEList(counts_keep)
y <- calcNormFactors(y)
y2 <- calcNormFactors(y2)

plotMDS(y) #examine processed data
plotMDS(y2)

#_____________________________________________________________________
design_matrix <-model.matrix(~ 0 + targets2$Factor) #create design matrix using targets file 
colnames(design_matrix) <- levels(targets$Factor)
#create design matrix

v <- voom(y,design_matrix, plot = TRUE)
v2 <- voom(y2,design_matrix, plot = TRUE)
fit <- lmFit(v)
fit2 <- lmFit(v2)

#create contrasts for specific hypothese 
cont.matrix <- makeContrasts(
                            M1VsV1 = V1 - M1, 
                             M1VsA1 = M1 - A1, 
                             V1VsA1 = V1 - A1,
                             M6VsV6 = M6 - V6, 
                             M6VsA6 = M6 - A6, 
                             V6VsA6 = V6 - A6,
                             M12VsV12 = M12 - V12, 
                             M12VsA12 = M12 - A12, 
                             V12VsA12 = V12 - A12,
                             levels= design_matrix
                            )

#____________ Fit model using different filtering criteria
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont2 <- contrasts.fit(fit2, cont.matrix)
fit.cont <- eBayes(fit.cont)
fit.cont2 <- eBayes(fit.cont2)
summa.fit <- decideTests(fit.cont, p.value=0.1, adjust.method = "BH")
summa.fit2 <- decideTests(fit.cont2, p.value=0.1, adjust.method = "BH")
summary(summa.fit)
summary(summa.fit2)

vennDiagram(summa.fit[,7:9], circle.col = c("Blue", "Red", "Green")) #Venn Diagram 12 Hour, up and down = TRUE 

##write  tables to working directory 
for(i in colnames(cont.matrix)){
  Sig.table <- topTable(fit.cont,coef= i,sort.by="p",n="Inf", adjust = "BH")
  Sig.table <- na.omit(Sig.table)
  write.csv(Sig.table, paste0(i,".csv", collapse = ""))
}


###DESeq2_____________________________________________________________ 

cdata <- data.frame(Group = factor(targets$Factor))
rownames(cdata) <- colnames(countDF)



describe <- function(cdata, countDF, x){
  
  deg <- countDF[rowSums(countDF > 3) > 2, ]
  d <- DESeqDataSetFromMatrix(deg, colData = cdata , design = ~ Group)
  d <- DESeq(d)
  vsd <- getVarianceStabilizedData(d)
  heatmap(cor(vsd), cexCol = .75, cexRow = ,75) 
  b <- prcomp(t(vsd))
                              
  return(list(d, b))
}


describe(d)

##calculate LFCs using padj = <= .1 criteria to filter the gene list, results displayed in a volcano plot
LFCsigtests <- function(d, x = c("Group","M1","V1")){
  
  comp1 <-as.data.frame(results(d, contrast=x))
  comp1$Gene <- rownames(comp1)
  rownames(comp1) <- NULL
  significant <- comp1[which(comp1$padj <=.1), ]
  plot(comp1$log2FoldChange, -log(comp1$padj), pch = 10, cex = .2, main = x, xlab = "Log2 Fold Change", ylab = "-log(adj p value)")
  points(significant$log2FoldChange, -log(significant$padj), col = "red", cex = .2, pch = 10 )
  write.csv(significant, paste0("DESeq2", x[2], "_", x[3], ".csv", collapse = ""))
  
  return(comp1)
  
} 

par(mfrow = c(3,3))
LFCsigtests(d, c("Group","M1","V1")) #generate plot for comparison, save file to working directory 
LFCsigtests(d, c("Group","M1","A1"))
LFCsigtests(d, c("Group","A1","V1"))

LFCsigtests(d, c("Group","M6","V6"))
LFCsigtests(d, c("Group","M6","A6"))
LFCsigtests(d, c("Group","A6","V6"))

LFCsigtests(d, c("Group","M12","V12"))
LFCsigtests(d, c("Group","M12","A12"))
LFCsigtests(d, c("Group","A12","V12"))


test <- c("M1VV1", "M1VA1", "A1VV1")

## Overlap: 1vs2, 1vs3, 2vs3, 1vs2vs3
Getcounts <- function(df,df2,df3) {
  
  data <-  data.frame( x = nrow(df1), y = nrow(df2), x = nrow(df3))
  One <-  intersect(df1$Gene, df2$Gene)
  Two <- intersect(df1$Gene, df3Gene)
  Three <- intersect(df2$Gene, df3$Gene)
  all <- intersect(intersect(df1$Gene, df2$Gene), df3$Gene)
  return(list(counts = data, Overlap = data.frame(One = One, Two = Two, Three = Three, all = all )))

}


###______________EdgeR______________________________________________ 
edgeRModel <- DGEList(counts = countDF )
edgeRModel <- edgeRModel$counts[rowSums(edgeRmodel$counts > 3) > 2, ]
edgeRModel <- calcNormFactors(edgeRModel, method = "TMM")
design.matrix <- model.matrix(~0 + Group, data = cdata) 
edgeRModel <- estimateGLMCommonDisp(edgeRmodel, design.matrix)
edgeRModel <- estimateGLMTrendedDisp(edgeRModel, design.matrix)
edgeRModel <- estimateGLMTagwiseDisp(edgeRModel, design.matrix)


GLM1 <- glmFit(edgeRModel, design.matrix)
colnames(coef(GLM1))
#"GroupA1"  "GroupA12" "GroupA6"  "GroupM1"  "GroupM12" "GroupM6"  "GroupV1"  "GroupV12" "GroupV6" 
 
contrastsEdge <- as.data.frame(makeContrasts(
                                             M1VsV1=  GroupV1 - GroupM1, #calculate contrasts
                                             M1VsA1 = GroupM1 - GroupA1, 
                                             V1VsA1 = GroupV1 - GroupA1,
                                             M6VsV6 = GroupM6 - GroupV6, 
                                             M6VsA6 = GroupM6 - GroupA6, 
                                             V6VsA6 = GroupV6 - GroupA6,
                                             M12VsV12 = GroupM12 - GroupV12, 
                                             M12VsA12 = GroupM12 - GroupA12, 
                                             V12VsA12 = GroupV12 - GroupA12,
                                             levels=design.matrix)
                                            )

##use contrast data frame and GLMfit to get all LFCs with threshold of adj p = .1, files written to wd
GetEdgeRSig <- function(ConFrame, GLMmodel) {
  
  
  for(i in colnames(contrastsEdge)){
    LRT <- glmLRT(GLM1, contrast = c(contrastsEdge[,i])) #specify custom contrasts 
    resultsOverview <- topTags(LRT, n = "inf", adjust.method="BH", sort.by = "PValue", p.value = .1)$table
    resultsOverview <- as.data.frame(na.omit(resultsOverview))
    write.csv(resultsOverview, file = paste0("EdgeR",i,".csv", collapse = ""))
    
    }
}

GetEdgeRSig(contrastsEdge, GLM1)



##______________________________________ROC curves 

#load true positive files for filtering; genes significant in all three approaches used in the 2013 paper
#variables are in excel sheets, so must be extracted first 
Hour12C_AVR <- read_excel(file.choose(), sheet = 1)
Hour12C_VIR <- read_excel(file.choose(), sheet = 2)
Hour12AVR_VIR <- read_excel(file.choose(), sheet = 3)


SignIndex <- c("Gene", "ER_Class.Sig", "ER_GLM.Sig", 
               "Cufflinks.Sig", "Total.Sig", "All.Sig")


Hour12C_AVR <- Hour12C_AVR[,SignIndex]
Hour12C_VIR <- Hour12C_VIR[,SignIndex]
Hour12AVR_VIR <- Hour12AVR_VIR[,SignIndex]

#conversion to .csv to minimize file complexity 
#write.csv(Hour12C_AVR, file = "Hour12C_AVR.csv", row.names=FALSE)
#write.csv(Hour12C_VIR, file = "Hour12C_VIR.csv", row.names=FALSE)
#write.csv(Hour12C_VIR, file = "Hour12AVR_VIR.csv", row.names=FALSE)

####--------------------------------------------------------laod files in for analyses
#Hour1C_AVR <- read.csv( file = "Hour1C_AVR.csv", header = TRUE, stringsAsFactors = FALSE)
#Hour1C_VIR <- read.csv( file = "Hour1C_VIR.csv", header = TRUE, stringsAsFactors = FALSE)
#Hour1AVR_VIR <- read.csv(file = "Hour1AVR_VIR.csv", header = TRUE, stringsAsFactors = FALSE )


#Hour6C_AVR <- read.csv( file = "Hour6C_AVR.csv", header = TRUE, stringsAsFactors = FALSE)
#Hour6C_VIR <- read.csv( file = "Hour6C_VIR.csv", header = TRUE, stringsAsFactors = FALSE)
#Hour6AVR_VIR <- read.csv(file = "Hour6AVR_VIR.csv", header = TRUE, stringsAsFactors = FALSE )

Hour12C_AVR <- read.csv( file = "Hour12C_AVR.csv", header = TRUE, stringsAsFactors = FALSE)
Hour12C_VIR <- read.csv( file = "Hour12C_VIR.csv", header = TRUE, stringsAsFactors = FALSE)
Hour12AVR_VIR <- read.csv(file = "Hour12AVR_VIR.csv", header = TRUE, stringsAsFactors = FALSE )

Voom_12C_AVR <- read.csv( file = file.choose(), header = TRUE, stringsAsFactors = FALSE)
Voom12C_VIR <- read.csv( file = file.choose(), header = TRUE, stringsAsFactors = FALSE)
Voom12AVR_VIR <- read.csv(file = file.choose(), header = TRUE, stringsAsFactors = FALSE )

GetROC <- function(DF, DF2, main, x) {
  
require(ROCR)
  
  colnames(DF2) <- c("Gene",	"logFC",	"AveExpr", 	"t",	"P.Value",	"adj.P.Val", "B")
  DF$Gene <- toupper(DF$Gene)
  DF2 <- DF2[which(DF2$adj.P.Val <= as.numeric(x)), ]
  DF <- DF[DF$Gene %in% DF2$Gene, ]
  preds <- prediction((DF2$adj.P.Val), DF$All.Sig)
  ROCcurve <- performance(preds, measure = "sens")
  plot(ROCcurve, main = paste(x, main), col = "red" )
}


par(mfrow = c(2,2))
for ( i in c(1, .1, .05, .01) ){
GetROC(DF = Hour12C_VIR, DF2 = Voom12C_VIR, "\nVoom: M12 vs. A12", x = i)
}





GetDiscoverRates <- function(DF, DF2, labels = NULL, package, cutoff ) {

  require(gridExtra)
  
  colnames(DF2) <- c("Gene",	"logFC",	"AveExpr", 	"t",	"P.Value",	"adj.P.Val", "B")
  DataTable = data.frame(TPR = "", FPR = "")
  
  if (is.null(labels) == TRUE) {
    message("No labels included, defaulting to TPR and FPR")
  }
  else {
  colnames(DataTable) <- labels 
  DF$Gene <- toupper(DF$Gene)
  }
  
  DF2 <- DF2[which(DF2$adj.P.Val <= as.numeric(cutoff)), ]
  DF <- DF[DF$Gene %in% DF2$Gene, ]
  TP <- DF[which(DF$All.Sig == 1), ]
  FP <- DF[which(DF$All.Sig == 0), ]
  TPRate <- as.numeric(nrow(TP)) / sum(as.numeric(nrow(FP)), as.numeric(nrow(TP)))
  FPRate <- as.numeric(nrow(FP)) / sum(as.numeric(nrow(FP)), as.numeric(nrow(TP)))
  DataTable$TPR <- round(TPRate, 2)
  DataTable$FPR <- round(FPRate, 2)
  DataTable <- rbind(DataTable, c(as.numeric(nrow(TP)), as.numeric(nrow(FP))))
  rownames(DataTable) <- c(paste("Package:", package), "Totals")
  g <- tableGrob(DataTable)
  grid.newpage()
  grid.draw(g)
  justify(g,"right", "top")

}

par(mfrow = c(2,2))
for ( i in c(.1, .05, .01)){
GetDiscoverRates(DF, DF2, cutoff = i, package = "VOOM") 
}


