setwd("C:/Users/Himanshu/Downloads")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("DESeq")

biocLite("pasilla")
install.packages("DESeq")
datafile = system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )
datafile
pasillaCountTable = read.table( datafile, header=TRUE, row.names=1 )
head( pasillaCountTable )

pasillaDesign = data.frame(
   row.names = colnames( pasillaCountTable ),
   condition = c( "untreated", "untreated", "untreated",
                    "untreated", "treated", "treated", "treated" ),
   libType = c( "single-end", "single-end", "paired-end",
                  "paired-end", "single-end", "paired-end", "paired-end" ) )
pasillaDesign

pairedSamples = pasillaDesign$libType == "paired-end"
countTable = pasillaCountTable[ , pairedSamples ]
condition = pasillaDesign$condition[ pairedSamples ]
head(countTable)
condition
library("DESeq")
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )
cds = estimateDispersions( cds )
str( fitInfo(cds) )
# This plot shows the difference btw sampling variance and true variance btw different genes
# using the position of the estimates relative to the regression line.
plotDispEsts( cds )

head( fData(cds) )

# Comparison:
res = nbinomTest( cds, "untreated", "treated" )
head(res)
# plot of log2fold change against mean normalised counts: 
# red genes are significant at 10% false discovery rate. (untreated vs treated)
plotMA(res)

# histogram of p values:
# 
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")

resSig = res[ res$padj < 0.1, ]
head( resSig[ order(resSig$pval), ] )

# looking at the most strongly down-regulated genes.
head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )

# looking at the most strongly up-regulated genes.
head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )

write.csv( res, file="My Pasilla Analysis Result Table.csv" )

ncu = counts( cds, normalized=TRUE )[ , conditions(cds)=="untreated" ]

plotMA(data.frame(baseMean = rowMeans(ncu),  log2FoldChange = log2( ncu[,2] / ncu[,1] )), col = "black")

# partially without replicates
cdsUUT = cds[ , 1:3]
pData( cdsUUT )

cdsUUT = estimateSizeFactors( cdsUUT )
cdsUUT = estimateDispersions( cdsUUT )
resUUT = nbinomTest( cdsUUT, "untreated", "treated" )

plotMA(resUUT)

# without any replicates.
cds2 = cds[ ,c( "untreated3", "treated3" ) ]
cds2 = estimateDispersions( cds2, method="blind", sharingMode="fit-only" )
res2 = nbinomTest( cds2, "untreated", "treated" )
plotMA(res2)
addmargins( table( res_sig = res$padj < .1, res2_sig = res2$padj < .1 ) )

head( pasillaCountTable )
pasillaDesign
cdsFull = newCountDataSet(pasillaCountTable, pasillaDesign)
cdsFull
cdsFull = estimateSizeFactors( cdsFull )
cdsFull = estimateDispersions( cdsFull )
plotDispEsts( cdsFull )

fit1 = fitNbinomGLMs( cdsFull, count ~ libType + condition )
fit0 = fitNbinomGLMs( cdsFull, count ~ libType )

str(fit1)

pvalsGLM = nbinomGLMTest( fit1, fit0 )
padjGLM = p.adjust( pvalsGLM, method="BH" )
padjGLM

tab1 = table( "paired-end only" = res$padj < .1, "all samples" = padjGLM < .1 )
addmargins( tab1 )
table(sign(fitInfo(cds)$perGeneDispEsts - fitInfo(cdsFull)$perGeneDispEsts))
trsf = function(x) log( (x + sqrt(x*x+1))/2 ) 
plot( trsf(fitInfo(cds)$perGeneDispEsts),  trsf(fitInfo(cdsFull)$perGeneDispEsts), pch=16, cex=0.45, asp=1) > abline(a=0, b=1, col="red3")
head(fit1)

cdsFullB = newCountDataSet( pasillaCountTable, pasillaDesign$condition )
cdsFullB = estimateSizeFactors( cdsFullB )
cdsFullB = estimateDispersions( cdsFullB )
resFullB = nbinomTest( cdsFullB, "untreated", "treated" )
tab2 = table(
   `all samples simple` = resFullB$padj < 0.1,
   `all samples GLM` = padjGLM < 0.1 )
addmargins(tab2)

# 5 Independent Filtering and multiple testing:
rs = rowSums ( counts ( cdsFull ))
rs
theta = 0.4
use = (rs > quantile(rs, probs=theta))
table(use)
cdsFilt = cdsFull[ use, ]

# 7 Data quality assessment:
# Heatmap of count table:
cdsFullBlind = estimateDispersions( cdsFull, method = "blind" )
vsdFull = varianceStabilizingTransformation( cdsFullBlind )
install.packages("RColorBrewer")
library("RColorBrewer")
install.packages("gplots")
library("gplots")
select = order(rowMeans(counts(cdsFull)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(5,3))
heatmap.2(counts(cdsFull)[select,], col = hmcol, trace="none", margin=c(10,6))
