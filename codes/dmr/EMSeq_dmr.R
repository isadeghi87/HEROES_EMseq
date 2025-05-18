

setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/EMseq/')
# Main analysis package
library("methylKit")
# Annotation package
library("genomation")
library("GenomicRanges")
library(bsseq)
library(magrittr)
library(readr)
library(dmrseq)

## metdata
# Define the list containing the bismark coverage files.
files <- list.files(path ='./results/nextflow/bismark/methylation_calls/methylation_coverage',full.names = T)

samples = samples = sub(".*OE0290-PED_(.*)_1_val_1_bismark_bt2_pe\\.deduplicated\\.bismark\\.cov\\.gz", "\\1", files)
metadata = read_delim('./data/HEROES-EMseq-plasma.csv')
metadata$samples = paste0(metadata$`Patient ID`,"_",metadata$`Sample Type`)
metadata$samples = gsub('OE0290-PED_','',metadata$samples)
metadata$diagnosis = metadata$`INFORM patient ID`
metadata$diagnosis = gsub('non-oncological control','control',metadata$diagnosis)
metadata$diagnosis[metadata$diagnosis != 'control'] = 'tumor'
clinic = read.delim('/Users/i439h/Documents/Heroes_EMseq_clinical.csv',header = T,sep = ';')
clinic$samples = paste0(clinic$Patient.ID,"_",clinic$Sample.Type)
clinic$samples = gsub('OE0290-PED_','',clinic$samples)
clinic$diagnosis = gsub('non-oncological control','control',clinic$Diagnosis)

## match diagnosis for samples
methyldata = data.frame(samples = samples, files = files)
methyldata$diagnosis = metadata$diagnosis[match(methyldata$samples,metadata$samples)]
methyldata = na.omit(methyldata)
methyldata$treatment = ifelse(methyldata$diagnosis=='control',0,1)


## get intersect of samples present
id = intersect(samples,metadata$samples)
methyldata = methyldata[methyldata$samples %in% id,]

## merge with clinics
merged = dplyr::left_join(methyldata,clinic,"samples")
cov = as.data.frame(table(merged$diagnosis.y))
colnames(cov) = c('diagnosis','frequency')
library(ggplot2)
bp = ggplot2::ggplot(cov,aes(x = diagnosis,y = frequency))+
  geom_bar(stat = 'identity',aes(fill = diagnosis))+
  labs(x = '',y = '')

# read the listed files into a methylRawList object making sure the other
# parameters are filled in correctly.
dirs = as.list(methyldata$files)
ids = as.list(methyldata$samples)

## create methylation object
myobj <- methRead(dirs,
                  sample.id=ids,
                  assembly="hg38", # replace with your genome assembly
                  treatment= methyldata$treatment, # define control (0) and treatment (1) if applicable
                  pipeline="bismarkCoverage",
                  mincov = 2)


# check number of samples
head(myobj)

# Get a histogram of the methylation percentage per sample
# Here for sample 1
getMethylationStats(myobj[[1]], plot=TRUE, both.strands=FALSE)

# Get a histogram of the read coverage per sample
getCoverageStats(myobj[[1]], plot=TRUE, both.strands=FALSE)

# Get percentile data by setting plot=FALSE
getCoverageStats(myobj[[1]], plot=T, both.strands=FALSE)
getCoverageStats(myobj[[2]], plot=T, both.strands=FALSE)


# ## filter 
myobj.filt <- filterByCoverage(myobj,
                               lo.count=2,
                               lo.perc=NULL,
                               hi.count=NULL,
                               hi.perc=99.9)
## normalize 
myobj.filt.norm <- normalizeCoverage(myobj, method = "median")
meth <- unite(myobj.filt.norm, destrand=F)
meth
saveRDS(meth,file = './data/processed/meth_unite.rds')
# get percent methylation matrix
pm=percMethylation(meth)

# Hierarchical clustering
hc <- hclust(dist(pm))
plot(hc)


# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)

# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
hist(sds, breaks = 100)

# keep only CpG with standard deviations larger than 2%
meth <- meth[sds > 2]

# meth_filt <- meth[sds > 2]

# This leaves us with this number of CpG sites
# nrow(meth_filt)

getCorrelation(meth,plot=TRUE)

## cluster samples 
pdf("./results/figures/dmr/cpg_methyl_clustering.pdf",width = 10,height = 20)
clusterSamples(meth, dist="correlation", method = 'ward', plot=TRUE)
dev.off()

## pca
PCASamples(meth, screeplot=TRUE)

pdf("./results/figures/dmr/pca_samples.pdf",width =10,height = 7)
PCASamples(meth)
dev.off()

## calculate DMRs
# Test for differential methylation... This might take a few minutes.
myDiff <- calculateDiffMeth(meth)
myDiff

# Simple volcano plot to get an overview of differential methylation
pdf('./results/figures/dmr/dmr_volcano.pdf',width = 8,height = 6)
plot(myDiff$meth.diff, -log10(myDiff$qvalue))
abline(v=0)
dev.off()


## get data
dat = getData(myDiff)
diffdat = subset(dat,pvalue<0.05& abs(meth.diff)>5)

# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hyper")
nrow(myDiff25p.hyper)

# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=5,qvalue=0.05,type="hypo")
nrow(myDiff25p.hypo)

# get all differentially methylated bases
myDiff25p <- getMethylDiff(myDiff,difference=5,qvalue=0.05)
nrow(myDiff25p)

# First load the annotation data; i.e the coordinates of promoters, TSS, intron and exons
refseq_anot <- readTranscriptFeatures("data/annotations/hg38_refseq_genes.bed")

# Annotate hypermethylated CpGs ("target") with promoter/exon/intron
# information ("feature"). This function operates on GRanges objects, so we # first coerce the methylKit object to GRanges.
myDiff25p.annot <- annotateWithGeneParts(target = as(myDiff25p,"GRanges"),
                                              feature = refseq_anot)

# Summary of target set annotation
myDiff25p.annot

# View the distance to the nearest Transcription Start Site; the target.row column in the output indicates the row number in the initial target set
dist_tss <- getAssociationWithTSS(myDiff25p.annot)
head(dist_tss)

# See whether the differentially methylated CpGs are within promoters,introns or exons; the order is the same as the target set
getMembers(myDiff25p.annot)


pdf(file = './results/figures/dmr/dmr_pie.pdf',width = 5,height = 5)
plotTargetAnnotation(myDiff25p.annot, main = "Differential Methylation Annotation")
dev.off()

# Adjust the distance for promoter definition, e.g., 2000 bp upstream
promoter_region <- promoters(refseq_anot, upstream=2000, downstream=0)

# Annotate DMRs with the expanded promoter regions
myDiff25p_promoter_annot <- annotateWithGeneParts(target = as(myDiff25p,"GRanges"),
                                                  feature = promoter_region)

# View the annotation results
plotTargetAnnotation(myDiff25p_promoter_annot)



# This can also be summarized for all differentially methylated CpGs

# Load the CpG info
cpg_anot <- readFeatureFlank("./data/annotations/hg38_cpg_islands.txt", feature.flank.name = c("CpGi", "shores"), flank=2000)
diffCpGann <- annotateWithFeatureFlank(as(myDiff25p,"GRanges"), feature = cpg_anot$CpGi, flank = cpg_anot$shores, feature.name = "CpGi", flank.name = "shores")

# See whether the CpG in myDiff25p belong to a CpG Island or Shore
head(getMembers(diffCpGann))

# This can also be summarized for all differentially methylated CpGs
pdf(file = './results/figures/dmr/dmr_cpg_pie.pdf',width = 5,height = 5)
plotTargetAnnotation(diffCpGann, main = "Differential Methylated Annotation in CpGs")
dev.off()

# Summarize the original object counts over a certain region, here the CpG Islands
# You can ignore the warnings here...
myobj_islands <- regionCounts(myobj, cpg_anot$CpGi)

# Filter the summarized counts by coverage
myobj_islands_filt <- filterByCoverage(myobj_islands,
                                       lo.count=2,
                                       lo.perc=NULL,
                                       hi.count=NULL,
                                       hi.perc=99.9)
# Perform simple normalization
myobj_islands_filt_norm <- normalizeCoverage(myobj_islands_filt, method = "median")
# Merge the samples again
meth_islands <- unite(myobj_islands_filt_norm, destrand=FALSE)

# Test for differential methylation... This might take a few minutes.
myDiff_islands <- calculateDiffMeth(meth_islands)

# Rank by significance
myDiff_islands <- myDiff_islands[order(myDiff_islands$pvalue),]

# get all differentially methylated CpG Islands
myDiff_islands_25p <- getMethylDiff(myDiff_islands,difference=10,value=0.05)

myDiff_islands_25p_ann <- annotateWithGeneParts(as((myDiff_islands_25p), "GRanges"), refseq_anot)
# View the distance to the nearest Transcription Start Site; the target.row column indicates the row number in myDiff_islands_25p
head(getAssociationWithTSS(myDiff_islands_25p_ann))

# Reconstruct original object, keeping a lower coverage this time
myobj_lowCov <- methRead(files,
                         sample.id=list("control", "tumor"),
                         assembly="hg38", # replace with your genome assembly
                         treatment=c(0,1), # define control (0) and treatment (1) if applicable
                         pipeline="bismarkCoverage",
                         mincov = 3)

# Group the counts
tiles <- tileMethylCounts(myobj_lowCov,win.size=1000,step.size=1000,cov.bases = 10)

# Inspect data
head(tiles[[1]])

## export as bedgraph for visualisation
bedgraph(myDiff25p, col.name = "meth.diff", file.name = "results/figures/dmr/diff_cpg_25p.bed")

