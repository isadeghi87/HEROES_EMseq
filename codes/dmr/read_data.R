

setwd('/Users/i439h/Library/Application Support/Mountain Duck/Volumes.noindex/home.localized/projects/Emseq_temp/')

# Main analysis package
library("methylKit")
# Annotation package
library("genomation")
library("GenomicRanges")
library(bsseq)
library(magrittr)
library(readr)
library(dmrseq)
options(future.globals.maxSize = 48 * 1024^3)

## metdata
# Define the list containing the bismark coverage files.

files1 <- list.files(path ='../EMseq/results/nextflow/bismark/methylation_calls/methylation_coverage',full.names = T)
exclude = paste0(c('I034-044','I034-034','2LB-128'),collapse = "|")
files1 = files1[grep(exclude,files1,invert = T)]

files2 <- list.files(path ='./results/nextflow/bismark/methylation_calls/methylation_coverage',full.names = T)

allfiles = c(files1,files2)
allfiles = allfiles[-11]
samples = stringr::str_extract(allfiles, "[A-Za-z0-9-]+_plasma-\\d{2}-\\d{2}")

annot = read.delim("~/Documents/emseq_annot.tsv")
## match diagnosis for samples
methyldata = data.frame(samples = samples, files = allfiles)
methyldata$diagnosis =ifelse(grepl('0LB',methyldata$samples),'control','tumor')
methyldata$treatment = ifelse(methyldata$diagnosis=='control',0,1)
methyldata$entity = annot$entity[match(methyldata$samples,annot$sample)]

# ## get intersect of samples present
# id = intersect(samples,metadata$samples)
# methyldata = methyldata[methyldata$samples %in% id,]

## merge with clinics
df =  as.data.frame(table(methyldata$entity))
colnames(df) = c("Entity","number")

library(ggplot2)
(pp = ggplot(df, aes(x = reorder(Entity, -number), y = number, fill = Entity)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = number), vjust = -0.5, size = 3) +
  labs(x = NULL, y = NULL,title = 'Number of Samples per Entity') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ))

ggsave(filename = "./results/QC/demographic/samples_barplot.pdf",plot = pp,width = 8,height = 4,units = 'in')
  

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
                  mincov = 1)


# check number of samples
head(myobj)

# # ## filter 
# myobj.filt <- filterByCoverage(myobj,
#                                lo.count=2,
#                                lo.perc=NULL,
#                                hi.count=NULL,
#                                hi.perc=99.9)
# ## normalize 
# myobj.filt.norm <- normalizeCoverage(myobj, method = "median")
# meth <- unite(myobj, destrand=F)
# 
# ## keep only main chromosomes
#chrs = paste0('chr',c(1:22,'X',"Y"))
# 
# # Step 1: Extract data from the methylBase object
# meth_data <- getData(meth)
# 
# # Step 2: Subset the data based on chromosomes
# meth_data_filt <- meth_data[meth_data$chr %in% chrs, ]
# 
# # Step 3: Optionally, reconstruct a new methylBase object from the filtered data
# meth_filt <- new("methylBase", meth_data_filt, sample.ids = meth@sample.ids, assembly = meth@assembly,
#                  context = meth@context, resolution = meth@resolution, treatment = meth@treatment,
#                  destranded = meth@destranded, numCs.index= meth@numCs.index, numTs.index = meth@numTs.index
#                  )
# 
# # Now meth_filt contains only the rows with the specified chromosomes
# 
# save(meth_filt, methyldata,file = './data/processed/normalized_methyl.rdata')

## regional
tiles = tileMethylCounts(myobj,win.size=1000,step.size=1000,cov.bases = 10)
reg_filt <- normalizeCoverage(tiles, method = "median")
meth_reg <- unite(reg_filt, destrand=F)

#  Extract data from the methylBase object
meth_data <- getData(meth_reg)

# ## keep only main chromosomes
chrs = paste0('chr',c(1:22,'X',"Y"))
#  Subset the data based on chromosomes
meth_data_filt <- meth_data[meth_data$chr %in% chrs, ]

#  Optionally, reconstruct a new methylBase object from the filtered data
meth_filt_reg <- new("methylBase", meth_data_filt, sample.ids = meth_reg@sample.ids, assembly = meth_reg@assembly,
                 context = meth_reg@context, resolution = meth_reg@resolution,treatment = meth_reg@treatment,
                 destranded = meth_reg@destranded, numCs.index= meth_reg@numCs.index, numTs.index = meth_reg@numTs.index)


save(meth_filt_reg, methyldata,file = './data/processed/normalized_methyl_region.rdata')
