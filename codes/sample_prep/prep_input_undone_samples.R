setwd("/home/i439h/projects/Emseq_temp/")

## read list of LB samples
lb_df = read.csv("./datasets/HEROES_AYA_LIQUID_BIOPSY.csv",skip = 1)
samples = na.omit(lb_df$`Liquid Biopsy ID`)
samples <- sub("-P.*", "", samples)

## get list of processed samples

dir = "/home/i439h/projects/EMseq/results/nextflow/bismark/methylation_calls/methylation_coverage/"
files = list.files(path = dir,all.files = T,full.names = F,recursive = T)
files2 =list.files(path = "./results/nextflow/bismark/methylation_calls/methylation_coverage/"
                   ,all.files = T,full.names = F,recursive = T)

allfiles = c(files,files2)

# Extract the desired ID parts using a regular expression
ids <- sub(".*_(\\w+-\\d{3}).*", "\\1", allfiles)

## get those not processed
not_done = setdiff(samples,ids)


## list fastq files
dir = '/omics/odcf/project/OE0290/pediatric_tumor/sequencing/whole_genome_bisulfite_sequencing/view-by-pid'
fastq = list.files(path = dir, pattern = '*.fastq.gz',full.names = T,recursive = T)
fastq.filt  = fastq[grep('plasma',fastq)] 

r1 = fastq.filt[grep('R1',fastq.filt)]
r2 = fastq.filt[grep('R2',fastq.filt)]
patient = str_split(r1,pattern = '/',n=12,simplify=T)[,10]
sample = str_split(r2,pattern = '/',n=12,simplify=T)[,11]
patient_sample = paste0(patient,'_',sample)
df= data.frame(sample=patient_sample,
               fastq_1=r1,
               fastq_2=r2)

new_ids = c('2LB-077_plasma-01-01',
            '2LB-197_plasma-03-01',
            '2LB-257_plasma-01-01',
            '2LB-269_plasma-01-01',
            '2LB-269_plasma-02-01',
            '2LB-304_plasma-01-01',
            '5LB-003_plasma-01-01',
            '5LB-017_plasma-01-01')
df_filt  = df[grep(paste0(new_ids,collapse = "|"),df$sample),] 

## Heroes samples 
hero_ids = c('H021-E7U1N6_plasma-03-01',
             'H021-PJCDQK_plasma-03-01')
dir = '/omics/odcf/project/OE0290/heroes-aya/sequencing/whole_genome_bisulfite_sequencing/view-by-pid/'
fastq2 = list.files(path = dir, pattern = '*.fastq.gz',full.names = T,recursive = T)
fastq2 = fastq2[grep('plasma',fastq2)]

r1 = fastq2[grep('R1',fastq2)]
r2 = fastq2[grep('R2',fastq2)]
patient = str_split(r1,pattern = '/',n=12,simplify=T)[,11]
sample = str_split(r2,pattern = '/',n=13,simplify=T)[,12]
patient_sample = paste0(patient,'_',sample)
df2= data.frame(sample=patient_sample,
               fastq_1=r1,
               fastq_2=r2)

allsamples = rbind(df_filt,df2)
allsamples$sample = gsub('OE0290-PED_','',allsamples$sample)
allsamples$sample = gsub('OE0290_HEROES-AYA_','',allsamples$sample)

readr::write_delim(allsamples,file='./datasets/sample_sheets/sampleSheet_undone.csv',delim = ',',col_names = T)
