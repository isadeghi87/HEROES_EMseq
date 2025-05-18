# drectory of samples 
setwd('/omics/odcf/project/OE0290/pediatric_tumor/sequencing/whole_genome_bisulfite_sequencing/view-by-pid')
setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/EMseq/data/raw_data/')
library(stringr)


# read metadata
metadata = read.csv('/home/i439h/projects/WGBS/data/Alignment_Quality_Control-OE0290_pediatric_tumor_2023-09-25.csv')
patient_list = list.dirs(recursive = F)
patient_list = gsub("./",'',patient_list)

## list fastq files
fastq = list.files(pattern = '*.fastq.gz',recursive = T)
fastq.filt  = fastq[grep('plasma',fastq)] 

df_sample = data.frame(patient,sample)

fastqs = paste0('/omics/odcf/project/OE0290/pediatric_tumor/sequencing/whole_genome_bisulfite_sequencing/view-by-pid/',fastq.filt)
r1 = fastqs[grep('R1',fastqs)]
r2 = fastqs[grep('R2',fastqs)]
patient = str_split(r1,pattern = '/',n=12,simplify=T)[,10]
sample = str_split(r2,pattern = '/',n=12,simplify=T)[,11]
patient_sample = paste0(patient,'_',sample)
df= data.frame(sample=patient_sample,
               fastq_1=r1,
               fastq_2=r2)
exclude = paste0(c('I034-044','I034-034'),collapse = "|")
df = df[grep(exclude,df$sample,invert=T),]

## read only plasma data 
# plasma  = readr::read_delim('/home/i439h/projects/EMseq/data/HEROES-EMseq-plasma.csv')
# id =  paste0(plasma$`Patient ID`,'_',plasma$`Sample Type`)
# df_plasma = df[grep(paste0(plasma$`Patient ID`,collapse = '|'),df$sample),  ]

# condition= read.csv('/home/i439h/projects/WGBS/data/Sequences_of_samples-OE0290_pediatric_tumor_2023-10-12.csv')
# condition = subset(condition,X %in% c( "Healthy Controls", "Tumor"))
# id = paste0(condition$Patient.ID,'_',condition$Sample.Type)

readr::write_delim(df,file='../sampleSheet.csv',delim = ',',col_names = T)

## new samples
new_ids = c('OE0290-PED_2LB-189','OE0290-PED_2LB-217','OE0290-PED_5LB-049')
samsheet = df
samsheet = samsheet[grep(paste(new_ids,collapse = '|'),samsheet$sample),]

readr::write_delim(samsheet,file='../sampleSheet_new_ids.csv',delim = ',',col_names = T)
  

