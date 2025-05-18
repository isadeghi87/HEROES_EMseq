setwd('/home/i439h/projects/EMseq/codes/nextflow/work/')

#library
library(ggplot2)
library(ggpubr)
library(tidyr)

f = list.files(path =getwd(), pattern = 'bismark_summary_report.txt',full.names = T,recursive = T,all.files = T)
f = f[grep('csf',f,invert = T)]
## the other qc files from Temp folder
dir2 = '/home/i439h/projects/Emseq_temp/codes/nextflow/work/'
f2= list.files(path =dir2, pattern = 'bismark_summary_report.txt',full.names = T,recursive = T,all.files = T)
dat = data.frame()

## the other qc files from Temp folder
dir2 = '/home/i439h/projects/Emseq_temp/codes/nextflow/work/'
f2= list.files(path =dir2, pattern = 'bismark_summary_report.txt',full.names = T,recursive = T,all.files = T)
allfiles = c(f2,f)

qcdata = data.frame()
for( file in allfiles){
  d = read.delim(file)
  print(d)
  qcdata = rbind(qcdata,d)
}

colnames(qcdata) = gsub("File",'sample',colnames(qcdata))
qcdata$sample = gsub('OE0290-PED_','',qcdata$sample)
qcdata$sample = gsub('_1_val_1_bismark_bt2_pe.bam','',qcdata$sample)
qcdata = qcdata[!duplicated(qcdata),]
qcdata = qcdata[!duplicated(qcdata$sample),]
exclude = paste0(c('I034-044','I034-034'),collapse = "|")
qc_filtered = qcdata[grep(exclude,qcdata$sample,invert = T),]

# Add a new column for the rate of methylated CpGs
qc_filtered$Methylated_CpG_Rate <- qc_filtered$Methylated.CpGs / (qc_filtered$Methylated.CpGs + qc_filtered$Unmethylated.CpGs)*100

n = nrow(qc_filtered)
#qc_filtered$sample = paste0('sample',1:n)
qc_filtered$sample = factor(qc_filtered$sample,levels = qc_filtered$sample)

#Plot 1: Bar plot of Total Reads
p1 = ggplot(qc_filtered, aes(x = Total.Reads,y = sample)) +
  geom_bar(stat = "identity", fill = "blue",color = 'darkgrey',width = 0.7) +
  scale_x_log10()+
  labs(title = "Total Reads per sample", y = "sample", x = "Total Reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw() # Removes the x-axis text

#Plot 2: Bar plot of aligned Reads
p2 = ggplot(qc_filtered, aes(x = Aligned.Reads,y = sample)) +
  geom_bar(stat = "identity", fill = "skyblue",color = 'darkgrey',width = 0.7) +
  scale_x_log10()+
  labs(title = "Aligned Reads per sample", y = "", x = "Aligned Reads") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw()+
  theme(axis.text.y = element_blank()) 

##plot3 methylated cpg rate barplot
p_meth_perc = ggplot(qc_filtered, aes(x = Methylated_CpG_Rate,y = sample)) +
  geom_bar(stat = 'identity', fill = "green", color = "darkgrey",width = 0.7) +
  labs(title = "Percentage of Methylated CpGs", y = "sample", x = "% methylation") +
  scale_x_continuous(limits = c(0,100))+
  theme_bw()


##plot4 number of methylated cpgs
p_meth = ggplot(qc_filtered, aes(x = Methylated.CpGs,y = sample)) +
  geom_bar(stat = 'identity',fill = "orange", color = "darkgrey",width = 0.7) +
  labs(title = "Methylated CpG sites", y = "", x = "methylated") +
  theme_bw()+
  theme(axis.text.y = element_blank())

##plot5 number of unmethylated cpgs
p5 = ggplot(qc_filtered, aes(x = Unmethylated.CpGs,y = sample)) +
  geom_bar(stat = 'identity',fill = "grey", color = "black",width = 0.7) +
  labs(title = "Unmethylated CpGs sites", y = "", x = "unmethylated") +
  theme_bw()+
  theme(axis.text.y = element_blank())

## coverage ####
cov1 = list.files(path ='/home/i439h/projects/EMseq/codes/nextflow/work/', pattern = 'general_stats.txt',full.names = T,recursive = T,all.files = T)
cov2 = list.files(path ='/home/i439h/projects/Emseq_temp/codes/nextflow/work/', pattern = 'general_stats.txt',full.names = T,recursive = T,all.files = T)
cov_files = c(cov2,cov1)

cov_dat =list()
for( file in cov_files){
  d = read.delim(file)
  if(ncol(d)==27){
  cov_dat = rbind(cov_dat,d)
    
  }
}
cov_dat$Sample = gsub('OE0290-PED_','',cov_dat$Sample)
colnames(cov_dat) = gsub('Bismark_mqc.generalstats.bismark.','',colnames(cov_dat))
colnames(cov_dat) = gsub('QualiMap_mqc.generalstats.qualimap.','',colnames(cov_dat))
colnames(cov_dat) = gsub('FastQC_mqc.generalstats.fastqc.','',colnames(cov_dat))
colnames(cov_dat) = gsub('Cutadapt_mqc.generalstats.cutadapt.','',colnames(cov_dat))
cov_dat = cov_dat[grep(exclude,cov_dat$sample,invert = T),]

mean_cov = cov_dat %>% dplyr::select(Sample,mean_coverage)
mean_cov = na.omit(mean_cov)
mean_cov = mean_cov[!duplicated(mean_cov),]
samp_names = qc_filtered$sample
mean_cov = mean_cov[grep(paste0(samp_names,collapse = '|'),mean_cov$Sample),]
id_match = match(qc_filtered$sample,mean_cov$Sample)
mean_cov = mean_cov[id_match,]
all(qc_filtered$sample==mean_cov$Sample)

#mean_cov$sample  =  paste0('sample',1:n)
mean_cov$sample = factor(mean_cov$Sample,levels = unique(mean_cov$Sample))

##plot6 coverage
p6 = ggplot(mean_cov, aes(x = mean_coverage,y = sample)) +
  geom_bar(stat = 'identity',fill = 'red', color = "black",width = 0.7) +
  labs(title = "Mean coverage per sample", y = "", x = "Mean coverage (x)") +
  theme_bw()

## distribution of x_pc
x_cov = cov_dat %>% dplyr::select(sample,contains('x_pc'))
x_cov= na.omit(x_cov)
x_cov= x_cov[!duplicated(x_cov),]
samp_names = qc_filtered$sample
x_cov= x_cov[grep(paste0(samp_names,collapse = '|'),x_cov$sample),]
id_match = match(qc_filtered$sample,x_cov$sample)
x_cov= x_cov[id_match,]
all(qc_filtered$sample==x_cov$Sample)

cov_long = pivot_longer(x_cov,cols = 2:6,names_to = 'Coverage',values_to = 'Percentage')
covs = unique(cov_long$Coverage)
cov_long$Coverage = factor(cov_long$Coverage,levels = covs)

# Generate Line Plot
p_line = ggplot(cov_long, aes(x = Coverage, y = Percentage, group = sample, color = sample)) +
  geom_line(size = 1,show.legend = T) +
  geom_point(size = 2,show.legend = F) +
  theme_minimal() +
  labs(title = "Coverage Trend Across Samples", x = "Coverage Depth", y = "Percentage (%)")


## total c ####
tot_c = cov_dat %>% dplyr::select(sample,total_c)
tot_c = na.omit(tot_c)

tot_c$sample = gsub("_1_val_1","",tot_c$sample)
tot_c= tot_c[!duplicated(tot_c),]
samp_names = unique(qc_filtered$sample)
tot_c= tot_c[grep(paste0(samp_names,collapse = '|'),tot_c$sample),]
id_match = match(qc_filtered$sample,tot_c$sample)
tot_c= tot_c[id_match,]
all(qc_filtered$sample==tot_c$Sample)

#tot_c$sample  =  paste0('sample',1:n)
tot_c$sample = factor(tot_c$sample,levels = tot_c$sample)

p_totalC = ggplot(tot_c, aes(x = total_c,y = sample)) +
  geom_bar(stat = 'identity',color = "black",width = 0.7,fill = 'orange') +
  labs(title = "Total Cs per sample", y = "", x = "Total Cs") +
  theme_bw()+
  theme(axis.text.y = element_blank())


### C to T coversions ####
# Load necessary library

# Step 1: Read the text file
setwd('~/Documents/qc_files/splitting_report/')
cfiles = list.files(path ='/home/i439h/projects/Emseq_temp/codes/nextflow/work/', pattern = 'deduplicated_splitting_report.txt',full.names = T,recursive = T,all.files = T)
cfiles = cfiles[grep('csf',cfiles,invert = T)]
cfiles2 = list.files(path ='/home/i439h/projects/EMseq/codes/nextflow/work/', pattern = 'deduplicated_splitting_report.txt',full.names = T,recursive = T,all.files = T)
allcfiles = c(cfiles, cfiles2)

samp_names = qc_filtered$sample
cfiles= cfiles[grep(paste0(samp_names,collapse = '|'),cfiles)]

cfdat = data.frame()
for( file in allcfiles){
  # Read the file content into a character vector
  text_data <- readr::read_lines(file)
  print(text_data)
  if(length(text_data)>10){
  # Step 2: Extract relevant lines using regular expressions
  cpg_line <- grep("Total C to T conversions in CpG context:", text_data, value = TRUE)
  c_to_t_cpg <- as.numeric(gsub("Total C to T conversions in CpG context:\\s*([0-9]+)", "\\1", cpg_line))
  
  chg_line <- grep("Total C to T conversions in CHG context:", text_data, value = TRUE)
  c_to_t_chg <- as.numeric(gsub("Total C to T conversions in CHG context:\\s*([0-9]+)", "\\1", chg_line))
  
  chh_line <- grep("Total C to T conversions in CHH context:", text_data, value = TRUE)
  c_to_t_chh <- as.numeric(gsub("Total C to T conversions in CHH context:\\s*([0-9]+)", "\\1", chh_line))
  
  #sample names
  sample_id <- gsub(".*/OE0290-PED_(.*)_1_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt", "\\1", file)
  
  # Step 3: Prepare the data for plotting
  conversion_data <- data.frame(
    sample = sample_id,
    Context = c("CpG", "CHG", "CHH"),
    Conversions = c(c_to_t_cpg, c_to_t_chg, c_to_t_chh)
  )

  cfdat = rbind(cfdat,conversion_data)
  }
  }

# Flag duplicates based on 'sample' and 'Conversions'
duplicate_flags <- !duplicated(cfdat[, c("sample", "Conversions")])
cfdat = cfdat[duplicate_flags,]

exclude = paste0(c('I034-044','I034-034'),collapse = "|")
cf_filtered = cfdat[grep(exclude,cfdat$sample,invert = T),]

# Stacked bar plot to visualize proportion of conversions for different contexts within each sample
ggplot(cf_filtered, aes(x = sample, y = Conversions, fill = Context)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Proportion of C to T Conversions by Context for Each Sample",
       x = "Sample",
       y = "Total C to T Conversions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Bar plot to compare conversions across different contexts for each sample
ggplot(cf_filtered, aes(x = sample, y = Conversions, fill = Context)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Total C to T Conversions by Context for Each Sample",
       x = "Sample",
       y = "Total C to T Conversions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Step 1: Calculate the average conversions for each context
avg_conversions <-cf_filtered %>%
  group_by(Context) %>%
  summarise(Average_Conversions = mean(Conversions))

# Print the averages
print(avg_conversions)

# Step 2: Create a box plot for the conversions by context
p_conv = ggplot(cf_filtered, aes(x = Context, y = Conversions, fill = Context)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16) +
  geom_point()+
  theme_minimal() +
  scale_y_log10()+
  labs(title = "Total C to T Conversions by Context",
       x = "Methylation Context",
       y = "C to T Conversions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_conv)

### merged plots

p_reads = ggarrange(p1,p_line,labels = c('A','B'))
p_meth = ggarrange(p_meth_perc, p_conv,labels = c('C','D'))

fig1 = p_reads / p_meth
setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/")
ggsave(filename = "./EMseq/results/figures/report/Figure1_QC.png",device = 'png',plot = fig1,width = 10,height = 8)
