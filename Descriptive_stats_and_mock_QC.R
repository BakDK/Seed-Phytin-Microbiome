# Quality checking the data
library(phyloseq)
library(tidyverse)
library(xlsx)
seed_phy<-readRDS("Syncom_phyloseq.rds")
sample_names(seed_phy)
#There are a lot of duplicates, due to the way Eurofins sent the data
# I will remove duplicate samples from the 8373 run. Mainly as a lot of these samples have very few reads, 

#Modify sample names
sample_names(seed_phy)<-gsub("_1","",sample_names(seed_phy))
sample_names(seed_phy)<-gsub("lib\\d+_","", sample_names(seed_phy))
sample_names(seed_phy)[sample_names(seed_phy)=="R7_T1_S_8373"]<-"R7S_T1_8373"
sample_names(seed_phy)[sample_names(seed_phy)=="R8_T1_S_8373"]<-"R8S_T1_8373"
sample_names(seed_phy)[sample_names(seed_phy)=="R8_T1_S_8361"]<-"R8S_T1_8361"
sample_names(seed_phy)[sample_names(seed_phy)=="Soil_2_T0_8373"]<-"Soil2_T0_8373"
sample_names(seed_phy)[sample_names(seed_phy)=="Soil_3_T0_8373"]<-"Soil3_T0_8373"
sample_names(seed_phy)[sample_names(seed_phy)=="R4T1_lib7223010304.fastq.gz"]<-"R4_T1_8222" #THe 8222 is arbitrary number to give the same format as other sample names
sample_names(seed_phy)[sample_names(seed_phy)=="Soil1T0_lib7223020304.fastq.gz"]<-"Soil1_T0_8222" #THe 8222 is arbitrary number to give the same format as other sample names
sample_names(seed_phy)[sample_names(seed_phy)=="Zymo_lib7223030304.fastq.gz"]<-"Mock1_T0_8222" #THe 8222 is arbitrary number to give the same format as other sample names

df<-str_split_fixed(sample_names(seed_phy),"_",3) %>% .[,-5] %>% as.data.frame()

rownames(df)<-sample_names(seed_phy)
colnames(df)<-c("Rep","Time","Run")
df$Sample_ID_long<-sample_names(seed_phy)
duplicated(df[,1:2])

df_red<-df %>% 
  group_by(Rep, Time) %>%
  mutate(num_dupli=n(),
         dup_id = Run) %>%
  ungroup() %>%
  filter(num_dupli == 1 |num_dupli >1 & !dup_id == 8373)

df_red$Sample_ID<-gsub('_[0-9]+',"",df_red$Sample_ID_long)
df_red$Sample_ID[df_red$Sample_ID %in% "Control"]<-"Control_4"
df_red$Sample_ID[df_red$Sample_ID %in% "Mock1_T0"]<-"Mock_T0"
#Now the phyloseq object is reduced based on this data frame

sample_names(seed_phy)
provdat<-sample_data(data.frame(Sample_ID_long=sample_names(seed_phy)))
sample_names(provdat)<-sample_names(seed_phy)
seed_phy2<-merge_phyloseq(seed_phy,provdat)
phyl_fin<-subset_samples( seed_phy2,Sample_ID_long %in% df_red$Sample_ID_long  )
sample_names(phyl_fin)<-df_red$Sample_ID
#Import metadata sheet
metadata<-read.xlsx("Metadata.xlsx", sheetIndex = 1)

metadata<-sample_data(metadata)
sample_names(metadata)<-metadata$Sample_ID
setdiff(sample_names(phyl_fin),metadata$Sample_ID)
setdiff(metadata$Sample_ID, sample_names(phyl_fin)) #This is correct. 8 soil samples, 2 T0 samples and three negative controls are not sequenced
sample_data(phyl_fin)<-metadata

phyl_fin0<-phyl_fin
phyl_fin<-prune_taxa(taxa_sums(phyl_fin)>0,phyl_fin )
# Total of 45 samples

phyl_fin<-subset_taxa(phyl_fin, Family != "Mitochondria") #Removes 23 ASV
phyl_fin<-subset_taxa(phyl_fin, Class != "Chloroplast") #Removes 0 ASV
phyl_fin<-subset_taxa(phyl_fin, Kingdom = "Bacteria") #Removes 0 ASV

#This gives 2463 ASVs


#Descriptive statistics of data and check of mock community

library(ampvis2)
phyl_ob<-phyl_fin

amp_ob<-amp_load(phyl_ob)
amp_ob$metadata

#Check control and mock samples
amp_ob %>% amp_subset_samples(Compartment %in% c("Control", "Mock")) %>%
  amp_heatmap(tax_aggregate = "Genus",
             # tax_add = c("Genus","Species"),
              tax_show = 25,
             round = 3)
# 92 ASVs present
#The theoretical values are:
#Lacobacillus fermentum: 18.4%
#Bacillus subtilis 17.4%
#S. aureus 15.5%
#Listeria monocytogenes 14.1%
#Salmonella enterica 10.4%
#E. coli 10.1%
#Enterococcus faecalis 9.9%
#P. aeruginosa 4.2%

#The sample primarily contains those 8 species. Other ASVs comprise <0.02% of reads. So some bias in the library prep
#But no worrying contamination. And only for one of the two mock samples (they are from two different runs).

amp_data<-amp_ob %>% amp_subset_samples(!Compartment %in% c("Control", "Mock"))
amp_data<-amp_data %>% amp_subset_samples(!(Compartment %in% "Soil" & Timepoint %in% "T1"))
# 41 samples, 2336 ASVs
amp_data
#Total reads: 1,833,595, min reads 19,714, max reads 87,676, median reads 42,835

amp_data$metadata$Treatment[is.na(amp_data$metadata$Treatment)]<-"Soil"
saveRDS(amp_data,"RDS_files/ampvis_obj.rds")

# Save phyloseq without these files
phyl_ob<-subset_samples(phyl_ob,!Compartment %in% c("Control", "Mock"))
phyl_ob<-subset_samples(phyl_ob,!(Compartment %in% "Soil" & Timepoint %in% "T1"))
phyl_fin<-prune_taxa(taxa_sums(phyl_ob)>0,phyl_ob )

saveRDS(phyl_fin,"RDS_files/phyloseq_final.rds")