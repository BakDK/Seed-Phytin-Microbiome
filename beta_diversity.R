# Beta diversity

library(phyloseq)
library(vegan)
library(tidyverse)
seed_phyl<-readRDS("RDS_files/phyloseq_final.rds")
asv_df<-data.frame((otu_table(seed_phyl)))
min(sample_sums(seed_phyl))

#Initial display of all samples in a PCA
full_bc_dist <-avgdist(asv_df,dmethod = "bray", sample= min(rowSums(asv_df))) #100 iterations is default         
               
library(ecodist)

pca_full<-pco(full_bc_dist)
# Determine how much the first two Principal components explain
PC_full<-sum(pca_full$values)
PC1_full<-pca_full$values[1]/PC_full*100 # 50.08%
PC2_full<-pca_full$values[2]/PC_full*100 # 29.38%
Full_bray_curtis_pcoa_df <- data.frame(pcoa1 = pca_full$vectors[,1], 
                                     pcoa2 = pca_full$vectors[,2])

Full_bray_curtis_pcoa_df$Sample_ID<-rownames(pca_full$vectors)

# Create a plot
inner_join(Full_bray_curtis_pcoa_df,sample_data(seed_phyl), by = c("Sample_ID" ="Sample_ID"))%>%
  ggplot(aes(x=pcoa1, y=pcoa2)) + geom_point(aes(shape = Timepoint, color = Treatment))+
  labs(x = "PCoA 1 50%" ,
       y = "PCoA 2 29%" ) 


#remove soil samples from T0
asv_df_plant<- asv_df %>% filter(!rownames(asv_df) %in% c("Soil1_T0","Soil2_T0","Soil3_T0")) 
plant_bc_dist <-avgdist(asv_df_plant,dmethod = "bray", sample= min(rowSums(asv_df_plant))) #100 iterations is default

plant_meta<-as_tibble(sample_data(seed_phyl)) %>% filter(!Sample_ID %in% c("Soil1_T0","Soil2_T0","Soil3_T0"))

# ======================= PCoA ===================

plot_pcoa_fb = function(data,data2,Var1, Var2){
  print("Computing PCoA")
  list_output = list()
    pcoa = pco(data2) #make ordination
  SUM_PCs = sum(pcoa$values)
  PC1 = pcoa$values[1]/SUM_PCs*100
  PC2 = pcoa$values[2]/SUM_PCs*100
  BC_df = data.frame(pcoa1 = pcoa$vectors[,1],
                     pcoa2 = pcoa$vectors[,2])
  BC_df$Sample_ID<-rownames(pcoa$vectors)
  list_output[[1]]=BC_df
  list_output[[2]]=pcoa
  list_output[[3]]= PC1
  

  print("Plotting")
 P1= inner_join(data,BC_df, by = c("Sample_ID"="Sample_ID"))%>%
    ggplot(aes(x = pcoa1, y = pcoa2))+geom_point(aes(shape = !!sym(Var1), color =!!sym(Var2)))+
    labs(x = paste(round(PC1,1),"%", sep =""), y = paste(round(PC2,1),"%",sep = ""))
  list_output[[4]] = P1
  return(list_output)
  
}

p4<-plot_pcoa_fb(plant_meta,plant_bc_dist, "Timepoint","Treatment")                            
p4[[1]]
p4[[4]]

# ================ Seed samples only ============


# ========= PERMANOVA ========
komsa<-data.frame(otu_table(plant_phy))
plant_phy<-subset_samples(seed_phyl, !Sample_ID %in% c("Soil1_T0","Soil2_T0","Soil3_T0"))
# First, rarefy to even depth
rarefied_physeq_list = list() # a list to store the rarefied phyloseq objects
for (i in 1:100) {
  rarefied_physeq = rarefy_even_depth(plant_phy, sample.size =  min(sample_sums(plant_phy)), rngseed = i)
  rarefied_physeq_list[[i]] <- rarefied_physeq
  cat("Loop number: \r", i, "out of 100")
}
otu_tables = lapply(rarefied_physeq_list, function(x) otu_table(x)) # regroup all the otu_tables in a list

cat("All the rarefied phyloseq objects have been created")
cat("Converting otu_tables to dataframes...")

for (i in 1:length(otu_tables)) { # convert all the otu_tables to dataframes
  otu_tables[[i]] = as.data.frame(otu_tables[[i]])
  otu_tables[[i]]$ASV = rownames(otu_tables[[i]])
  cat("Loop number: \r", i, "out of 100")
}

cat("Dataframes created")
cat("Concatenating dataframes...")

df_rare = as.data.frame(otu_tables[[1]]) # create a dataframe with the first otu_table for later
for (i in 2:length(otu_tables)) { # Concatenate all the dataframes vertically
  df_rare = bind_rows(df_rare, as.data.frame(otu_tables[[i]]))
  cat("Loop number: \r", i, "out of 100")
}

cat("Rarefied dataframe created")
cat("Calculating mean and median...")

# Regroup rows by NAME and calculate the mean and median
df_rare_mean = df_rare %>% group_by(ASV) %>% summarise(across(everything(), list(mean = mean)))
df_rare_median = df_rare %>% group_by(ASV) %>% summarise(across(everything(), list(median = median)))

names(df_rare_mean) = sub("_mean", "", names(df_rare_mean)) # Remove mean and median to column names
names(df_rare_median) = sub("_median", "", names(df_rare_median))

df_rare_mean <- df_rare_mean %>%
  tibble::column_to_rownames("ASV")

df_rare_median <- df_rare_median %>%
  tibble::column_to_rownames("ASV")

df_rare_mean = otu_table(df_rare_mean, taxa_are_rows = FALSE)
df_rare_median = otu_table(df_rare_median, taxa_are_rows = FALSE)

TAX<-tax_table(plant_phy)
META_Sam<-sample_data(plant_phy)
Rarefied_physeq_mean = phyloseq(df_rare_mean, TAX, META_Sam) # create another phyloseq object with rarefied data as OTU table
Rarefied_physeq_median = phyloseq(df_rare_median, TAX, META_Sam)

otu_data = round(otu_table(Rarefied_physeq_mean))
# In case where we have float number in mean (for richness), we round the number to have an integer
Rarefied_physeq_mean_rounded <- phyloseq(otu_data, TAX, META_Sam)

saveRDS(Rarefied_physeq_mean, "RDS_files/physeq_rarefied.rds")
saveRDS(Rarefied_physeq_median, "RDS_files/physeq_rarefied_median.rds")
saveRDS(Rarefied_physeq_mean_rounded, "RDS_files/physeq_rarefied_mean_rounded.rds")


#==========================PERMANOVA=========================
Rarefied_physeq_mean_rounded<-readRDS("RDS_files/physeq_rarefied_mean_rounded.rds")
asv_plant_df<-data.frame(otu_table(plant_phy))
plant_bc_dist <-avgdist(asv_plant_df,dmethod = "bray", sample= min(rowSums(asv_plant_df))) #100 iterations is default

set.seed(1)
adonis2(plant_bc_dist~Treatment*Timepoint, data = data.frame(sample_data(plant_phy)), permutations = 999)

# Test of seed samples only

fro_phyl<-subset_samples(Rarefied_physeq_mean_rounded, Compartment %in% "Seed")
fro_phyl2<-prune_taxa(taxa_sums(fro_phyl)>0,fro_phyl)
asv_fro_df<-data.frame(otu_table(fro_phyl2))
fro_bc_idst<-avgdist(asv_fro_df, dmethod = "bray", sample= min(rowSums(asv_fro_df))) #27,291

set.seed(1)
adonis2(fro_bc_idst~Treatment*Timepoint, data = data.frame(sample_data(fro_phyl)), by = "terms",permutations = 999)

fro_meta<-as_tibble(sample_data(fro_phyl))

p3<-plot_pcoa_fb(fro_meta,fro_bc_idst, "Timepoint","Treatment")                            
x11(width = 5, height = 5)
seed_pcoa_plot_f<-p3[[4]]+theme_bw() +theme(panel.grid = element_blank(), strip.background = element_blank())+
  scale_color_manual(values= c("#111184","darkgrey"), guide = "none")+theme(legend.position = "bottom")+
  scale_shape_manual(values = c(15,16,17), labels = c("0 DAS","2 DAS", "6 DAS"))

ggsave("Test_plots/PCoA_seedonly.png",seed_pcoa_plot_f)
dev.off()
#Import from "alpha_div_rarefaction_curves.R"
Shannon_plot2<-readRDS("RDS_files/Shannon_plot.rds")

library(ggpubr)

x11(width = 7, height = 4)
ggarrange(Shannon_plot2,seed_pcoa_plot_f, labels = c("A","B"))
ggsave("Test_plots/SHannon_and_PCOA.png",dpi = 300)
dev.off()
