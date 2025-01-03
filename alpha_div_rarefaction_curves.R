# Calculate alpha diversity
library(phyloseq)
library(ggpubr)
seed_phyl<-readRDS("RDS_files/phyloseq_final.rds")

seed_phyl_with_seqs<-seed_phyl
taxa_names(seed_phyl_with_seqs)<-refseq(seed_phyl)

taxonomic_table<-tax_table(seed_phyl_with_seqs)
taxonomic_table<-as.data.frame((taxonomic_table))

taxonomic_table$ASV<-taxa_names(seed_phyl)
write.csv(taxonomic_table,"Seed_tax_table.csv")
# Export the ASV table
asv_table<-otu_table(seed_phyl)
write.csv(asv_table, "Seed_asv_table.csv")
sample_sums(seed_phyl)
min(sample_sums(seed_phyl)) #19,714
#Calculating Shannon index by rarefying 100 times 

# Create a matrix with all the samples as rows and 100 empty columns
shan_ra<-matrix(NA,nsamples(seed_phyl),100)

# Then rarefy once and then calculate  Shannon diversity, repeat 100 times. Do not set rngseed as the values will the be the same
for(i in 1:100){
  {Full16_rare<-rarefy_even_depth(seed_phyl, sample.size = min(sample_sums(seed_phyl)))
  resultat<-plot_richness(Full16_rare ,x = "Sample_ID", measures = "Shannon")
  shan_ra[,i]<-resultat$data$value}
}

print(shan_ra)

#Calculate mean of the 100 rarefied Shannon indices
shan_ra_me<-apply(shan_ra,MARGIN =  1, mean)

shan_mean_met<-cbind(sample_data(seed_phyl),Shannon =shan_ra_me)

#  Plot the data
shan_mean_met$Treatment[is.na(shan_mean_met$Treatment)]<-"Soil"
shan_mean_met$Days[shan_mean_met$Days %in% 14]<-21
#Create a subdata set for stats
shan_seed_mean<-shan_mean_met %>% filter(!Treatment %in% "Soil")

stats_shan<-compare_means(Shannon~Treatment, group.by = "Days", shan_seed_mean, method = "t.test")
#Display only the significant p-value

Shannon_plot<-
shan_mean_met %>% 
  ggplot(aes(x = Days, y = Shannon, color = Treatment))+
  geom_point(position =position_dodge(width = 1))+theme_bw()+
  labs(color = "Sample")+theme(legend.position = "bottom")+
  xlab("Days after sowing (DAS)")+
  theme(panel.grid = element_blank(), strip.background = element_blank())+
  guides(color = guide_legend(nrow = 1))+ scale_color_manual(values=c("#111184","#90EE90","darkgrey"),
                                                             labels = c("Control","Soil", "Sterilized"))+
   geom_text(x = 1.2, y = 4, label = paste( "p =" , stats_shan$p.adj[stats_shan$Days == 0]),
             size = 3, color = "black")
Shannon_plot

ggsave("Test_plots/Shannon.png")
saveRDS(Shannon_plot,"RDS_files/Shannon_plot.rds") # Use this to combine with beta_diversity plot

#=============================================================================
  #Rarefaction curves
library(ampvis2)
amp_ob<-readRDS("RDS_files/ampvis_obj.rds")
amp_ob$metadata


comb_rarecurves<-amp_ob %>% 
  amp_rarecurve(color_by = "Treatment",
                facet_by = "Timepoint") +
  labs(y = "Number of ASVs") + scale_color_manual(values = c("#E36C14","#1c4eaa","#666667"))

x11(width = 5, height = 4)
comb_rarecurves
ggsave("Test_plots/Rarefaction_curves.png")


