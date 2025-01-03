# Plotting heatmap (Figure 6)
library(ampvis2)
library(tidyverse)

amp_ob<-readRDS("RDS_files/ampvis_obj.rds")

amp_ob$metadata$Days[amp_ob$metadata$Days %in% "14"]<-21
amp_ob$metadata$Treatment[amp_ob$metadata$Treatment %in% "NS"]<-"Control"
amp_ob$metadata$Treatment[amp_ob$metadata$Treatment %in% "Sterile"]<-"Sterilized"

seed_heat_v1<-
amp_ob %>% amp_subset_samples(Compartment %in% "Seed") %>%
  amp_heatmap(group_by = "Treatment",
              facet_by = "Days",
              tax_aggregate = "Genus",
              tax_show = 15,
              
              round = 1,
              color_vector = c("wheat","darkgreen"))+ 
  theme(strip.background = element_rect(fill = "white"))


PEP_subset<-amp_ob %>% 
  amp_subset_taxa(tax_vector = c("Erwinia","Pseudomonas","Pantoea"), normalise = TRUE)

PEP_data<-PEP_subset %>%amp_heatmap(group_by = "Sample_ID",
                        tax_aggregate = "Genus",
            round = 1,
            textmap = TRUE, normalise = FALSE) %>% t() %>% cbind(.,Sample_ID=rownames(.)) %>%
  merge(.,PEP_subset$metadata, by = "Sample_ID")

#Turn into long format
PEP_data_long<-PEP_data %>% pivot_longer(cols= Pseudomonas:Erwinia, names_to = "Genus",values_to = "Abundance")

# plot timely development of the 3 genera
PEP_data_long$Abundance<-as.numeric(PEP_data_long$Abundance)

PEP_data_long$Days[PEP_data_long$Days %in% "14"]<-21
raw_pep_plot<-ggplot(PEP_data_long, aes(x=Days, y = log10(Abundance), color = Treatment)) + geom_point() + 
  facet_wrap(~Genus,  scales = "free_y")

new_PEP_plot<-raw_pep_plot + theme_bw()+
  theme(strip.text = element_text(size = 11,face = "italic"),
        panel.grid = element_blank(), strip.background = element_blank(),
        panel.background = element_blank()) + labs(y = "Log10 Relative Abundance", x = "Days after sowing (DAS)")+ 
  scale_color_manual(values=c("#111184","#90EE90","darkgrey"),
                     labels = c("Control","Soil", "Sterilized"))+
  scale_x_continuous(breaks = c(0,2,6,21))

#Combine with the heatmap
library(ggpubr)
x11(width = 12, height = 7)
ggarrange(seed_heat_v1,new_PEP_plot,ncol = 1,labels = c("A","B"))
ggsave("Test_plots/Heatmap_and_PEP.png", dpi = 300) 

#Which Pseudomonas ASVs are found?
PEP_subset %>%
  amp_subset_taxa(tax_vector = "Pseudomonas", normalise = FALSE) %>%
  amp_heatmap(group_by = "Treatment",
              facet_by = "Days",
              tax_aggregate = "OTU",
              tax_add = "Genus",
              round = 2,
              tax_show = 15, normalise = FALSE)
#Need to follow up on this (20-12-2024)              


## Supplementary heat map of roots only

heat_root<-amp_ob %>% amp_subset_samples(Compartment %in% "Root") %>%
  amp_heatmap(group_by = "Treatment",
              facet_by = "Treatment",
              tax_aggregate = "Genus",
              tax_show = 20,
              tax_add = "Family",
              round = 1) +theme(axis.text.x = element_blank(),
                                axis.ticks.x = element_blank())
x11(width = 8, height = 7)
heat_root
ggsave("Test_plots/Heatmap_roots_sup.png", dpi = 300)
dev.off()
soil_heat_v1 <-amp_ob %>% amp_subset_samples(Compartment %in% "Soil") %>%
  amp_subset_taxa(tax_vector = "Pantoea",normalise = TRUE) %>%
  amp_heatmap(tax_aggregate = "Genus",normalise = FALSE, round = 3)
