# import the phyloseq object created in the "Descriptive_stats_and_mock_QC.R"
# Use corncob to test differential abundant taxa

library(phyloseq)
library(magrittr)
library(corncob)
seed_phyl<-readRDS("RDS_files/phyloseq_final.rds")

# Remove the control samples, and remove the soil samples from T1.
# This gives 3 soil samples at T0, 8 T0 seeds, 10 T1 seeds, 10 T2 seeds and 10 T3 roots -> 41 samples. 
seed2<-seed_phylo %>% subset_samples(!Compartment %in% "Control")
seed<-seed2 %>% subset_samples(!(Compartment %in% "Soil" & Time %in% "T1"))
seed<-prune_taxa(taxa_sums(seed)>0,seed)

#test at To
seed_t0<-subset_samples(seed_phylo, Time %in% "T0" & Compartment %in% "Seed")  %>% prune_taxa(taxa_sums(.)>0,.) %>% tax_glom("Genus")
# 106 genera
sample_data(seed_t0)

Genus_t0_analysis <-differentialTest(formula = ~ Treatment,
                                     phi.formula = ~ Treatment,
                                     formula_null = ~ 1,
                                     phi.formula_null = ~ Treatment,
                                     data = seed_t0,
                                     test= "Wald",
                                     boot = FALSE,
                                     fdr_cutoff = 0.05)

plot(Genus_t0_analysis)
Genus_t0_analysis$significant_models
t0_plot<-plot(Genus_t0_analysis)

#test at T1
seed_t1<-subset_samples(seed_phylo, Time %in% "T1" & Compartment %in% "Seed")  %>% prune_taxa(taxa_sums(.)>0,.) %>% tax_glom("Genus")
# 226 genera
sample_data(seed_t1)
Genus_t1_analysis <-differentialTest(formula = ~ Treatment,
                                     phi.formula = ~ Treatment,
                                     formula_null = ~ 1,
                                     phi.formula_null = ~ Treatment,
                                     data = seed_t1,
                                     test= "Wald",
                                     boot = FALSE,
                                     fdr_cutoff = 0.05)

t1_plot<-plot(Genus_t1_analysis)
Genus_t1_analysis$significant_models
str(t1_plot)
t1_plot$data

#test at T2
seed_t2<-subset_samples(seed_phylo, Time %in% "T2" & Compartment %in% "Seed")  %>% prune_taxa(taxa_sums(.)>0,.) %>% tax_glom("Genus")
# 233 genera
seed_t2

Genus_t2_analysis <-differentialTest(formula = ~ Treatment,
                                     phi.formula = ~ Treatment,
                                     formula_null = ~ 1,
                                     phi.formula_null = ~ Treatment,
                                      test= "Wald",
                                     boot = FALSE,
                                     data = seed_t2,
                                     fdr_cutoff = 0.05)
plot(Genus_t2_analysis)
Genus_t2_analysis$significant_models
t2_plot<-plot(Genus_t2_analysis)


root_t3<-subset_samples(seed_phylo, Time %in% "T3" )  %>% prune_taxa(taxa_sums(.)>0,.) %>% tax_glom("Genus")
sample_data(root_t3)
Genus_t3_analysis <-differentialTest(formula = ~ Treatment,
                                     phi.formula = ~ Treatment,
                                     formula_null = ~ 1,
                                     phi.formula_null = ~ Treatment,
                                     test= "Wald",
                                     boot = FALSE,
                                     data = root_t3,
                                     fdr_cutoff = 0.05)
Genus_t3_analysis$significant_taxa
# No significant differences! Also expected after inspection of the heatmaps

## Combine the data from the plots in order to make something better

t0_data<-t0_plot$data
t0_data$time<-"T0"
t1_data<-t1_plot$data
t1_data$time<-"T1"
t2_data<-t2_plot$data
t2_data$time<-"T2"

all_data<-rbind(t0_data,t1_data,t2_data)
library(stringi)
etall<-lapply(strsplit(all_data$taxa,split = "_"),`[`)
taxa_list<-do.call(rbind,etall)
taxa_list<-taxa_list[,-1]
colnames(taxa_list)<-c("Phylum","Class","Order","Family","Genus")
taxa_for_plot<-cbind(all_data,taxa_list)

library(ggplot2)

differential_plot<-taxa_for_plot %>%
  ggplot(aes(x=x, y = Genus)) +geom_point()+
  facet_wrap(.~time,scales = "free") + geom_vline(xintercept = 0)+
  geom_errorbar(aes(xmin = xmin, xmax = xmax),width = 0.1)
ggsave("Plots/differential_plot_first_attempt.png",differential_plot, width = 12, heigh = 7)
#Plot used for supplementary material
taxa_for_plot %>%
  ggplot(aes(x=x, y = Genus)) +geom_point()+
   geom_vline(xintercept = 0)+
  geom_errorbar(aes(xmin = xmin, xmax = xmax),width = 0.1)
