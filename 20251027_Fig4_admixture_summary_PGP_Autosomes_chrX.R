
### Reading Q file from RFMIX ver2

path<-"/Users/vborda/Documents/Public_data/INS_data/data+RefPops/"
key<-read.table(paste0(path,"FULL_1KGP_HGDP_PGP_IDs.txt"), head=FALSE, fill = TRUE)

name.file<-"Local_ancestry/4ANC/INS_LDGH_LAI_4ANC_"

for(i in 1:22) { 
  filename <- paste0(path,name.file,"chr",i,".rfmix.Q")
  nam <- paste("df", i, sep = "")
  assign(nam, read.table(filename,comment.char="", head=TRUE,skip=1)  )
}

library(data.table)
tablesum<-rbindlist(list(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18,df19,df20,df21,df22))[, lapply(.SD, sum, na.rm = TRUE), by = X.sample]

tablesum$AFR<-tablesum$AFR/22
tablesum$EUR<-tablesum$EUR/22
tablesum$NAT<-tablesum$NAT/22
tablesum$EAS<-tablesum$EAS/22

colnames(key)<-c("ID","Population")
colnames(tablesum)[1]<-"ID"


demo.prop<-merge(key,tablesum,by.x = "ID",by.y = "ID")

write.table(demo.prop,paste0(path,"Admixture_proportions_INS_LDGH.props"),sep="\t", na="NA",col.names = TRUE, quote = F,row.names = FALSE)


########################################
##### March 2025
########################################

autopropnamefile<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_K4.txt"
xpropnamefile<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_chrX_K4.txt"

autopropnamefile5<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_K5.txt"
xpropnamefile5<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_chrX_K5.txt"
xpropnamefile8<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_proportions_chrX_K8.txt"

autoprop<-read.table(autopropnamefile,header=FALSE,sep="")
xprop<-read.table(xpropnamefile,header=FALSE,sep="\t")

#autoprop5<-read.table(autopropnamefile5,header=FALSE,sep="")
#xprop5<-read.table(xpropnamefile5,header=FALSE,sep="\t")
#xprop8<-read.table(xpropnamefile8,header=FALSE,sep="\t")

colnames(autoprop)<-c("Population","id","AFR","EUR","EAS","NAT")
colnames(xprop)<-c("Population","id","AFR","EUR","EAS","NAT")

#colnames(autoprop6)<-c("Population","id","Amazonian","African","South Andean","North Andean","East Asian","European")

library(dplyr)
library(reshape2)
library(ggplot2)

autoprop$Population[autoprop$Population == "Afro_des"] <- "El Carmen"
autoprop$Population[autoprop$Population == "Ancash"] <- "Huaraz"
autoprop$Population[autoprop$Population == "Puno"] <- "Juliaca"

xprop$Population[xprop$Population == "Afro_des"] <- "El Carmen"
xprop$Population[xprop$Population == "Ancash"] <- "Huaraz"
xprop$Population[xprop$Population == "Puno"] <- "Juliaca"

regions<-c("Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","Lima","El Carmen","Arequipa","Moquegua","Tacna","Ayacucho","Cusco","Juliaca")

autoprop<-autoprop[autoprop$Population %in% regions,]
xprop<-xprop[xprop$Population %in% regions,]

xprop$Population<- factor(xprop$Population, levels = regions)

autoprop$Population<- factor(autoprop$Population, levels = c("Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","Lima",
                                                             "El Carmen","Arequipa","Moquegua","Tacna","Ayacucho","Cusco","Juliaca"))


automeanprops <- autoprop %>%
  group_by(Population) %>%
  summarise(across(c(AFR, EUR, NAT, EAS), 
                   list(mean = mean, sd= sd),
                        #se = ~sd(.) / sqrt(n())),
                   .names = "{col}_Autosomes_{fn}"))


#automeanprops<-autoprop %>%
#  group_by(Population) %>%
#  summarise_at(vars(AFR, EUR,NAT,EAS), list(name = mean))

xmeanprops<-xprop %>%
  group_by(Population) %>%
  summarise(across(c(AFR, EUR, NAT, EAS), 
               list(mean = mean, sd= sd,se = ~sd(.) / sqrt(n())),
               .names = "{col}_X_{fn}"))



#colnames(automeanprops)<-c("Population","AFR_A","EUR_A","NAT_A","EAS_A")
#colnames(xmeanprops)<-c("Population","AFR_X","EUR_X","NAT_X","EAS_X")
#colnames(xprop)<-c("Population","id","AFR_X","EUR_X","EAS_X","NAT_X")


#props.full<-merge(xprop,autoprop)

props.table<-merge(xmeanprops,automeanprops)

colnames(xprop)<-c("Population","id","AFR_X","EUR_X","EAS_X","NAT_X")
props.full<-merge(xprop,autoprop)

write.csv(props.table,"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/202503_admixture_mean_proportions_auto_chrX_K4.txt", quote=FALSE,row.names = FALSE)
#write.csv(totest,"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/admixture_mean_proportions_auto_chrX_K4.txt", quote=FALSE,row.names = FALSE)


library(dplyr)
library(tidyr)

# Define ancestries to test
ancestries <- c("AFR", "EUR", "NAT", "EAS")

# Initialize an empty data frame to store results
wilcoxon_results <- data.frame()

# Loop through each ancestry and perform Wilcoxon test per population
for (anc in ancestries) {
  results <- props.full %>%
    group_by(Population) %>%
    summarise(
      p_value = wilcox.test(.data[[anc]], .data[[paste0(anc, "_X")]], 
                            paired = TRUE, exact = FALSE)$p.value,
      .groups = "drop"
    ) %>%
    mutate(Ancestry = anc)  # Add ancestry as a column
  
  # Append results to the data frame
  wilcoxon_results <- bind_rows(wilcoxon_results, results)
}

# Reshape the data to have ancestries as columns and populations as rows
wilcoxon_results_wide <- wilcoxon_results %>%
  pivot_wider(names_from = Ancestry, values_from = p_value)

# Save results as a CSV file
write.csv(wilcoxon_results_wide,
          "/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/wilcoxon_results_wide.csv", row.names = FALSE)

# Print first few rows to confirm
head(wilcoxon_results_wide)






library(dplyr)

props.table <- props.table %>%
  mutate(
    AFR_Z = (AFR_Autosomes_mean - AFR_X_mean) / sqrt(AFR_X_sd^2 + AFR_Autosomes_sd^2),
    EUR_Z = (EUR_Autosomes_mean - EUR_X_mean) / sqrt(EUR_X_sd^2 + EUR_Autosomes_sd^2),
    NAT_Z = (NAT_Autosomes_mean - NAT_X_mean) / sqrt(NAT_X_sd^2 + NAT_Autosomes_sd^2),
    EAS_Z = (EAS_Autosomes_mean - EAS_X_mean) / sqrt(EAS_X_sd^2 + EAS_Autosomes_sd^2)
  )

write.csv(props.table,"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/proportions/202503_admixture_mean_proportions_auto_chrX_K4.txt", quote=FALSE,row.names = FALSE)


#### Plotting both autosomal and chromosome x



library(ggExtra)

p <- ggplot(props.full, aes(x = AFR, y = AFR_X, color = Population, size = Population)) + 
  geom_point(alpha = 0.7) + 
  theme_minimal() + 
  theme(legend.position = "none")
ggMarginal(p, type = "boxplot") 


props.full2<-melt(props.full)
props.full2$Population<-as.character(props.full2$Population)
#props.full2$Population[props.full2$Population == "Afro_des"] <- "Ica"
props.full2$color[props.full2$variable == 'AFR' | props.full2$variable == 'AFR_X' ] <- 'blue'
props.full2$color[props.full2$variable == 'EUR' | props.full2$variable == 'EUR_X' ] <- 'darkorchid' #'firebrick'
props.full2$color[props.full2$variable == 'EAS' | props.full2$variable == 'EAS_X' ] <- 'gold4'
props.full2$color[props.full2$variable == 'NAT' | props.full2$variable == 'NAT_X' ] <- 'lightgreen'


props.full2$Ancestry[props.full2$variable == 'AFR' | props.full2$variable == 'AFR_X' ] <- 'African'
props.full2$Ancestry[props.full2$variable == 'EUR' | props.full2$variable == 'EUR_X' ] <- 'European'
props.full2$Ancestry[props.full2$variable == 'EAS' | props.full2$variable == 'EAS_X' ] <- 'East Asian'
props.full2$Ancestry[props.full2$variable == 'NAT' | props.full2$variable == 'NAT_X' ] <- 'Indigenous American'

props.full2$chromosome[props.full2$variable == 'AFR_X' | props.full2$variable == 'EUR_X'|
                       props.full2$variable == 'NAT_X' | props.full2$variable == 'EAS_X' ] <- 'Chromosome X'
props.full2$chromosome[props.full2$variable == 'AFR' | props.full2$variable == 'EUR'|
                       props.full2$variable == 'NAT' | props.full2$variable == 'EAS' ] <- 'Autosomes'



props.full2$xmin[props.full2$chromosome== 'Chromosome X' ] <- 0
props.full2$xmax[props.full2$chromosome== 'Chromosome X' ] <- 0.5
props.full2$xmin[props.full2$chromosome== 'Autosomes' ] <- 0.5
props.full2$xmax[props.full2$chromosome== 'Autosomes' ] <- 1

props.full2$colorback[props.full2$chromosome== 'Chromosome X' ] <- "white"
props.full2$colorback[props.full2$chromosome== 'Autosomes' ] <- "black"


#ggplot(data = props.full2, aes(x=Population, y=value)) + geom_boxplot(aes(fill=variable))

source("/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/shift_legend.R")
library(ggplot2)
library(grid)

p<-ggplot(data = props.full2, aes(x=Population, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_manual(values=c("dodgerblue1","darkorchid","gold2","lightgreen","dodgerblue1","darkorchid","gold2","lightgreen")) +
  #scale_fill_manual(values=c("dodgerblue1","firebrick","gold2","lightgreen","dodgerblue1","firebrick","gold2","lightgreen")) +
  facet_wrap( ~ Population, scales="free") + ylab("Genome-wide Ancestry Proportion") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ xlab ("Population") + guides(fill = guide_legend(title.position = "top",
                                                                                        label.position = "bottom",
                                                                                        ncol=4)) +
  theme(legend.direction = "horizontal") + labs(fill='Ancestry') + 
  annotate("rect", xmin=0.5, xmax=1, ymin=0, ymax=Inf, alpha=0.2, fill="black")
  #annotate("rect", xmin=0.5, xmax=1, ymin=0, ymax=Inf, alpha=0.1, fill="purple")


props.full2$variable <- factor(props.full2$variable, 
                               levels = unique(props.full2$variable),
                               labels = rep(c("African", "European", "Indigenous American", "East Asian"), 2) # Adjust mapping
)



ggplot(data = props.full2, aes(x=Population, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_manual(
    values = c("dodgerblue1", "darkorchid", "gold2", "lightgreen"),
    #values = c("dodgerblue1", "firebrick", "gold2", "lightgreen"),darkorchid
    labels = c("African", "European", "Indigenous American", "East Asian") # Custom labels
  ) +
  facet_wrap(~ Population, scales = "free") +
  ylab("Genome-wide Ancestry Proportion") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + 
  xlab("Population") + 
  guides(fill = guide_legend(
    title = "Ancestry",
    title.position = "top",
    label.position = "bottom",
    ncol = 4
  )) + 
  theme(legend.direction = "horizontal") + 
  annotate("rect", xmin = 0.5, xmax = 1, ymin = 0, ymax = Inf, alpha = 0.2, fill = "black")



path<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/"
  jpeg(filename = paste0(path,"2025Oct_Contrasting_Autosome_ChromosomeX_K4.jpg"),width = 24, height = 18, units = "cm", pointsize =12,  res = 300) 
  #jpeg(filename = paste0(path,"Contrasting_Autosome_ChromosomeX_K4_black.jpg"),width = 24, height = 18, units = "cm", pointsize =12,  res = 300) 
  grid.draw(shift_legend(p))
dev.off()


library(ggplot2)
library(gridExtra)
library(grid)

# Create the plot without legend
base_plot <- ggplot(data = props.full2, aes(x=Population, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_manual(
    values = c("dodgerblue1", "darkorchid", "gold2", "lightgreen"),
    #values = c("dodgerblue1", "firebrick", "gold2", "lightgreen"),
    labels = c("African", "European", "Indigenous American", "East Asian")
  ) +
  facet_wrap(~ Population, scales = "free") +
  ylab("Genome-wide Ancestry Proportion") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + # Moves legend to bottom
  labs(fill = "Ancestry") 

# Extract legend
legend <- cowplot::get_legend(base_plot)

# Generate 13 plots (assuming you store them in a list)
plots <- replicate(13, base_plot, simplify = FALSE) # Replace with actual plots

# Arrange plots and legend
grid.arrange(
  grobs = c(plots, list(legend)),
  ncol = 4, # Adjust layout
  nrow = 4, # 3 rows for plots, 1 for legend
  heights = c(rep(1, 3), 0.2) # Give the last row less space for the legend
)


df<-props.full2
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(patchwork)

# Define ancestry colors
ancestry_colors <- c("Indigenous American" = "lightgreen", 
                     "African" = "dodgerblue1", 
                     "European" = "darkorchid", #"firebrick"
                     "East Asian" = "gold2")

# Function to create a plot for each population
plot_population <- function(pop_name, data) {
  ggplot(data %>% filter(Population == pop_name) %>%
           mutate(chromosome = factor(chromosome)),
         aes(x = chromosome, y = value, fill = Ancestry)) +
    geom_jitter(aes(color = Ancestry),
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                size = 1.2, alpha = 0.7) +  # Draw points first
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.5) +  # Semi-transparent boxes in front
    scale_fill_manual(values = ancestry_colors) +
    scale_color_manual(values = ancestry_colors) +
    labs(title = pop_name, x = "Chromosome", y = "Ancestry Proportion") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = rel(0.8)),
      axis.title = element_text(size = rel(0.8)),
      axis.text = element_text(size = rel(0.8)),
      strip.text = element_text(size = rel(0.8))
    ) +
    geom_vline(xintercept = 1.5, linetype = "dashed", color = "black")
}



# Get unique populations
populations <- unique(df$Population)

# Generate individual plots
plots <- lapply(populations, plot_population, data = df)

# Add empty plots to make a 4x4 grid (13 populations → need 3 more blank spaces)
while (length(plots) < 16) {
  plots <- c(plots, list(ggplot() + theme_void())) # Blank plot for spacing
}

# Arrange plots in a 4x4 grid and collect legend
final_plot <- wrap_plots(plots, ncol = 4) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") # Ensure legend is placed at the bottom

# Display final figure
path<-"/Users/vborda/Documents/Public_data/INS_data/ADMIXTURE_results/"
jpeg(filename = paste0(path,"2025Oct_Contrasting_Autosome_ChromosomeX_K4.jpg"),width = 24, height = 18, units = "cm", pointsize =12,  res = 300) 
#jpeg(filename = paste0(path,"Contrasting_Autosome_ChromosomeX_K4_black.jpg"),width = 24, height = 18, units = "cm", pointsize =12,  res = 300) 
final_plot
dev.off()

#### HIGH RESOLUTION

# Function to create a plot for each population
plot_population <- function(pop_name, data) {
  ggplot(data %>% filter(Population == pop_name) %>%
           mutate(chromosome = factor(chromosome)),
         aes(x = chromosome, y = value, fill = Ancestry)) +
    geom_jitter(aes(color = Ancestry),
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                size = 1.2, alpha = 0.7) +  # Draw points first
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.5) +  # Semi-transparent boxes in front
    scale_fill_manual(values = ancestry_colors) +
    scale_color_manual(values = ancestry_colors) +
    labs(title = pop_name, x = "Chromosome", y = "Ancestry Proportion") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = rel(0.8)),
      axis.title = element_text(size = rel(0.6)),
      axis.text = element_text(size = rel(0.6)),
      strip.text = element_text(size = rel(0.6))
    ) +
    geom_vline(xintercept = 1.5, linetype = "dashed", color = "black")
}



# Get unique populations
populations <- unique(df$Population)

# Generate individual plots
plots <- lapply(populations, plot_population, data = df)

# Add empty plots to make a 4x4 grid (13 populations → need 3 more blank spaces)
while (length(plots) < 16) {
  plots <- c(plots, list(ggplot() + theme_void())) # Blank plot for spacing
}

# Arrange plots in a 4x4 grid and collect legend
final_plot <- wrap_plots(plots, ncol = 4) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")



# Save as TIFF (high-res) at A4
path<-"/Users/vborda/Documents/Public_data/INS_data/Manuscript/R2_Figures/"
tiff(filename = paste0(path, "Figure_4_Autosome_ChromosomeX_K4.tiff"),
     width = 20, height = 15, units = "cm",
     res = 300, pointsize = 12)

final_plot

dev.off()
