### Plotting FST

library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(Hmisc)
library(reshape2)

fst<-read.table("/Users/vborda/Documents/Public_data/INS_data/FST/20250809_FST_INS_PEL.fst.summary",header = TRUE,comment.char = "")
colnames(fst)<-c("Population1","Population2","Hudson_Fst")

fst$Population1[fst$Population1 == "Afro_des"] <- "El Carmen"
fst$Population1[fst$Population1 == "Ancash"] <- "Huaraz"
fst$Population1[fst$Population1 == "Puno"] <- "Juliaca"
fst$Population2[fst$Population2 == "Ica"] <- "El Carmen"
fst$Population2[fst$Population2 == "Ancash"] <- "Huaraz"
fst$Population2[fst$Population2 == "Puno"] <- "Juliaca"

fstnonIca <- fst %>%
  filter(Population1 != "El Carmen" & Population2 != "El Carmen")

regions<-c("Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","PEL","Lima","El Carmen","Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Juliaca")
regions2<-c("Iquitos","Tumbes","Lambayeque","Trujillo","Huaraz","PEL","Lima","Arequipa","Ayacucho","Moquegua","Tacna","Cusco","Juliaca")

fst$Population1 <- factor(fst$Population1, levels=regions)
fst$Population2 <- factor(fst$Population2, levels=regions)
fstnonIca$Population1 <- factor(fstnonIca$Population1, levels=regions2)
fstnonIca$Population2 <- factor(fstnonIca$Population2, levels=regions2)


fst_fixed <- fst %>%
  rowwise() %>%
  mutate(
    # numeric factor codes to compare order
    p1_num = as.numeric(Population1),
    p2_num = as.numeric(Population2),
    # Swap if Population2 comes earlier in order
    pop1_new = ifelse(p1_num <= p2_num, as.character(Population1), as.character(Population2)),
    pop2_new = ifelse(p1_num <= p2_num, as.character(Population2), as.character(Population1))
  ) %>%
  ungroup() %>%
  select(-Population1, -Population2, -p1_num, -p2_num) %>%
  rename(Population1 = pop1_new, Population2 = pop2_new)

fst_unique <- fst_fixed %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(as.character(Population1), as.character(Population2))), collapse = ",")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

fst_unique$Population1 <- factor(fst_unique$Population1, levels=regions)
fst_unique$Population2 <- factor(fst_unique$Population2, levels=regions)

fst_tri <- fst_unique %>%
  filter(as.numeric(Population1) <= as.numeric(Population2))

fstnonIca <- fst_tri %>%
  filter(Population1 != "El Carmen" & Population2 != "El Carmen")

fst3<-fstnonIca %>%  filter(!Population2=='Juliaca')
#fst3<-fstnonIca %>%  filter(!Population1=='Tumbes' | !Population2=='Juliaca')

ALL<-ggplot(fst_tri, aes(x = Population1, y = Population2, fill = Hudson_Fst)) +
  geom_raster() +
  scale_fill_distiller(palette = "Spectral") +
  theme_minimal() +
  labs(x = "Urban Population 1", y = "Urban Population 2", fill = "Hudson's Fst") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

NONICA<-ggplot(fstnonIca, aes(x = Population1, y = Population2, fill = Hudson_Fst)) +
  geom_raster() +
  scale_fill_distiller(palette = "Spectral") +
  theme_minimal() +
  labs(x = "Urban Population 1", y = "Urban Population 2", fill = "Hudson's Fst") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

PLOT3<-ggplot(fst3, aes(x = Population1, y = Population2, fill = Hudson_Fst)) +
  geom_raster() +
  scale_fill_distiller(palette = "Spectral") +
  theme_minimal() +
  labs(x = "Urban Population 1", y = "Urban Population 2", fill = "Hudson's Fst") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#PLOT3<-ggplot(subset3, aes(Population1, Population2, fill= Hudson_Fst)) + 
#  geom_raster() +
#  scale_fill_distiller(palette = "Spectral") +
#  theme_minimal() + xlab("Urban Population 1") +  ylab("Urban Population 2") +
#  theme(panel.grid = element_blank(),
#        axis.text = element_text(size = 6),
#        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +  labs(fill = "Hudson's Fst")

library("ggplot2") 
library("grid") 
library("gridExtra") 
library("cowplot") 
library(ggpubr)

ggarrange(ALL,ggplot() + theme_void(),NONICA, ggplot() + theme_void(),PLOT3,                      # First row with scatter plot
          labels = c("A","", "B","","C"),                  # Second row with box and dot plots
          ncol = 1, heights = c(1, 0.1, 1,0.1,1)       # Labels of the scatter plot
)

ggsave("/Users/vborda/Documents/Public_data/INS_data/Manuscript/R1_Figures/Fig_S3_202508_FST.jpg",
       width = 18, height = 24, units = 'cm', dpi = 300)








library(tidyverse)
library(Hmisc)
library(reshape2)

fst %>% 
  arrange(Population1) %>%
  group_by(Population2) %>%
  filter(row_number() >= which(Population1 == Population2)) %>%
  ggplot(aes(x=Population1, y=Population2, fill=HUDSON_FST)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle=90, hjust=TRUE)) +
  xlab("") + 
  ylab("") +
  geom_text(aes(label=regions), size=txtsize)



