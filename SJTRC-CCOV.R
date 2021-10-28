

# SJTRC-CCOV
# R version 4.1.0

library(tidyverse) 
library(ggpubr)
library(rstatix)
library(corrplot)
library(ggsci)
library(vegan)
library(factoextra)
require(RColorBrewer)
require(ComplexHeatmap)
require(circlize)
require(digest)
require(cluster)
require(factoextra) 

demo <- readxl::read_xlsx("CCOV-SARS2 Ab data 5-25-21 Combined for Li CLEAN.xlsx", sheet = "Demographics") 
sample <- readxl::read_xlsx("CCOV-SARS2 Ab data 5-25-21 Combined for Li CLEAN.xlsx", sheet = "CCOV") %>% 
  mutate(infected=ifelse(`Collection day (post diagnosis)`!="N/A",1,0))
severity <- readxl::read_xlsx("CCOV-SARS2 Ab data 5-25-21 Combined for Li CLEAN.xlsx", sheet = "COVID infection details") %>% 
  mutate(SJTRC.severity=`Revised SJTRC severity (5-25-21)`,
         SJTRC.severity.duration=as.numeric(ifelse(`Revised Duration of reported symptoms 4-22-21 (days)`=="Missing",NA,`Revised Duration of reported symptoms 4-22-21 (days)`))) 
x.list <- readxl::read_xlsx("List of subject with a baseline sample MM.xlsx") %>% select(`Exclude these from both severity and hCOV boost analysis`) %>% 
  mutate(exclude=`Exclude these from both severity and hCOV boost analysis`) %>% filter(!is.na(exclude))%>% pull(exclude)

# Baseline
baseline.sample <- sample %>% filter(`Sample type`=="Baseline") %>% arrange(SJTRCID, `Collection day (post enrollment)`)
multi.baseline.ID <- baseline.sample %>% group_by(SJTRCID) %>% summarise(Freq=n()) %>% filter(Freq>1) 
multi.baseline <- baseline.sample %>% filter(SJTRCID %in% multi.baseline.ID$SJTRCID) 
multi.baseline %>% select(SJTRCID,`Sample ID`,`Sample type`,`Collection day (post enrollment)`,`Collection day (post diagnosis)`) 
baseline.sample <- sample %>% filter(`Sample type`=="Baseline") %>% group_by(SJTRCID) %>% arrange(SJTRCID, `Collection day (post enrollment)`) %>% filter(row_number()==1)
baseline.sample %>% filter(SJTRCID %in% multi.baseline.ID$SJTRCID) 
length(unique(baseline.sample$SJTRCID))      #1202

# Demo + Sample
baseline.alldata <- baseline.sample %>% left_join(demo) %>% 
  mutate(race=gsub("^(.*?),.*", "\\1",ifelse(Race=="American Indian and Alaska Native", "Other", Race))) %>% 
  mutate(age=factor(cut(`Age at enrollment`, breaks=c(0,43,1000), right=F, labels=c("Younger","Older")))) %>% 
  mutate(patient.contact=factor(`Patient contact`), age.enroll=`Age at enrollment`) %>% 
  mutate(Gender=factor(Gender)) %>% 
  mutate(race=factor(ifelse(race=="99999",NA,race))) %>% 
  mutate(infection=factor(infected,levels=c("0","1"), labels=c("Not infected","Infected"))) 
cov.list <- baseline.alldata %>% select(starts_with(c("OC43","NL63","229E","HKU1"))) 
cols <- names(cov.list[,-1])
baseline.alldata[paste0("log10 ", cols)] <- log10(baseline.alldata[cols])

outlist <- baseline.alldata %>% select(SJTRCID, `Sample ID`, infection)
table(outlist$infection)

# Data for Plots
baseline.selected <- baseline.alldata %>% 
  select(SJTRCID,`Age at enrollment`,starts_with(c("OC43","NL63","229E","HKU1"))) %>% 
  select(SJTRCID,`Age at enrollment`,ends_with(c("IgG","IgM","IgA"))) 
baseline.data.long <- gather(baseline.selected[,-2], CoV.abs, titers, -SJTRCID) %>% 
  left_join(select(baseline.alldata, SJTRCID, infection, Gender, patient.contact, `Age at enrollment`, race)) %>% 
  mutate(age=factor(cut(`Age at enrollment`, breaks=c(0,43,1000), right=F, labels=c("Younger","Older")))) %>% 
  mutate(age.enroll=`Age at enrollment`) %>% 
  mutate(logtiter=log10(titers)) %>% 
  mutate(ccov=word(CoV.abs, 1), type=word(CoV.abs, 2)) %>% 
  mutate(type=factor(type,levels = c("IgG","IgM","IgA")), ccov=factor(ccov,levels = c("OC43", "HKU1", "229E","NL63"))) %>% 
  arrange(type, ccov) 
orders <- unique(baseline.data.long$CoV.abs)
baseline.data.long$CoV.abs <- factor(baseline.data.long$CoV.abs, levels = orders)

# correlation analysis
cor.matrix <- cor(baseline.selected[,-c(1,2)], method="spearman", use = "pairwise.complete.obs")
pval <- psych::corr.test(baseline.selected[,-c(1,2)], method="spearman", adjust="bonferroni")$p
tiff(file = "Fig 1E.tiff", width=2000, height=2000, res=300)
corrplot(cor.matrix, method="color", type="upper", 
         p.mat=pval, tl.col = 'black', addCoef.col = "white", 
         col=colorRampPalette(c("blue","yellow","red"))(20),
         diag = F, mar=c(0,0,0,0), tl.srt = 45,title = "", cl.cex = 1.025)
dev.off()

# Age
stat.test <- baseline.data.long[complete.cases(baseline.data.long$logtiter),] %>%
  group_by(CoV.abs) %>%
  wilcox_test(logtiter ~ age) %>%  
  mutate(p.adj=p.adjust(p, method="bonferroni")) %>% 
  add_significance("p.adj") %>% 
  add_xy_position(x = "age", dodge = 1, step.increase = .1)
ggboxplot(baseline.data.long[complete.cases(baseline.data.long$logtiter),], 
          x = "age", y = "logtiter", fill="age", palette = "jco", facet.by = "CoV.abs", ncol=4, scales="free") +
  scale_y_continuous("Log10 Titers",expand = expansion(add = 0.2)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif") +
  theme_pubr() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        strip.text = element_text(size=15),
        legend.position = "none",
        panel.border = element_rect(color ="black", fill=NA)) 
ggsave("Fig 2A.tiff", width = 12, height = 8, dpi = 300)

# Patient contact
stat.test <- baseline.data.long[complete.cases(baseline.data.long$logtiter),] %>%
  group_by(CoV.abs) %>%
  wilcox_test(logtiter ~ patient.contact) %>%  
  mutate(p.adj=p.adjust(p, method="bonferroni")) %>% 
  add_significance("p.adj") %>% 
  add_xy_position(x = "patient.contact", dodge = 1, step.increase = .1)
ggboxplot(baseline.data.long[complete.cases(baseline.data.long$logtiter),], 
          x = "patient.contact", y = "logtiter", fill="patient.contact", palette = "jco", facet.by = "CoV.abs", ncol=4, scales="free") +
  scale_y_continuous("Log10 Titers",expand = expansion(add = 0.2)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif") +
  theme_pubr() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        strip.text = element_text(size=15),
        legend.position = "none",
        panel.border = element_rect(color ="black", fill=NA))  
ggsave("Fig 2B.tiff", width = 12, height = 8, dpi = 300)

# Gender
stat.test <- baseline.data.long[complete.cases(baseline.data.long$logtiter),] %>%
  group_by(CoV.abs) %>%
  wilcox_test(logtiter ~ Gender) %>%  
  mutate(p.adj=p.adjust(p, method="bonferroni")) %>% 
  add_significance("p.adj") %>%
  add_xy_position(x = "Gender", dodge = 1, step.increase = .1)
ggboxplot(baseline.data.long[complete.cases(baseline.data.long$logtiter),], 
          x = "Gender", y = "logtiter", fill="Gender", palette = "jco", facet.by = "CoV.abs", ncol=4, scales="free") +
  scale_y_continuous("Log10 Titers",expand = expansion(add = 0.2)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif") +
  theme_pubr() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        strip.text = element_text(size=15),
        legend.position = "none",
        panel.border = element_rect(color ="black", fill=NA)) 
ggsave("Fig S1A.tiff", width = 12, height = 8, dpi = 300)

# race
race.complete.data <- baseline.data.long[complete.cases(baseline.data.long$logtiter),] %>% filter(!is.na(race))
stat.test <- race.complete.data %>% 
  group_by(CoV.abs) %>%
  wilcox_test(logtiter ~ race) %>%  
  mutate(p.adj=p.adjust(p, method="bonferroni")) %>% 
  add_significance("p.adj") %>%
  add_xy_position(x = "race", dodge = 1, step.increase = .1)
ggboxplot(race.complete.data, 
          x = "race", y = "logtiter", fill="race", palette = "jco", facet.by = "CoV.abs", ncol=4, scales="free") +
  scale_y_continuous("Log10 Titers",expand = expansion(add = 0.2)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif") +
  theme_pubr() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        strip.text = element_text(size=15),
        legend.position = "none",
        panel.border = element_rect(color ="black", fill=NA)) 
ggsave("Fig S1B.tiff", width = 12, height = 8, dpi = 300)

# infection
stat.test <- baseline.data.long[complete.cases(baseline.data.long$logtiter),] %>%
  group_by(CoV.abs) %>%
  wilcox_test(logtiter ~ infection) %>%  
  mutate(p.adj=p.adjust(p, method="bonferroni")) %>% 
  add_significance("p.adj") %>% 
  add_xy_position(x = "infection", dodge = 1, step.increase = .1)
ggboxplot(baseline.data.long[complete.cases(baseline.data.long$logtiter),], 
          x = "infection", y = "logtiter", fill="infection", palette = "jco", facet.by = "CoV.abs", ncol=4, scales="free") +
  scale_y_continuous("Log10 Titers",expand = expansion(add = 0.2)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif") +
  theme_pubr() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        strip.text = element_text(size=15),
        legend.position = "none",
        panel.border = element_rect(color ="black", fill=NA)) 
ggsave("Fig 4.tiff", width = 12, height = 8, dpi = 300)

# Baseline
baseline.sample <- sample %>% filter(!SJTRCID %in% x.list) %>% filter(`Collection day (post diagnosis)`!="N/A") %>% filter(`Sample type`=="Baseline") 
multi.baseline.ID <- baseline.sample %>% group_by(SJTRCID) %>% summarise(Freq=n()) %>% filter(Freq>1); multi.baseline.ID
baseline.sample <- baseline.sample %>% group_by(SJTRCID) %>% arrange(SJTRCID, `Collection day (post enrollment)`) %>% filter(row_number()==1) 
baseline.sample %>% nrow()
# Acute
acute.sample <- sample %>% filter(!SJTRCID %in% x.list) %>% filter(`Collection day (post diagnosis)`!="N/A") %>% filter(`Sample type`=="Acute")
multi.acute.ID <- acute.sample %>% group_by(SJTRCID) %>% summarise(Freq=n()) %>% filter(Freq>1); multi.acute.ID
acute.sample %>% nrow()
# Convalescent
convalescent.sample <- sample %>% filter(!SJTRCID %in% x.list) %>% filter(`Collection day (post diagnosis)`!="N/A") %>% filter(`Sample type`=="Convalescent")
multi.convalescent.ID <- convalescent.sample %>% group_by(SJTRCID) %>% summarise(Freq=n()) %>% filter(Freq>1); multi.convalescent.ID
convalescent.sample <- convalescent.sample %>% group_by(SJTRCID) %>% arrange(SJTRCID, `Collection day (post enrollment)`) %>% filter(row_number()==1) 
convalescent.sample %>% nrow()
# baseline 2+ data 
all.sample <- rbind(baseline.sample,acute.sample,convalescent.sample) %>% ungroup
check.table <- table(all.sample$SJTRCID,all.sample$`Sample type`) %>% rbind() %>% 
  data.frame() %>% rownames_to_column("SJTRCID") %>%  mutate(nsample=Baseline+Acute+Convalescent)
id.full <- check.table %>% filter(nsample==3) %>% select(SJTRCID); nrow(id.full) 
id.b2plus <- check.table %>% filter(Baseline==1 & nsample>1) %>% select(SJTRCID); nrow(id.b2plus) 
id.12 <- check.table %>% filter(Baseline==1 & Acute==1) %>% select(SJTRCID); nrow(id.12) 
id.13 <- check.table %>% filter(Baseline==1 & Convalescent==1) %>% select(SJTRCID); nrow(id.13) 
sample.long <- all.sample %>% left_join(check.table) %>% 
  filter(Baseline==1 & nsample>1) %>% 
  select(SJTRCID,`Sample type`,ends_with(c("IgG","IgM","IgA"))) %>% 
  gather(Abs, titers, -c(SJTRCID, `Sample type`)) %>% 
  mutate(Status=factor(`Sample type`,levels=c("Baseline","Acute","Convalescent"))) %>% 
  mutate(logtiter=log10(titers)) %>% 
  select(-c(titers,`Sample type`)) %>% 
  filter(!is.na(logtiter))
sample.test <- spread(sample.long, Status, logtiter)

# 1 and 2
sample.prep <- sample.test %>% filter(!is.na(Baseline) & !is.na(Acute)) %>% mutate(pc=(10^(Acute-Baseline)-1)*100) 
summary(sample.prep$pc)
sample.data <- sample.prep %>% select(SJTRCID, Abs, pc) %>% spread(Abs, pc) %>% 
  filter(complete.cases(.)) %>% 
  select(SJTRCID, starts_with(c("OC43","HKU1","229E","NL63"))) %>% 
  select(SJTRCID, ends_with(c("IgG","IgM","IgA"))) 
groups <- data.frame(cols=names(sample.data)[-1]) %>% 
  mutate(ccov=factor(word(cols, 1),levels = c("OC43", "HKU1", "229E","NL63")), 
         type=factor(word(cols, 2),levels=c("IgG","IgM","IgA")))
sample.fin <- sample.data %>% column_to_rownames(var="SJTRCID") %>% as.matrix() 
missing.list <- id.12 %>% filter(!SJTRCID %in% sample.data$SJTRCID)
all.sample %>% filter(SJTRCID %in% missing.list$SJTRCID) %>% t() %>% data.frame() %>% 
  rownames_to_column() %>% 
  mutate(n=row_number(), g=ifelse(n>3,0,NA)) %>% 
  filter(!complete.cases(.)) %>% 
  select(-n,-g)
ht_opt(legend_title_gp = gpar(fontsize = 15),legend_labels_gp = gpar(fontsize = 15))
mycols <- c("cyan","white","magenta", "red")
mybreaks <- c(-50, 0, 100, 300)
col = list(hCoV = c("OC43" = brewer.pal(12, 'Paired')[9], "HKU1" = brewer.pal(12, 'Paired')[10],
                    "229E" = brewer.pal(12, 'Paired')[7],"NL63" = brewer.pal(12, 'Paired')[8]))
anno <- HeatmapAnnotation(hCoV = groups$ccov, col = col, annotation_name_side = "left", annotation_name_gp= gpar(fontsize = 15))
ht1 <- Heatmap(sample.fin, col = colorRamp2(mybreaks, mycols),
               top_annotation = anno,
               cluster_columns = FALSE,
               column_split=groups$type,
               column_title_gp = gpar(fontsize = 16), 
               show_row_dend = FALSE, 
               show_row_names = FALSE,
               show_column_names = FALSE,
               clustering_distance_rows = "euclidean", 
               clustering_method_rows = "ward.D2",
               heatmap_legend_param = list(at=mybreaks, title = "Percent Change"))
tiff(file = "Fig 3A.tiff", width=10*150, height=20*150, res = 300)
draw(ht1, merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()

# 1 and 3
sample.prep <- sample.test %>% filter(!is.na(Baseline) & !is.na(Convalescent)) %>% 
  filter(!SJTRCID %in% c("5D110532","88D43EFC")) %>% 
  mutate(pc=(10^(Convalescent-Baseline)-1)*100) 
summary(sample.prep$pc)
sample.data <- sample.prep %>% select(SJTRCID, Abs, pc) %>% spread(Abs, pc) %>% 
  filter(complete.cases(.)) %>% 
  select(SJTRCID, starts_with(c("OC43","HKU1","229E","NL63"))) %>% 
  select(SJTRCID, ends_with(c("IgG","IgM","IgA"))) 
groups <- data.frame(cols=names(sample.data)[-1]) %>% 
  mutate(ccov=factor(word(cols, 1),levels = c("OC43", "HKU1", "229E","NL63")), 
         type=factor(word(cols, 2),levels=c("IgG","IgM","IgA")))
sample.fin <- sample.data %>% column_to_rownames(var="SJTRCID") %>% as.matrix() 
mycols <- c("cyan","white","magenta", "red")
mybreaks <- c(-50, 0, 100, 300)
col = list(hCoV = c("OC43" = brewer.pal(12, 'Paired')[9], "HKU1" = brewer.pal(12, 'Paired')[10],
                    "229E" = brewer.pal(12, 'Paired')[7],"NL63" = brewer.pal(12, 'Paired')[8]))
anno <- HeatmapAnnotation(hCoV = groups$ccov, col = col, annotation_name_side = "left", annotation_name_gp= gpar(fontsize = 15))
ht2 <- Heatmap(sample.fin, col = colorRamp2(mybreaks, mycols),
               top_annotation = anno,
               cluster_columns = FALSE,
               column_split=groups$type,
               column_title_gp = gpar(fontsize = 16), 
               show_row_dend = FALSE, 
               show_row_names = FALSE,
               show_column_names = FALSE,
               clustering_distance_rows = "euclidean", 
               clustering_method_rows = "ward.D2",
               heatmap_legend_param = list(at=mybreaks, title = "Percent Change"))
tiff(file = "Fig 3B.tiff", width=10*150, height=20*150, res = 300)
draw(ht2, merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()

# figure 5A
baseline.sample <- sample %>% filter(`Sample type`=="Baseline") %>% arrange(SJTRCID, `Collection day (post enrollment)`)
multi.baseline.ID <- baseline.sample %>% group_by(SJTRCID) %>% summarise(Freq=n()) %>% filter(Freq>1) 
multi.baseline <- baseline.sample %>% filter(SJTRCID %in% multi.baseline.ID$SJTRCID) 
baseline.sample <- sample %>% filter(`Sample type`=="Baseline") %>% group_by(SJTRCID) %>% arrange(SJTRCID, `Collection day (post enrollment)`) %>% filter(row_number()==1)
baseline.alldata <- baseline.sample %>% 
  filter(`Collection day (post diagnosis)`!="N/A") %>% 
  filter(!SJTRCID %in% x.list) %>% 
  left_join(demo) %>% 
  left_join(severity) %>% 
  mutate(race=gsub("^(.*?),.*", "\\1",ifelse(Race=="American Indian and Alaska Native", "Other", Race))) %>% 
  mutate(severity=ifelse(!SJTRC.severity %in% c("1","2","2.5","3.5","4"),NA, ifelse(SJTRC.severity %in% c("2.5","3.5","4"), 1, 0)),
         severity_5=ifelse(!SJTRC.severity %in% c("1","2","2.5","3.5","4"),NA, as.numeric(SJTRC.severity))) %>% 
  mutate(patient.contact=factor(`Patient contact`), age.enroll=`Age at enrollment`) %>% 
  mutate(Gender=factor(Gender)) %>% 
  mutate(race=factor(ifelse(race=="99999",NA,race))) %>% 
  mutate(white=factor(ifelse(is.na(race),NA,ifelse(race=="White","White","Non-White")))) %>% 
  mutate(severity_f=factor(severity,c("1","0"),labels=c("Severe/Critical","Asymptomatic/Mild/Moderate"))) %>%
  mutate(severity_5_f=factor(severity_5,c("1","2","2.5","3.5","4"),labels=c("1","2","3","4","5"))) 
baseline.fin <- baseline.alldata %>% filter(!is.na(severity)) %>%  mutate(bmi=as.numeric(ifelse(BMI=="N/A", NA, BMI)))
cov.list <- baseline.fin %>% select(starts_with(c("OC43","NL63","229E","HKU1"))) 
cols <- names(cov.list[,-1])
baseline.fin[cols] <- log10(baseline.fin[cols])
mat <- baseline.fin %>% 
  select(SJTRCID,ends_with(c("IgA","IgM","IgG"))) %>% 
  remove_rownames %>% column_to_rownames(var="SJTRCID") %>% as.matrix() %>% t()
mycols<-colorRampPalette(c("white","yellow","red"), space="rgb")(100)
mybreaks <- seq(0.5, 3, length.out = 100)
col.Sex <- brewer.pal(8, 'Set2')
age.range <- summary(baseline.fin$age.enroll)[c(1,3,6)] %>% as.vector()
race.col <- brewer.pal(9, 'Set3')
duration.range <- summary(baseline.fin$SJTRC.severity.duration)[c(1,3,6)] %>% as.vector()
col = list(
  Age = colorRamp2(age.range,brewer.pal(n=3, name="Purples")), 
  Sex = c("Female" = col.Sex[6], "Male" = col.Sex[1]),
  Race = c("White" = race.col[2], "Black" = race.col[4], "Asian" = race.col[5], "Other" =race.col[7]),
  `Patient Contact` = c("Direct" = "black", "Indirect" = "orange", "None"="lightgray"),
  `Days of Symptom` = colorRamp2(duration.range,c("cyan","white","magenta")), 
  `Severity of Symptoms` = c("1" = "lightblue1", "2" = "deepskyblue", "3"="lightpink", "4" = "deeppink", "5" = "brown"),
  `Severity of Symptoms (grouped)` = c("Severe/Critical" = "red", "Asymptomatic/Mild/Moderate" = "blue"))
ht_opt(legend_title_gp = gpar(fontsize = 12),legend_labels_gp = gpar(fontsize = 12))
anno <- HeatmapAnnotation(
  `Severity of Symptoms (grouped)` = baseline.fin$severity_f,
  `Severity of Symptoms` = baseline.fin$severity_5_f,
  `Days of Symptom` =  baseline.fin$SJTRC.severity.duration,
  Age = baseline.fin$age.enroll,
  Sex = baseline.fin$Gender,
  Race = baseline.fin$race,
  `Patient Contact`  = baseline.fin$`Patient contact`, 
  col = col,
  annotation_name_side = "left"
)
ht <- Heatmap(mat, name="Log10 Antibody Titers", 
              col = colorRamp2(mybreaks, mycols), 
              row_names_side = "left", row_dend_side = "left",
              top_annotation = anno,
              clustering_distance_rows = "spearman",
              clustering_method_rows = "ward.D2",
              column_split=baseline.fin$severity_f,
              column_title = " ",
              clustering_distance_columns = "euclidean", 
              show_column_names = FALSE) 

tiff(file = "Fig 5A.tiff", width=16*300, height=6*300, res = 300)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# figure 5C
baseline.sample <- sample %>% filter(`Sample type`=="Baseline") %>% group_by(SJTRCID) %>% arrange(SJTRCID, `Collection day (post enrollment)`) %>% filter(row_number()==1)
baseline.alldata <- baseline.sample %>% 
  filter(`Collection day (post diagnosis)`!="N/A") %>% 
  filter(!SJTRCID %in% x.list) %>% 
  left_join(demo) %>% 
  left_join(severity) %>% 
  mutate(race=gsub("^(.*?),.*", "\\1",ifelse(Race=="American Indian and Alaska Native", "Other", Race))) %>% 
  mutate(severity=ifelse(!SJTRC.severity %in% c("1","2","2.5","3.5","4"),NA, ifelse(SJTRC.severity %in% c("2.5","3.5","4"), 1, 0)),
         severity_5=ifelse(!SJTRC.severity %in% c("1","2","2.5","3.5","4"),NA, as.numeric(SJTRC.severity))) %>% 
  mutate(patient.contact=factor(`Patient contact`), age.enroll=`Age at enrollment`) %>% 
  mutate(Gender=factor(Gender)) %>% 
  mutate(race=factor(ifelse(race=="99999",NA,race))) %>% 
  mutate(white=factor(ifelse(is.na(race),NA,ifelse(race=="White","White","Non-White")))) %>% 
  mutate(severity_f=factor(severity,c("1","0"),labels=c("Severe/Critical","Asymptomatic/Mild/Moderate"))) %>%
  mutate(severity_5_f=factor(severity_5,c("1","2","2.5","3.5","4"),labels=c("1","2","3","4","5"))) 
baseline.fin <- baseline.alldata %>% filter(!is.na(severity)) %>% 
  mutate(bmi=as.numeric(ifelse(BMI=="N/A", NA, BMI)))
cov.list <- baseline.fin %>% select(starts_with(c("OC43","NL63","229E","HKU1"))) 
cols <- names(cov.list[,-1])
baseline.fin[paste0("log10 ", cols)] <- log10(baseline.fin[cols])
baseline.pca.data <- baseline.fin %>% select(SJTRCID,severity_f, severity_5_f,SJTRC.severity.duration, starts_with(c("log10"))) %>% ungroup()
PCA.alldata <- baseline.pca.data %>% select(SJTRCID, severity_f, severity_5_f, SJTRC.severity.duration, starts_with(c("log10"))) %>% column_to_rownames(var="SJTRCID") 
PCA.abs <- PCA.alldata %>% select(starts_with(c("log10")))
basic_plot <- fviz_pca_ind(prcomp(PCA.abs, scale=TRUE), geom.ind = "point")   
basic_plot_out <- basic_plot$data %>% mutate(SJTRCID=name) %>% select(-name)
PCA_plot <- baseline.pca.data %>% left_join(basic_plot_out) %>% data.frame()
ggplot(data = PCA_plot, aes(x, y)) +
  geom_vline(xintercept=0, lty=2)+
  geom_hline(yintercept=0, lty=2)+
  stat_ellipse(geom = "polygon", aes(x, y, fill=severity_f), alpha=0.1, show.legend = FALSE) +
  geom_point(aes(x, y, col=severity_5_f, size=SJTRC.severity.duration)) +
  scale_fill_manual(values = c("red","blue")) +
  scale_color_manual(values =  c("lightblue1", "deepskyblue", "lightpink", "deeppink", "brown")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16)) +
  labs(x=basic_plot$labels$x, y=basic_plot$labels$y, size="Symptom Duration (Days)", color="SJTRC Severity") +
  ggtitle("PCA: Baseline IgA, IgM and IgG")  
ggsave("Fig 5C.tiff", width = 8, height = 6, dpi = 300)

# add duration
baseline.selected <- baseline.fin %>% select(SJTRCID, severity_f, severity_5_f, SJTRC.severity.duration, bmi, patient.contact, age.enroll, Gender, race, starts_with(c("log10"))) %>% 
  mutate(duration_c=factor(cut(SJTRC.severity.duration, breaks=c(-1,0,7,30, 1000), right=T, 
                               labels = c("No symptoms (0 days)", "Symptoms resolved in 7 days", 
                                          "Symptoms resolved in month","Symptoms not resolved in 1 month")),
                           levels=c("No symptoms (0 days)", "Symptoms resolved in 7 days", 
                                    "Symptoms resolved in month","Symptoms not resolved in 1 month")))
baseline.data.long <- baseline.selected %>% 
  select(SJTRCID, starts_with(c("log10"))) %>% 
  gather(CoV.abs, titers, -SJTRCID) %>% 
  left_join(select(baseline.selected, SJTRCID, severity_f, severity_5_f, SJTRC.severity.duration, duration_c)) %>% 
  mutate(ccov=word(CoV.abs, -2), type=word(CoV.abs, -1)) %>% 
  mutate(CoV.abs=paste0(ccov," ",type)) %>% 
  mutate(type=factor(type,levels = c("IgG","IgM","IgA")), ccov=factor(ccov,levels = c("OC43", "HKU1", "229E","NL63"))) %>% 
  arrange(type, ccov) 
glimpse(baseline.data.long)
orders <- unique(baseline.data.long$CoV.abs)
baseline.data.long$CoV.abs <- factor(baseline.data.long$CoV.abs, levels = orders)
stat.test <- baseline.data.long[complete.cases(baseline.data.long$titers),] %>%
  group_by(CoV.abs) %>%
  wilcox_test(titers ~ severity_f) %>%  
  mutate(p.adj=p.adjust(p, method="bonferroni")) %>% 
  add_significance("p.adj") %>% 
  add_xy_position(x = "severity_f", dodge = 1, step.increase = .1)
ggboxplot(baseline.data.long[complete.cases(baseline.data.long$titers),], 
          x = "severity_f", y = "titers", fill="severity_f", palette = "jco", facet.by = "CoV.abs", ncol=4, scales="free") +
  scale_y_continuous("Log10 Titers",expand = expansion(add = 0.2)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif") +
  theme_pubr() +
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        strip.text = element_text(size=15),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        panel.border = element_rect(color ="black", fill=NA)) 
ggsave("Fig 5B.tiff", width = 12, height = 8, dpi = 300)
stat.test <- baseline.data.long[complete.cases(baseline.data.long$titers),] %>%
  group_by(CoV.abs) %>%
  wilcox_test(titers ~ duration_c) %>%  
  mutate(p.adj=p.adjust(p, method="bonferroni")) %>% 
  add_significance("p.adj") %>% 
  add_xy_position(x = "duration_c", dodge = 1, step.increase = .1)
ggboxplot(baseline.data.long[complete.cases(baseline.data.long$titers),], 
          x = "duration_c", y = "titers", fill="duration_c", palette = "jco", facet.by = "CoV.abs", ncol=4, scales="free") +
  scale_y_continuous("Log10 Titers",expand = expansion(add = 0.2)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif") +
  theme_pubr() +
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        strip.text = element_text(size=15),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        panel.border = element_rect(color ="black", fill=NA)) 
ggsave("Fig S4.tiff", width = 12, height = 8, dpi = 300)






















































