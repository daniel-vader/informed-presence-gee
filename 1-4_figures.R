# Generate descriptive figures
# Author: Daniel Vader
library(ggridges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(viridis)
#library(hrbrthemes)

source("H:/projects/cmv_sot/1-5_functions.R")

cmvdat <- loaddat(sel="cmv")

# Density ridges histogram test ################################################
ggplot(cmvdat[as.numeric(cmvdat$id) < 30,], aes(y=id, x=time, fill=id)) +
  geom_density_ridges(alpha=0.6, stat="binline", bins=30) +
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  scale_fill_viridis(discrete=T)

# Density plots ################################################################
d.demo <- readRDS("Z:/18-015433_CMV in SOT_Downes/05 Data/dan_analytic_data/1_adat_2021-05-20.rds") %>%
  mutate(blood_test_result = factor(blood_test_result, 
                                    levels=c(0,1), 
                                    labels=c("Negative", "Positive")),
         blood_cmv_test = factor(blood_cmv_test, 
                                 levels=c(0,1), 
                                 labels=c("Not tested", "Tested")),
         failure = factor(failure01, 
                          levels=c(0,1), 
                          labels=c("No failure", "Failure")),
         drany = ifelse(cmv_donor_recip_sens > 0, 1, 0),
         drstatus = ifelse(drany > 0, "At Risk", "Low Risk"),
         sot_type = ifelse(sot_type == 0, "Kidney",
                           ifelse(sot_type == 1, "Liver",
                                  ifelse(sot_type == 3 | sot_type == 4, "Lung", "Heart"))),
         ) %>%  
  filter(follow_day > 13)

# Proportion of positive tests
p1 <- ggplot(d.demo, 
             aes(x=follow_day, fill=blood_cmv_test, after_stat(count))) + 
  geom_density(alpha=0.6, 
               position="fill", 
               outline.type = "lower",
               bw="sj") + 
  scale_fill_grey(start=1, end=0.2) +
  ylab("Probability") +
  xlab("Follow-up day") +
  ylim(0,0.05) +
  #labs(fill="Organ", color="Organ") +
  ggtitle("Probability of being tested") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

p1s <- p1 + facet_grid(rows=vars(sot_type), cols=vars(drstatus)) +
  ylim(0,.1) +
  theme(plot.title = element_blank())

# Proportion of positive tests
d.demo.t <-  d.demo %>% filter(blood_cmv_test == "Tested")
p2 <- ggplot(d.demo.t, aes(x=follow_day, fill=blood_test_result, after_stat(count))) + 
  geom_density(alpha=0.6, 
               position="fill", 
               outline.type = "lower",
               bw="sj") + 
  scale_fill_grey(start=1, end=0.2) +
  ylim(0,0.3) +
  ylab("Probability") +
  xlab("Follow-up day") +
  labs(fill="Blood Test Result") +
  ggtitle("Probability of positive test when tested") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

p2s <- p2 + facet_grid(rows=vars(sot_type), cols=vars(drstatus)) + 
  ylim(0,.8) +
  theme(plot.title = element_blank())

tiff(filename = "H:/projects/cmv_sot/figures/testing-patterns.tiff", units="in", width=4, height=6, res=300)
ggpubr::ggarrange(p1,p2, ncol=1,nrow=2, labels = c("A", "B"))
dev.off()

tiff(filename = "H:/projects/cmv_sot/figures/tp-stratified1.tiff", units="in", width=4, height=6, res=300)
p1s
dev.off()

tiff(filename = "H:/projects/cmv_sot/figures/tp-stratified2.tiff", units="in", width=4, height=6, res=300)
p2s
dev.off()

# Heat maps ####################################################################
cmvdat2 <- loaddat(sel="none") 

# Arrange subjects by dropout time to make graph easier to read
cmvdat3 <- cmvdat2 %>%
  group_by(id) %>%
  summarize(endt = max(endt)) %>%
  arrange(endt) %>%
  mutate(id2 = row_number()) %>%
  right_join(cmvdat2, by="id")

# Check dist
scmv <- cmvdat3 %>% 
  group_by(id) %>% 
  summarize(drany = drany[1], sot_type = sot_type[1])

scmv.t <- table(scmv$sot_type, scmv$drany)
  
# Fill in LTFU time points with NA
cmvdat4 <- cmvdat3 %>%
  mutate(drstatus = ifelse(drany > 0, "HR", "LR"),
         sot_type = ifelse(sot_type == 0, "Kidney",
                           ifelse(sot_type == 1, "Liver",
                           ifelse(sot_type == 3 | sot_type == 4, "Lung", "Heart"))),
         drsot = ifelse(drstatus == "HR", 
                        ifelse(sot_type == "Kidney", 1,
                               ifelse(sot_type == "Liver", 3,
                                      ifelse(sot_type == "Lung", 5, 7))),
                        ifelse(sot_type == "Kidney", 2,
                               ifelse(sot_type == "Liver", 4,
                                      ifelse(sot_type == "Lung", 6, 8)))),
         drsotf = factor(drsot, levels=c(1,2,3,4,5,6,7,8),
                        labels=c(paste0("AR Kidney\n(N=", scmv.t[1,2], ")"), 
                                 paste0("LR Kidney\n(N=", scmv.t[1,1], ")"),
                                 paste0("AR Liver\n(N=", scmv.t[2,2], ")"),
                                 paste0("LR Liver\n(N=", scmv.t[2,1], ")"),
                                 paste0("AR Lung\n(N=", scmv.t[4,2] + scmv.t[5,2], ")"),
                                 paste0("LR Lung\n(N=", scmv.t[4,1] + scmv.t[5,1], ")"),
                                 paste0("AR Heart\n(N=", scmv.t[3,2] + scmv.t[6,2], ")"),
                                 paste0("LR Heart\n(N=", scmv.t[3,1] + scmv.t[6,1], ")"))
                        )
         ) %>%
  complete(id2, time) %>%
  group_by(id2) %>%
  mutate(drany = max(drany, na.rm=T),
         sot_type_collapse = sot_type_collapse[1],
         sot_type = sot_type[1],
         drstatus = drstatus[1],
         drsotf = drsotf[1]
         )
  
# Graph
heatmap <- ggplot(cmvdat4, aes(y=as.factor(id2), x=time, fill=test_2wk, color="")) +
  geom_tile() +
  #coord_flip(expand=F) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_blank()
        ) +
  labs(fill = "Tests/\n2 weeks") +
  ylab("Patient") +
  xlab("Time (days)") +
  scale_colour_manual(values=NA) +              
  guides(color=guide_legend("Lost to\nfollow-up", override.aes=list(fill="grey40"))) +
  scale_fill_viridis_c(na.value = "grey40", option="magma")

tiff("H:/projects/cmv_sot/figures/heatmap_all.tiff", 
     res=300, units = "in", width=5, height=5)
heatmap
dev.off()


# Graph by risk group and organ type
tiff("H:/projects/cmv_sot/figures/heatmap_organdr.tiff", 
     res=300, units = "in", width=5, height=6)
heatmap + #coord_flip(expand=F) +
  facet_grid(space="free", scales="free", rows=vars(drsotf)) +
  theme(strip.text.y.right = element_text(angle = 0),
        legend.position = "right")
dev.off()




################################################################################

## Heatmap for just high risk [Kevin's paper] ##
heatmap <- ggplot(cmvdat3[cmvdat3$drany == 1,], aes(y=as.factor(id2), x=time, fill=test_2wk, color="")) +
  geom_tile() +
  #coord_flip(expand=F) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.x = element_blank()
  ) +
  labs(fill = "Tests/\n2 weeks") +
  ylab("Subject") +
  xlab("Time (days)") +
  scale_colour_manual(values=NA) +              
  guides(colour=guide_legend("Lost to\nfollow-up", override.aes=list(colour="grey30"))) +
  scale_fill_viridis_c(na.value = "grey30", option="magma")


tiff("H:/projects/cmv_sot/figures/heatmap_organ-hr-only.tiff", 
     res=300, units = "in", width=5, height=5)
heatmap + #coord_flip(expand=F) +
  facet_grid(space="free", scales="free", rows=vars(sot_type)) +
  theme(strip.text.y.right = element_text(angle = 0),
        legend.position = "right")
dev.off()

# Table 1 ######################################################################
library(tidyverse)
library(furniture)

source("H:/projects/cmv_sot/1-5_functions.R")

cmvdat.tab <- loaddat(sel="cmv") %>% 
  group_by(id) %>%
  summarize(
    sex = sex[1],
    age = age[1],
    reth = race_cat2[1],
    sot_type = sot_type[1],
    drany = drany[1],
    cmvtests = as.numeric(sum(blood_cmv_test)),
    cmv = sum(cmv),
    steroid = sum(steroid_2wk)
  ) %>%
  mutate(sex = ifelse(sex == 1, "Female", "Male"),
         reth = case_when(reth == 1 ~ "White",
                          reth == 3 ~ "Black",
                          reth %in% c(2,9) ~ "Other"),
         sot_type = ifelse(sot_type == 0, "Kidney",
                           ifelse(sot_type == 1, "Liver",
                                  ifelse(sot_type == 3 | sot_type == 4, 
                                         "Lung", "Heart"))),
         drstatus = ifelse(drany > 0, "High Risk", "Low Risk"),
         anycmv = ifelse(cmv > 0, "Yes", "No"))

cvars <- c("sex", "reth", "sot_type", "drstatus", "anycmv")
nvars <- c("age", "cmvtests", "steroid")
tabout <- tibble(var=character(), cat=character(),
                 N = numeric(), Percent = numeric(), upper=numeric())

for(i in 1:length(cvars)){
  t1 <- table(cmvdat.tab[, cvars[i]])
  tp1 <- prop.table(t1)
  for(j in 1:length(t1)){
    tabout <- add_row(tabout, 
                      var=cvars[i], cat=names(t1)[j],
                      N=t1[j], Percent=tp1[j])
  }
}

for(i in 1:length(nvars)){
  d <- cmvdat.tab %>% pull(nvars[i])
  m1 <- median(d)
  lq <- quantile(d, .25)
  uq <- quantile(d, .75)
  tabout <- add_row(tabout,
                    var=nvars[i], N=m1, Percent=lq, upper=uq)
  
}

write.csv(tabout, "H:/projects/cmv_sot/tables/t1.csv")
