library(dplyr)
library(vegan)
library(MASS)
library(mgcv)
library(AER)
library(grid)
library(gridExtra)
library(ggplot2)
library(ggeffects)
library(data.table)
library(tidyr)
library("tidyverse")
library(pairwiseAdonis)
library(openxlsx)

# All resulting figures were later improved using corelDRAW.
#generated csv's were later expanded in excel to create the appendices. 
#PDF files were glued and title pages added.

#setWD

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#preparing variables

flora_dbmean <- read.table("flora_plotlevel_final.csv", header = TRUE, sep = ";", dec = ".", na.strings="NA", strip.white=TRUE)
trem_cat <- read.table("trems_names_explained.csv", header = TRUE, sep = ";", dec = ".", na.strings="NA", strip.white=TRUE)
flora_dbmean$habitat <- as.factor(flora_dbmean$habitat)

#changing the names of columns

setnames(flora_dbmean, old = trem_cat$trem.no,
         new = trem_cat$name)


#preparing categories for columns

trem_forms <- spread(trem_cat, form, name)
trem_forms <- trem_forms[,-c(1,2)]


#calculating richnesses for trem forms

for(i in c(1:length(colnames(trem_forms)))) {
  varname <- paste0("richness_",colnames(trem_forms[i])) 
  group <- as.vector(trem_forms[[i]])
  group <- group[!is.na(group)] 
  flora_dbmean[[varname]] <- rowSums(flora_dbmean[,group] > 0)
}


#### Community differences between the plots ####

emmeans_plot <- function(model, axisy, title) {
  model_emms <- predict_response(model,terms=c("habitat"),margin="empirical")
  plotplease <- ggplot(model_emms, aes(x, predicted, col=x, fill=x)) +
    geom_point(cex=3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width=0.2,linewidth=1) +
    theme_classic() +
    scale_color_manual(values=c("#81b29a","#e07a5f","#3d405b")) +
    xlab("habitat type") +
    ylab(axisy) +
    ggtitle(title)
  print(plotplease)
}

flora_dbmean$habitat <- relevel(flora_dbmean$habitat, c("unaffected"))
flora_dbmean$habitat <- factor(flora_dbmean$habitat, levels = c("unaffected", "unmanaged_disturbance", "logged"))
flora_dbmean <- flora_dbmean %>% mutate(sum.trees.dens = sum.trees/491*10000)
flora_dbmean <- flora_dbmean %>% mutate(dead.trees.dens = dead.tree.no/491*10000)
gam1 <- gam(sum.trees.dens ~ habitat, data = flora_dbmean)
gam2 <- gam(stand.age ~ habitat, data = flora_dbmean)
gam3 <- gam(dead.trees.dens ~ habitat, data = flora_dbmean)
gam4 <- gam(dead.tree.share ~ habitat, data = flora_dbmean)

    #DBH (single tree level)
set.seed(1)
flora_db_circ <- read.table("trems-dbh.csv", header = TRUE, sep = ";", dec = ".", na.strings="NA", strip.white=TRUE)
flora_db_circ$habitat <- factor(flora_db_circ$habitat, levels = c("unaffected", "unmanaged_disturbance", "logged"))
flora_db_circ$plot_ID <- as.factor(flora_db_circ$plot_ID)
gam5 <- gam(DBH ~ habitat + s(plot_ID,bs="re"), data = flora_db_circ)

    #living trees density

flora_dbmean_livtrees <- flora_dbmean
flora_dbmean_livtrees <- flora_dbmean_livtrees %>%
  mutate(liv.trees.no = sum.trees - dead.tree.no) %>%
  mutate(liv.trees.dens = liv.trees.no/491*10000)
flora_dbmean_livtrees$habitat <- relevel(flora_dbmean_livtrees$habitat, c("unaffected"))
gam6 <- gam(liv.trees.dens ~ habitat, data = flora_dbmean_livtrees)

    #Summary tables

model_summaries <- list()
for (i in c(1:3,5,6)) {
gamno <- get(paste0("gam", i))
output <- capture.output(summary(gamno), file=NULL,append=FALSE) 
output[1] <- paste0(gamno$formula[2])
plot.new()
grid.table(output)
model_sum <- recordPlot()
model_summaries[[i]] <- model_sum
}


    #Plots are printed here, then manually saved for further processing in CorelDRAW

emmeans_plot(gam1,"total tree density (per ha)","total tree density (per ha)")
emmeans_plot(gam2,"stand age","stand age")
emmeans_plot(gam3,"dead trees density (per ha)","dead trees density (per ha)")
emmeans_plot(gam5, "DBH","DBH")
emmeans_plot(gam6,"density of living trees (per ha)","density of living trees (per ha)")

#significance levels of the results

gam1_pred <- ggeffects::predict_response(gam1,terms=c("habitat"))
ggeffects::test_predictions(gam1_pred)
total_tree_dens_pred <- as.data.frame(test_predictions(predict_response(gam1,terms=c("habitat"))))
total_tree_dens_pred$name <- "total tree density"
gam2_pred <- predict_response(gam2,terms=c("habitat"))
stand.age_pred <- as.data.frame(test_predictions(predict_response(gam2,terms=c("habitat"))))
stand.age_pred$name <- "stand.age"
gam3_pred <- predict_response(gam3,terms=c("habitat"))
dead.trees.dens_pred <- as.data.frame(test_predictions(predict_response(gam3,terms=c("habitat"))))
dead.trees.dens_pred$name <- "dead.trees.dens"
gam5_pred <- predict_response(gam5,terms=c("habitat"))
dbh_pred <- as.data.frame(test_predictions(predict_response(gam5,terms=c("habitat"),margin="empirical")))
dbh_pred$name <- "DBH"
gam6_pred <- predict_response(gam6,terms=c("habitat"))
liv.trees.dens_pred <- as.data.frame(test_predictions(predict_response(gam6,terms=c("habitat"))))
liv.trees.dens_pred$name <- "liv.trees.dens"

total_table <- rbind(total_tree_dens_pred,liv.trees.dens_pred, dead.trees.dens_pred,stand.age_pred,dbh_pred)


#tree composition differences

flora_dbmean_perm <- flora_dbmean

flora_dbmean_perm <- flora_dbmean_perm %>%
  mutate(across(6:14, ~ . / sum.trees, .names = "{.col}_share"))

        community_data <- flora_dbmean_perm[, c(75:83)]
        env_data <- flora_dbmean_perm[, c(5,62)]
set.seed(1)
pairwise.adonis2(community_data ~ habitat,
                 data=flora_dbmean_perm,
                 sim.function = "vegdist",
                 sim.method = "bray",
                 p.adjust.m = "bonferroni",
                 perm = 9999)     


#### GAMMs for TreM richnesses ####

trems_richness_db <- flora_dbmean


trems_richness_db$habitat <- relevel(trems_richness_db$habitat, c("unaffected"))
emms_total_rich <- data.frame()
table_list_rich <- list()
plot_list_rich <- list()

for(i in c(63,66:72)) { 
  g1<-gam(trems_richness_db[,i] ~  
            habitat, data=trems_richness_db, family=poisson, method="ML")
  g2<- gam(trems_richness_db[,i] ~ 
             habitat + stand.age, data=trems_richness_db, family=poisson, method="ML")
  aicnostand <- round(AIC(logLik(g1)), 2)
  aicstand <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
  } else {
    g_good <- g1 }
  summary(g_good)
  colnames(trems_richness_db[i])
  model_emms <- predict_response(g_good,terms=c("habitat"))
  emms2 <- as.data.frame(test_predictions(model_emms))
  emms2$name <- colnames(trems_richness_db[i])
  emms_total_rich <- rbind(emms_total_rich, emms2)
  plotplease <- ggplot(model_emms, aes(x, predicted, col=x, fill=x)) +
    geom_point(cex=3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width=0.2,linewidth=1) +
    theme_classic() +
    scale_color_manual(values=c("#81b29a","#e07a5f","#3d405b")) +
    xlab("habitat type") +
    ylab("response") +
    ggtitle(colnames(trems_richness_db[i])) 
  j = i - 62
  plot_list_rich[[j]] <- plotplease
  output<-capture.output(summary(g_good), file=NULL,append=FALSE) 
  output[1] <- colnames(trems_richness_db[i])
  output[4] <- paste("AIC no.stand.age", aicnostand, "    ", "AIC stand.age", aicstand)
  plot.new()
  grid.table(output)
  model_sum <- recordPlot()
  table_list_rich[[j]] <- model_sum
}

plot_list_rich <- plot_list_rich[!sapply(plot_list_rich,is.null)]
plot_list_rich2 <- plot_list_rich
glist <- lapply(plot_list_rich2, ggplotGrob)
ggsave("plots__richnesses.pdf", marrangeGrob(glist, nrow = 4, ncol = 2), height=9, width=9, dpi = 300)



#### GAMMs - Abundance of individual TREMS with total sum >30 ####


trems_db_individual <- flora_dbmean

modelledspecies <- c()
for(i in 15:61) { if(is.numeric(trems_db_individual[,i]) == TRUE) { if(sum((trems_db_individual[,i]), na.rm=TRUE) > 30){
  modelledspecies <- c(modelledspecies, colnames(trems_db_individual)[i]) }}}
trems_db_individual3 <- cbind(trems_db_individual[,c(1:14,62)], filter(trems_db_individual[,modelledspecies]))

trems_db_individual3$habitat <- relevel(trems_db_individual3$habitat, c("unaffected"))
emms_total_individual <- data.frame()
table_list_individual <- list()
plot_list_individual <- list()
set.seed(1)
for(i in c(16:31)) { 
  g1<-gam(trems_db_individual3[,i] ~  
            habitat, data=trems_db_individual3, family=nb, method="ML")
  g2<- gam(trems_db_individual3[,i] ~ 
             habitat + stand.age, data=trems_db_individual3, family=nb, method="ML")
  aicnostand <- round(AIC(logLik(g1)), 2)
  aicstand <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
  } else {
    g_good <- g1 }
  summary(g_good)
  colnames(trems_db_individual3[i])
  model_emms <- predict_response(g_good,terms=c("habitat"))
  emms2 <- as.data.frame(test_predictions(model_emms))
  emms2$name <- colnames(trems_db_individual3[i])
  emms_total_individual <- rbind(emms_total_individual, emms2)
  plotplease <- ggplot(model_emms, aes(x, predicted, col=x, fill=x)) +
    geom_point(cex=3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width=0.2,size=1) +
    theme_classic() +
    scale_color_manual(values=c("#81b29a","#e07a5f","#3d405b")) +
    xlab("habitat type") +
    ylab("response") +
    ggtitle(colnames(trems_db_individual3[i])) 
  j = i - 15
  plot_list_individual[[j]] <- plotplease
  output<-capture.output(summary(g_good), file=NULL,append=FALSE) 
  output[1] <- colnames(trems_db_individual3[i])
  output[4] <- paste("AIC no.stand.age", aicnostand, "    ", "AIC stand.age", aicstand)
  plot.new()
  grid.table(output)
  model_sum <- recordPlot()
  table_list_individual[[j]] <- model_sum
  
}

plot_list_individual <- plot_list_individual[!sapply(plot_list_individual,is.null)]
plot_list_individual2 <- plot_list_individual
glist <- lapply(plot_list_individual2, ggplotGrob)
ggsave("plots_singletrems_all.pdf", marrangeGrob(glist, nrow = 6, ncol = 3), height=14, width=18, dpi = 300) 

dev.off()
plot_list_individual2_sign <- plot_list_individual2[c(1:10,14,15,16)]
plot_list_individual2_sign2 <- lapply(plot_list_individual2_sign, ggplotGrob)
ggsave("singletrems_significant.pdf", marrangeGrob(plot_list_individual2_sign2, nrow = 6, ncol = 3), height=12, width=12, dpi = 300)


#1 pdf of all model summaries:

summaries_total <- list()
model_summaries <- model_summaries[!sapply(model_summaries,is.null)]
summaries_total <- c(model_summaries, table_list_rich, table_list_individual)
summaries_total <- summaries_total[!sapply(summaries_total,is.null)]

title_page <- list("spirittheyaregonespirittheyvevanishedisthebestalbumever")
summaries_total <- c(title_page, summaries_total)

pdf("Appendix B.pdf", width = 14, height = 10)

for (i in 1:length(summaries_total)) {
  if (i == 1) {
    # Special formatting for the title page
    plot.new()  # Create a blank page
    
    # Large title for "Appendix B"
    text(0.5, 0.7, "Appendix B", cex = 3, font = 2)  # Bigger, bold text
    
    # Smaller text for the rest of the title
    text(0.5, 0.55, "Outputs of all GAM/GAMM models:", cex = 1.5, font = 1)
    text(0.5, 0.48, "(Pages 1-5) Stand structure", cex = 1.5, font = 1)
    text(0.5, 0.42, "(Pages 6-13) TreM richnesses", cex = 1.5, font = 1)
    text(0.5, 0.36, "(14-29) Single TreMs", cex = 1.5, font = 1)
  } else {
    # Regular model summary printing
    print(summaries_total[[i]])
  
  
  # Add page number slightly below the top
    k <- i - 1
  mtext(paste("Page", k), side = 3, line = -2, adj = 0.49, cex = 1.5)
}}

dev.off()

total_emmeans_table <- bind_rows(total_table, emms_total_rich, emms_total_individual)
total_emmeans_table <- total_emmeans_table[,c("name",setdiff(colnames(total_emmeans_table), "name"))]
total_emmeans_table <- total_emmeans_table %>% relocate(z, .before = t)
total_emmeans_table[is.na(total_emmeans_table)] <- "-"
write.xlsx(total_emmeans_table, "total_emmeans_table.xlsx", overwrite = TRUE)
