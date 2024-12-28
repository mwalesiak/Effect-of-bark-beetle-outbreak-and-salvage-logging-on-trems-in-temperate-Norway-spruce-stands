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

# All resulting figures were later improved using corelDRAW.
#generated csv's were later expanded in excel to create the appendices. 
#PDF files were glued and title pages added.

#setWD

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#preparing variables

flora_dbmean <- read.table("flora_plotlevel_final.csv", header = TRUE, sep = ";", dec = ".", na.strings="NA", strip.white=TRUE)
trem_cat <- read.table("trems_names_explained.csv", header = TRUE, sep = ";", dec = ".", na.strings="NA", strip.white=TRUE)
View(flora_dbmean)
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

View(flora_dbmean)


#### Community differences between the plots ####

emmeans_plot <- function(model, axisy, title) {
  model_emms <- predict_response(model,terms=c("habitat"))
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
flora_dbmean <- flora_dbmean %>% mutate(sum.trees.dens = sum.trees/491*10000)
flora_dbmean <- flora_dbmean %>% mutate(dead.trees.dens = dead.tree.no/491*10000)
gam1 <- gam(sum.trees.dens ~ habitat, data = flora_dbmean)
gam2 <- gam(stand.age ~ habitat, data = flora_dbmean)
gam3 <- gam(dead.trees.dens ~ habitat, data = flora_dbmean)
gam4 <- gam(dead.tree.share ~ habitat, data = flora_dbmean)


    #circumference (single trees level)
set.seed(1)
flora_db_circ <- read.table("trems-circumferences.csv", header = TRUE, sep = ";", dec = ".", na.strings="NA", strip.white=TRUE)
flora_db_circ$habitat <- relevel(as.factor(flora_db_circ$habitat), c("unaffected"))
flora_db_circ$plot_ID <- as.factor(flora_db_circ$plot_ID)
gam5 <- gam(circumference ~ habitat + s(plot_ID,bs="re"), data = flora_db_circ)

    #living trees density

flora_dbmean_livtrees <- flora_dbmean
flora_dbmean_livtrees <- flora_dbmean_livtrees %>%
  mutate(liv.trees.no = sum.trees - dead.tree.no) %>%
  mutate(liv.trees.dens = liv.trees.no/491*10000)
View(flora_dbmean_livtrees)
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

model_summaries <- model_summaries[!sapply(model_summaries,is.null)]

dev.off()
pdf("model_outputs-structure.pdf")
pdf.options(width = 14, height = 10)
for (i in 1:length(model_summaries)){
  print(model_summaries[[i]])
}
dev.off()

    #Plots are printed here, then manually saved for further processing in CorelDRAW

emmeans_plot(gam1,"total tree density (per ha)","total tree density (per ha)")
emmeans_plot(gam2,"stand age","stand age")
emmeans_plot(gam3,"dead trees density (per ha)","dead trees density (per ha)")
emmeans_plot(gam5, "tree circumference","tree circumference")
emmeans_plot(gam6,"density of living trees (per ha)","density of living trees (per ha)")

#significance levels of the results

gam1_pred <- predict_response(gam1,terms=c("habitat"))
total_tree_dens_pred <- as.data.frame(test_predictions(predict_response(gam1,terms=c("habitat"))))
total_tree_dens_pred$name <- "total tree density"
gam2_pred <- predict_response(gam2,terms=c("habitat"))
stand.age_pred <- as.data.frame(test_predictions(predict_response(gam2,terms=c("habitat"))))
stand.age_pred$name <- "stand.age"
gam3_pred <- predict_response(gam3,terms=c("habitat"))
dead.trees.dens_pred <- as.data.frame(test_predictions(predict_response(gam3,terms=c("habitat"))))
dead.trees.dens_pred$name <- "dead.trees.dens"
gam5_pred <- predict_response(gam5,terms=c("habitat"))
circumference_pred <- as.data.frame(test_predictions(predict_response(gam5,terms=c("habitat"))))
circumference_pred$name <- "circumference"
gam6_pred <- predict_response(gam6,terms=c("habitat"))
liv.trees.dens_pred <- as.data.frame(test_predictions(predict_response(gam6,terms=c("habitat"))))
liv.trees.dens_pred$name <- "liv.trees.dens"

total_table <- rbind(total_tree_dens_pred,liv.trees.dens_pred, dead.trees.dens_pred,stand.age_pred,circumference_pred)

write.csv(total_table, "structure_differences.csv")


#### GAMMs for TreM richnesses ####

flora_dbmean2 <- flora_dbmean
View(flora_dbmean2)

flora_dbmean2$habitat <- relevel(flora_dbmean2$habitat, c("unaffected"))
emms_total <- data.frame()
table_list <- list()
plot_list <- list()

for(i in c(63,66:72)) { 
  g1<-gam(flora_dbmean2[,i] ~  
            habitat, data=flora_dbmean2, family=poisson, method="ML")
  g2<- gam(flora_dbmean2[,i] ~ 
             habitat + stand.age, data=flora_dbmean2, family=poisson, method="ML")
  aicnostand <- round(AIC(logLik(g1)), 2)
  aicstand <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
  } else {
    g_good <- g1 }
  summary(g_good)
  colnames(flora_dbmean2[i])
  model_emms <- predict_response(g_good,terms=c("habitat"))
  emms2 <- as.data.frame(test_predictions(model_emms))
  emms2$name <- colnames(flora_dbmean2[i])
  emms_total <- rbind(emms_total, emms2)
  plotplease <- ggplot(model_emms, aes(x, predicted, col=x, fill=x)) +
    geom_point(cex=3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width=0.2,linewidth=1) +
    theme_classic() +
    scale_color_manual(values=c("#81b29a","#e07a5f","#3d405b")) +
    xlab("habitat type") +
    ylab("response") +
    ggtitle(colnames(flora_dbmean2[i])) 
  j = i - 62
  plot_list[[j]] <- plotplease
  output<-capture.output(summary(g_good), file=NULL,append=FALSE) 
  output[1] <- colnames(flora_dbmean2[i])
  output[4] <- paste("AIC habitatsolo", aicnostand, "    ", "AIC stand.age", aicstand)
  plot.new()
  grid.table(output)
  model_sum <- recordPlot()
  table_list[[j]] <- model_sum
}

dev.off()
pdf("model_outputs_richnesses.pdf")
pdf.options(width = 14, height = 10)
for (i in 1:length(table_list)){
  print(table_list[[i]])
}
dev.off()

write.csv(emms_total, "pairwise_results_richnesses.csv")

plot_list <- plot_list[!sapply(plot_list,is.null)]
plot_list2 <- plot_list
glist <- lapply(plot_list2, ggplotGrob)
ggsave("plots__richnesses.pdf", marrangeGrob(glist, nrow = 4, ncol = 2), height=9, width=9, dpi = 300)

#### GAMMs - Abundance of individual TREMS with total sum >30 ####
View(flora_dbmean)
modelledspecies <- c()
for(i in 15:61) { if(is.numeric(flora_dbmean[,i]) == TRUE) { if(sum((flora_dbmean[,i]), na.rm=TRUE) > 30){
  modelledspecies <- c(modelledspecies, colnames(flora_dbmean)[i]) }}}
flora_dbmean3 <- cbind(flora_dbmean[,c(1:14,62)], filter(flora_dbmean[,modelledspecies]))
View(flora_dbmean3)

flora_dbmean3$habitat <- relevel(flora_dbmean3$habitat, c("unaffected"))
emms_total <- data.frame()
table_list <- list()
plot_list <- list()

for(i in c(16:31)) { 
  g1<-gam(flora_dbmean3[,i] ~  
            habitat, data=flora_dbmean3, family=nb, method="ML")
  g2<- gam(flora_dbmean3[,i] ~ 
             habitat + stand.age, data=flora_dbmean3, family=nb, method="ML")
  aicnostand <- round(AIC(logLik(g1)), 2)
  aicstand <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
  } else {
    g_good <- g1 }
  summary(g_good)
  colnames(flora_dbmean3[i])
  model_emms <- predict_response(g_good,terms=c("habitat"))
  emms2 <- as.data.frame(test_predictions(model_emms))
  emms2$name <- colnames(flora_dbmean3[i])
  emms_total <- rbind(emms_total, emms2)
  plotplease <- ggplot(model_emms, aes(x, predicted, col=x, fill=x)) +
    geom_point(cex=3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width=0.2,size=1) +
    theme_classic() +
    scale_color_manual(values=c("#81b29a","#e07a5f","#3d405b")) +
    xlab("habitat type") +
    ylab("response") +
    ggtitle(colnames(flora_dbmean3[i])) 
  j = i - 15
  plot_list[[j]] <- plotplease
  output<-capture.output(summary(g_good), file=NULL,append=FALSE) 
  output[1] <- colnames(flora_dbmean3[i])
  output[4] <- paste("AIC habitatsolo", aicnostand, "    ", "AIC stand.age", aicstand)
  plot.new()
  grid.table(output)
  model_sum <- recordPlot()
  table_list[[j]] <- model_sum
}

dev.off()
pdf("model_outputs_singletrems.pdf")
pdf.options(width = 14, height = 10)
for (i in 1:length(table_list)){
  print(table_list[[i]])
}
dev.off()

write.csv(emms_total, "pairwise_results_singletrems.csv")

plot_list <- plot_list[!sapply(plot_list,is.null)]
plot_list2 <- plot_list
glist <- lapply(plot_list2, ggplotGrob)
ggsave("plots_singletrems_all.pdf", marrangeGrob(glist, nrow = 6, ncol = 3), height=14, width=18, dpi = 300) 

dev.off()
plot_list2_sign <- plot_list2[c(1:11,14,15)]
plot_list2_sign2 <- lapply(plot_list2_sign, ggplotGrob)
ggsave("singletrems_significant.pdf", marrangeGrob(plot_list2_sign2, nrow = 6, ncol = 3), height=12, width=12, dpi = 300)

dev.off()

