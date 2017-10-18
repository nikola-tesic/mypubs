library(knitr)
library(rmarkdown)
library(ggplot2)
library(dplyr)
#
options(bitmapType = "cairo")
opts_chunk$set(tidy = TRUE, highlight = T,fig.align = "center", fig.height = 6, fig.width = 10, message = F, error = F, warning = F, bootstrap.show.code = FALSE)
#
data <- read.table("/Users/vrjcc/Desktop/results_sastanak.txt", skip = 1)
data$V2 <- as.character(data$V2)
data$V4 <- as.character(data$V4)
data$V5 <- data$V5/abs(data$V5)
data$V5[data$V5 == 1] <- "Na"
data$V5[data$V5 == -1] <- "Add"
data$V6 <- data$V6/abs(data$V6)
data$V6[data$V6 == 1] <- "Na"
data$V6[data$V6 == -1] <- "Mut"
data$changes <- paste(data$V5, data$V6)
#
data$V9 <- as.character(data$V9)
#
#######################################################
data$TP <- apply(data[, c(2, 4)], 1, function(x) {
  v <- grep(x[1], x[2], ignore.case = T)
  if (length(v) == 0) {
    v <- 0
  }
  return(v)
})
data$TPmirna <- apply(data[, c(2, 4)], 1, function(x) {
  h1 = unlist(strsplit(x[1], split = "-"))[1:3]
  h2 = unlist(strsplit(x[2], split = "-"))[1:3]
  v <- grep(paste0(h1, collapse = "-"), paste0(h2, collapse = "-"), ignore.case = T)
  if (length(v) == 0) {
    v <- 0
  }
  return(v)
})
#################################################
# 3p i 5p cuts

data$p3 <- lapply(data$V1, function(x) {
  y = unlist(strsplit(as.character(x), split = "_"))
  k = as.integer(unlist(strsplit(y[4], split = ":")))
  v <- k[1]
  return(v)
})

data$p5 <- lapply(data$V1, function(x) {
  y = unlist(strsplit(as.character(x), split = "_"))
  k = as.integer(unlist(as.integer(unlist(strsplit(y[4], split = ":")))))
  v <- k[2]
  return(v)
})

data$trim <- apply(data[,c(14,15)],1,function(x){
  if (x[1] !=0 && x[2] == 0 ){v<-'5p'}
  else if (x[1] ==0 && x[2] != 0){v<-'3p'}
  else if (x[1] !=0 && x[2] !=0) {v<-'both'}
  else {v<-'no change'}
  return (v)
})
#head(data)

###########################################
# Count of trim effects. Mind that we used the same file for every tool so we get the same numbers

dt2 = data %>% group_by(V10, trim) %>% summarise(total = n()) %>% as_data_frame()
ggplot(dt2, aes(x = V10, y = total, fill = trim)) + geom_bar(stat = "identity") + 
  geom_text(aes(label = total), vjust = -1) + theme_bw() + labs(x = "") + 
  ylim(0, max(dt2$total) + 2000) + scale_fill_brewer("trim end", palette = "Set1") + 
  facet_wrap(~trim, scales = "free_y", ncol = 1) + theme(axis.text.x = element_text(angle = 90))
##########################################
# Mapped yes/no and trim effects

ggplot(data, aes(V10, fill = trim)) + geom_bar() + theme_bw() + labs(x = "") + 
facet_wrap(~V3) + scale_fill_brewer("trim changes", palette = "Set1") + theme(axis.text.x = element_text(angle = 90))
##########################################
# Mapped yes/no  --- Proportion of mapped and no-mapped sequences

dt = data %>% group_by(V10, V3) %>% summarise(total = n()) %>% as_data_frame()
ggplot(dt, aes(x = V10, y = total, fill = V3)) + geom_bar(stat = "identity") + 
  geom_text(aes(label = total), vjust = -1) + theme_bw() + labs(x = "") + 
  ylim(0, max(dt$total) + 2000) + scale_fill_brewer("mapped", palette = "Set1") + 
  facet_wrap(~V3, scales = "free_y", ncol = 1) + theme(axis.text.x = element_text(angle = 90))
#########################################
# Boxplot with mapped reads

ggplot(data, aes(V10, V7, fill = V3)) + geom_boxplot() + theme_bw() + labs(x = "", 
y = "length") + scale_fill_brewer("mapped", palette = "Set1") + theme(axis.text.x = element_text(angle = 90))

########################################
# Changes

ggplot(data, aes(V10, fill = changes)) + geom_bar() + theme_bw() + labs(x = "") + 
facet_wrap(~V3) + scale_fill_brewer("changes", palette = "Set1") + theme(axis.text.x = element_text(angle = 90))

#######################################
# Specifity at precursor level - How many were assigned to the correct miRNA using the precursor name
#as the true positive. Red would be ???not correct??? and blue ???correct???. This is only considering mapped sequences.

dt = data %>% filter(V3 == "yes") %>% group_by(V10, TP) %>% summarise(total = n()) %>% 
  as_data_frame()
ggplot(dt, aes(x = V10, y = total, fill = factor(TP))) + geom_bar(stat = "identity") +
geom_text(aes(label = total), vjust = -1, size = 3) + theme_bw() + labs(x = "") + 
ylim(0, max(dt$total) + 200) + scale_fill_brewer(guide = FALSE, "correct", 
palette = "Set1") + theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~TP, 
scales = "free_y") + ggtitle("Specificity at precursor level")

#######################################
# Specifity at mature miRNA level - Same logic than before but with the miRNA names in case the tool gives
#the miRNA names, if not, only the three first field in the name are used as annotation 
#(i.e hsa-let-7a-1 will ignore any character beyond 7a). This will increase the number of TP,
#since many miRNAs has multiple precursors being the same mature miRNA at the end.

dt = data %>% filter(V3 == "yes") %>% group_by(V10, TPmirna) %>% summarise(total = n()) %>% 
  as_data_frame()
ggplot(dt, aes(x = V10, y = total, fill = factor(TPmirna))) + geom_bar(stat = "identity") + 
  geom_text(aes(label = total), vjust = -1, size = 3) + theme_bw() + labs(x = "") + 
  ylim(0, max(dt$total) + 200) + scale_fill_brewer(guide = FALSE, "correct", 
  palette = "Set1") + theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~TPmirna, 
  scales = "free_y") + ggtitle("Specificity at mature miRNA level")
########################################
# Aligned on hairpin

ggplot(subset(data, data$V9 %in% c('yes','no')), aes(V10, fill = trim)) + geom_bar() + theme_bw() + labs(x = "") + 
  facet_wrap(~V9) + scale_fill_brewer("Aligned", palette = "Set1") + theme(axis.text.x = element_text(angle = 90))

##########################################
# Aligned yes/no  --- Proportion of mapped and no-mapped sequences

dt = subset(data,data$V9 %in% c('yes','no')) %>% group_by(V10, V9) %>% summarise(total = n()) %>% as_data_frame()
ggplot(dt, aes(x = V10, y = total, fill = V9)) + geom_bar(stat = "identity") + 
  geom_text(aes(label = total), vjust = -1) + theme_bw() + labs(x = "") + 
  ylim(0, max(dt$total) + 2000) + scale_fill_brewer("aligned", palette = "Set1") + 
  facet_wrap(~V9, scales = "free_y", ncol = 1) + theme(axis.text.x = element_text(angle = 90))

########################################
