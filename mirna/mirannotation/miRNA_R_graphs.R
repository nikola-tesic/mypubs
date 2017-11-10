library(knitr)
library(rmarkdown)
library(ggplot2)
library(dplyr)

data <- read.table("/Users/vrjcc/Documents/git_repo/mypubs/mirna/mirannotation/results.txt", skip = 1)
data$V2 <- as.character(data$V2)
data$V4 <- as.character(data$V4)
data$V5 <- data$V5/abs(data$V5)
data$V5[data$V5 == 1] <- "Na"
data$V5[data$V5 == -1] <- "Add"
data$V6 <- data$V6/abs(data$V6)
data$V6[data$V6 == 1] <- "Na"
data$V6[data$V6 == -1] <- "Mut"
data$changes <- paste(data$V5, data$V6)
data$V9 <- as.character(data$V9)
data$V10 <- as.character(data$V10)

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
data$V10 <- lapply(data$V10,function(x){
  k1 = unlist(strsplit(x,'_'))
  if (k1[1] == 'razer3') {
    k3 = unlist(strsplit(k1[4],''))
    k5 = c(k1[1],k1[2],k3[1])
    v = paste(k5,collapse = '_')
  }
  else if (k1[1] == 'star') {
    k3 = unlist(strsplit(k1[3],''))
    k5 = c(k1[1],k1[2],k3[1])
    v = paste(k5,collapse = '_')
  }
  else{
    v = x
  }
  return (v)
})
data$groupup <-lapply(data$V10,function(x){
  k1 = unlist(strsplit(x,'_'))
  k2 = c(k1[2],k1[3])
  v = paste(k2,collapse = '_')
  return(v)
})
data$groupup <- as.character(data$groupup)
data$V10 <- as.character(data$V10)

colnames(data) <- c("name","miRNA","mapped","hairpin","add","mut","read_length","quality","aligned","tool","changes","TP","TPmirna",'p3','p5','trim','groupup')

# subsets from data

dt_3p_1 = subset(data,data[17]=='3p_1')
dt_3p_2 = subset(data,data[17]=='3p_2')
dt_5p_1 = subset(data,data[17]=='5p_1')
dt_5p_2 = subset(data,data[17]=='5p_2')

# Mapped

mapped_yes_no_plot <- function(dt_sub){
  
  dt = dt_sub %>% group_by(tool, mapped) %>% summarise(total = n()) %>% as_data_frame()
  #
  str_name = deparse(substitute(dt_sub))
  if (str_name != 'data'){ 
  str_name = paste(unlist(strsplit(deparse(substitute(dt_sub)),"_"))[2:3],collapse="_")
  }
  #
  ggplot(dt, aes(x = tool, y = total, fill = mapped)) + geom_bar(stat = "identity") + 
    ggtitle("Mapped yes/no : ",str_name)+ geom_text(aes(label = total), vjust = -1) + theme_bw() + labs(x = "") + 
    ylim(0, max(dt$total) + 2000) + scale_fill_brewer("mapped", palette = "Set1") + 
    facet_wrap(~mapped, scales = "free_y", ncol = 1) + theme(axis.text.x = element_text(angle = 90))

}
mapped_yes_no_plot(data)
#
mapped_yes_no_plot(dt_3p_1)
mapped_yes_no_plot(dt_3p_2)
mapped_yes_no_plot(dt_5p_1)
mapped_yes_no_plot(dt_5p_2)

# Mapped - Box Plot

mapped_box_plot <- function(dt_sub){
  
  #
  str_name = deparse(substitute(dt_sub))
  if (str_name != 'data'){ 
    str_name = paste(unlist(strsplit(deparse(substitute(dt_sub)),"_"))[2:3],collapse="_")
  }
  #
  ggplot(dt_sub, aes(tool, read_length, fill = mapped)) + geom_boxplot() + ggtitle('Mapped box plot : ',str_name)+ theme_bw() +
    labs(x = "", y = "length") +
    scale_fill_brewer("mapped", palette = "Set1") + theme(axis.text.x = element_text(angle = 90))
}
mapped_box_plot(data)
#
mapped_box_plot(dt_3p_1)
mapped_box_plot(dt_3p_2)
mapped_box_plot(dt_5p_1)
mapped_box_plot(dt_5p_2)

# Changes
changes_plot <- function(dt_sub){
  
  #
  str_name = deparse(substitute(dt_sub))
  if (str_name != 'data'){ 
    str_name = paste(unlist(strsplit(deparse(substitute(dt_sub)),"_"))[2:3],collapse="_")
  }
  #
ggplot(dt_sub, aes(tool, fill = changes)) + geom_bar() + ggtitle('Changes : ',str_name)+ theme_bw() + theme_bw() + labs(x = "") + 
  facet_wrap(~mapped) + scale_fill_brewer("changes", palette = "Set1") + theme(axis.text.x = element_text(angle = 90))
}
changes_plot(data)
#
changes_plot(dt_3p_1)
changes_plot(dt_3p_2)
changes_plot(dt_5p_1)
changes_plot(dt_5p_2)

#Specifity at precursor level - How many were assigned to the correct miRNA using the precursor name
#as the true positive. Red would be ???not correct??? and blue ???correct???. This is only considering mapped sequences.

spec_precursor_plot <- function(data){
  
  #
  str_name = deparse(substitute(data))
  if (str_name != 'data'){ 
    str_name = paste(unlist(strsplit(deparse(substitute(data)),"_"))[2:3],collapse="_")
  }
  #
  
  dt = data %>% filter(mapped == "yes") %>% group_by(tool, TP) %>% summarise(total = n()) %>% as_data_frame()

  ggplot(dt, aes(x = tool, y = total, fill = factor(TP))) + geom_bar(stat = "identity") +
    geom_text(aes(label = total), vjust = -1, size = 3) + theme_bw() + labs(x = "") + 
    ylim(0, max(dt$total) + 200) + scale_fill_brewer(guide = FALSE, "correct", 
    palette = "Set1") + theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~TP, scales = "free_y") + ggtitle("Specificity at precursor level : ", str_name)
}

spec_precursor_plot(data)
#
spec_precursor_plot(dt_3p_1)
spec_precursor_plot(dt_3p_2)
spec_precursor_plot(dt_5p_1)
spec_precursor_plot(dt_5p_2)

# Specifity at mature miRNA level - Same logic than before but with the miRNA names in case the tool gives
#the miRNA names, if not, only the three first field in the name are used as annotation 
#(i.e hsa-let-7a-1 will ignore any character beyond 7a). This will increase the number of TP,
#since many miRNAs has multiple precursors being the same mature miRNA at the end.

mature_mirna_level_plot <- function(data){
  #
  str_name = deparse(substitute(data))
  if (str_name != 'data'){ 
    str_name = paste(unlist(strsplit(deparse(substitute(data)),"_"))[2:3],collapse="_")
  }
  #
  dt = data %>% filter(mapped == "yes") %>% group_by(tool, TPmirna) %>% summarise(total = n()) %>% as_data_frame()
  
  ggplot(dt, aes(x = tool, y = total, fill = factor(TPmirna))) + geom_bar(stat = "identity") + 
    geom_text(aes(label = total), vjust = -1, size = 3) + theme_bw() + labs(x = "") + 
    ylim(0, max(dt$total) + 200) + scale_fill_brewer(guide = FALSE, "correct", 
    palette = "Set1") + theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~TPmirna, 
    scales = "free_y") + ggtitle("Specificity at mature miRNA level : ",str_name)
}

mature_mirna_level_plot(data)
#
mature_mirna_level_plot(dt_3p_1)
mature_mirna_level_plot(dt_3p_2)
mature_mirna_level_plot(dt_5p_1)
mature_mirna_level_plot(dt_5p_2)


# Aligned yes/no  --- Proportion of mapped and no-mapped sequences

aligned_yes_no <- function(data){
  
  str_name = deparse(substitute(data))
  if (str_name != 'data'){ 
    str_name = paste(unlist(strsplit(deparse(substitute(data)),"_"))[2:3],collapse="_")
  }
  
  dt = subset(data,aligned %in% c('yes','no')) %>% group_by(tool, aligned) %>% summarise(total = n()) %>% as_data_frame()

  ggplot(dt, aes(x = tool, y = total, fill = aligned)) + geom_bar(stat = "identity") + 
    geom_text(aes(label = total), vjust = -1) + theme_bw() + labs(x = "") + 
    ylim(0, max(dt$total) + 2000) + scale_fill_brewer("aligned", palette = "Set1") + 
    facet_wrap(~aligned, scales = "free_y", ncol = 1) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Aligned : ", str_name)
}

aligned_yes_no(data)
#
aligned_yes_no(dt_3p_1)
aligned_yes_no(dt_3p_2)
aligned_yes_no(dt_5p_1)
aligned_yes_no(dt_5p_2)
