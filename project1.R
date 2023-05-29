#CS461
#title: Project 1
#author: Grant Davis
#date: 10/19/2022

knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(tidyverse)

#(1)
# Loading data using function `load()`
load("GSE9782.RData")

treated <- which(group == "treated")
control <- which(group != "treated")

#(2)
# log base 2 of the expression
log2Data <- log(x = dataGSE9782, base = 2)

# a function to calculate the difference in mean
calc_mean_diff <- function(x, treated, control) {
  mean(x[treated]) - mean(x[control])
}
# used apply to call the function
logFC <- apply(log2Data, 1, calc_mean_diff, treated, control)

#(3)
# Function to calculate p-value
calc_p_value <- function(x, treated, control) {
  t.test(x[treated], x[control])$p.value
}
# Function to calculate t-score
calc_t_score <- function(x, treated, control) {
  t.test(x[treated], x[control])$statistic
}

PValue <- apply(log2Data, MARGIN = 1, FUN = calc_p_value, treated, control)
TScore <- apply(log2Data, MARGIN = 1, FUN = calc_t_score, treated, control)

# rownames used as gene ids
geneIds <- rownames(log2Data)

df <- data.frame(
  row.names = NULL,
  "GeneID" = geneIds,
  "PValue" = PValue,
  "TScore" = TScore,
  "LogFC" = logFC
)

# printing the head of data frame
head(df)
save(df, file = "DE-results.rds")

#(4)
pdf(file = "volcano.pdf",
    width = 10,
    height = 10)
plot(
  x = df$LogFC,
  y = -log10(df$PValue),
  xlab = 'logFC',
  ylab = '-log10(p-value)',
  main = "Volcano plot",
  col = ifelse(abs(df$LogFC) > 1 & df$PValue < 0.05, 'red', 'black'),
  xlim = c(-2, 2)
)
abline(h = -log10(0.05), col = "red")
abline(v = -1, col = "blue")
abline(v = 1, col = "blue ")
dev.off()

#(5)
# Function to calculate p-value with Wilcoxon test
calc_p_value_wilcox <- function(x, treated, control) {
  wilcox.test(x[treated], x[control])$p.value
}

# Function to calculate t-score using Wilcoxon test
calc_t_score_wilcox <- function(x, treated, control) {
  wilcox.test(x[treated], x[control])$statistic
}

# used apply to call the function
PValueWx <-
  apply(log2Data, MARGIN = 1, FUN = calc_p_value_wilcox, treated, control)
TScoreWx <-
  apply(log2Data, MARGIN = 1, FUN = calc_t_score_wilcox, treated, control)

# Finding the genes that are significant at 5% based on both types of tests (t-test and Wilcoxon)
dfPvalue <-
  data.frame(row.names = geneIds,
             "PValueWx" = PValueWx,
             "PValueTs" = PValue)

significantGenes <-
  row.names(dfPvalue[which(dfPvalue$PValueWx < 0.05 &
                             dfPvalue$PValueTs < 0.05),])

# Printing the length of significant genes
length(significantGenes)

# Printing the first 10 significant genes
print(significantGenes[1:10])
save(significantGenes, file = "Wilcoxon-results.rds")

#(6)
pdf(file = "hist.pdf",
    width = 10,
    height = 10)
hist(log2Data, probability = T)
dev.off()

#(7)
pdf(file = "boxplot.pdf",
    width = 10,
    height = 10)
maxGene_df <- data.frame(group = group,
                         value = log2Data[which(abs(logFC) == max(abs(logFC))),])
maxGene_df %>% ggplot(aes(x = group, y = value)) + geom_boxplot() + theme_classic()
dev.off()

