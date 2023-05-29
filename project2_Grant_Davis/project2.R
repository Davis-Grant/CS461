#CS461
#title: Project 2
#author: Grant Davis
#date: 12/09/2022

knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
library(ggplot2)
library(tidyverse)

library(GEOquery)
library(stringr)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)


#Question 1
dataset <- "GSE19804"
gsets <- getGEO(dataset, GSEMatrix = TRUE, getGPL = TRUE)
gset <- gsets[[1]]
expr <- exprs(gset)

pdata <- pData(gset)
control <- rownames(pdata[grep("Lung Normal", pdata$title), ])
cancer <- rownames(pdata[grep("Lung Cancer", pdata$title), ])

 
# Question 2
calc_mean_diff <- function(x, cancer, control) {
  mean(x[cancer]) - mean(x[control])
}
# Function to calculate p-value
calc_p_value <- function(x, cancer, control) {
  t.test(x[cancer], x[control])$p.value
}
# Function to calculate t-score
calc_t_score <- function(x, cancer, control) {
  t.test(x[cancer], x[control])$statistic
}

# used apply to call the function
logFC <- apply(expr, 1, calc_mean_diff, cancer, control)
PValue <- apply(expr, MARGIN = 1, FUN = calc_p_value, cancer, control)
TScore <- apply(expr, MARGIN = 1, FUN = calc_t_score, cancer, control)

# rownames used as gene ids
geneIds <- rownames(expr)

df <- data.frame(
  row.names = NULL,
  "GeneID" = geneIds,
  "PValue" = PValue,
  "TScore" = TScore,
  "LogFC" = logFC
)

save(df, file = "DE-results.rds")

#Question 3
pdf(file = "volcano.pdf",
    width = 10,
    height = 10)
plot(
  x = df$LogFC,
  y = -log10(df$PValue),
  xlab = "logFC",
  ylab = "-log10(p-value)",
  main = "Volcano plot",
  col = ifelse(abs(df$LogFC) > 1 & df$PValue < 0.05, "red", "black")
)
abline(h = -log10(0.05), col = "red")
abline(v = -1, col = "blue")
abline(v = 1, col = "blue")
dev.off()

# Question 4

calc_t_score <- function(x, cancer, control) {
  t.test(x[cancer], x[control])$statistic
}

TScore <- apply(expr, MARGIN = 1, FUN = calc_t_score, cancer, control)
t_OBSERVED <- data.frame(TScore)

save(t_OBSERVED, file = "t-observed.rds")

# Question 5
permuteExpr <- expr
t_NULL_DISTRIBUTION <- lapply(c(1:10), function(i) {
  message(i)
  colnames(permuteExpr) <- sample(colnames(expr))
  apply(permuteExpr, MARGIN = 1, FUN = calc_t_score, cancer, control)
}) %>%
  do.call(what = cbind) %>%
  as.data.frame()

pT <- lapply(rownames(t_OBSERVED), function(r) {
  sum(abs(t_NULL_DISTRIBUTION[r, ]) > abs(t_OBSERVED[r, ])) / ncol(t_NULL_DISTRIBUTION)
}) %>% unlist()

save(pT, file = "p-empirical-t-score.rds")

# Question 6
calc_e_score <- function(x, cancer, control) {
  dist(rbind(mean(x[cancer]), mean(x[control])), method = "euclidean")
}

eSCORE <- apply(expr, MARGIN = 1, FUN = calc_e_score, cancer, control)
e_OBSERVED <- data.frame(eSCORE)

e_NULL_DISTRIBUTION <- lapply(c(1:10), function(i) {
  message(i)
  samples <- sample(ncol(expr), replace = F)
  controls <- samples[1:length(control)]
  cancers <- samples[(length(control) + 1):length(samples)]
  apply(
    expr,
    MARGIN = 1,
    FUN = calc_e_score,
    cancer = cancers,
    control = controls
  )
}) %>%
  do.call(what = cbind) %>%
  as.data.frame()

pE <- lapply(rownames(e_OBSERVED), function(r) {
  mean(abs(e_NULL_DISTRIBUTION[r, ]) > abs(e_OBSERVED[r, ]))
}) %>% unlist()
 
save(pE, file = "p-empirical-euclidean.rds")

# Question 7
pdf(file = "hist-pT.pdf",
    width = 10,
    height = 10)

par(mfrow = c(1, 2))
hist(pT, breaks = 10, col = "red", main = "pT")
dev.off()

pdf(file = "hist-pE.pdf",
    width = 10,
    height = 10)
hist(pE, breaks = 10, col = "blue", main = "pE")
dev.off()
 
# Question 8
cor(pT, pE)

