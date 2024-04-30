library("dplyr")
library(TMExplorer)
library(caTools)
library(factoextra)
library(ROCR)

res = TMExplorer::queryTME(geo_accession = "GSE131907", sparse = TRUE)
data = counts(res[[1]])
assign <- as.data.frame(cbind(colnames(data), res[[1]]$label))
rm(res)
rownames(assign) <- assign$V1
assign <- subset(assign, select = -c(V1))
assign <- dplyr::filter(assign, grepl("ndothelial", assign$V2))

data <- data[,grepl("LUNG", colnames(data))]
data <- data[, colnames(data) %in% rownames(assign)]

lungs_normal = t(data[,grepl("LUNG_N", colnames(data))])
lungs_disease = t((data[,grepl("LUNG_T", colnames(data))]))
rm(data)
lungs_total <- rbind(lungs_normal, lungs_disease)
if (length(lungs_total[,1]) > 1000) {
  lungs_total <- lungs_total[sample((1:length(lungs_total[,1])), size = 1000, replace = F), ]
}
vars <- apply(lungs_total, MARGIN = 2, FUN = var)
sorted_variances <- sort(vars, decreasing = TRUE, index.return = TRUE)$ix[1:2000]

lungs_normal <- as.matrix(lungs_normal[, sorted_variances])
lungs_normal <- as.data.frame(lungs_normal)
lungs_normal$state <- "Normal"
lungs_disease <- as.matrix(lungs_disease[, sorted_variances])
lungs_disease <- as.data.frame(lungs_disease)
lungs_disease$state <- "Tumor"
rm(lungs_total)

lungs <- rbind(lungs_normal, lungs_disease)
lungs <- merge(lungs, assign, by = 0)
rownames(lungs) <- lungs$Row.names
lungs <- subset(lungs, select = -c(Row.names))

if(length(lungs[ , 1]) > 2000) {
  lungs_sub <- dplyr::sample_n(lungs, 2000)
} else {
  lungs_sub <- lungs
}

pc_lung_sub <- prcomp(subset(lungs_sub, select = -c(state, V2)), scale = TRUE, center = TRUE)
fviz_pca_ind(pc_lung_sub, label = "none", habillage = lungs_sub$state)
pc_results <- get_pca_ind(pc_lung_sub)$coord
pc_results <- as.data.frame(pc_results)
pc_results$state <- lungs_sub$state

top <- as.data.frame(get_pca_var(pc_lung_sub)$contrib)
top <- top[order(top$Dim.1, decreasing = TRUE),]
print(rownames(top)[1:5])

pc_results <- pc_results %>% mutate(
  state = case_when(
    state == "Tumor" ~ 1,
    state == "Normal" ~ 0
  )
)

split <- sample.split(pc_results, SplitRatio = 0.8)
train_reg <- subset(pc_results, split == "TRUE")
test_reg <- subset(pc_results, split == "FALSE")

logistic_model <- glm(state ~ Dim.1 + Dim.2 + Dim.3 + Dim.4 + Dim.5 + Dim.6 + Dim.7 + Dim.8 + Dim.9 + Dim.10,
                      data = train_reg,
                      family = "binomial")
summary(logistic_model)

PredLR <- predict(logistic_model, test_reg, type = "response")
lgPredObj <- prediction(PredLR, test_reg$state)
lgPerfObj <- performance(lgPredObj, "tpr","fpr")
AUC <- performance(lgPredObj, measure = "auc")
AUC <- AUC@y.values[[1]]
print(AUC)
plot(lgPerfObj, main = "ROC Curve", col = 2, lwd = 2)
abline(a = 0,b = 1,lwd = 2,lty = 3,col = "black")

