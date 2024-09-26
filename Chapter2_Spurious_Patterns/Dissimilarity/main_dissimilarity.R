source('Dissimilarity/dissimilarity_measure.R')

#germany
#load modified and raw gt data
gt_data_germany_final <- readRDS("Data/GT_Data/germany/gt_data_germany_final.RDS")
germany_modified1 <- readRDS("Data/GT_Data/modified/germany_modified1.RDS")
germany_modified2 <- readRDS("Data/GT_Data/modified/germany_modified2.RDS")

dissimilarity_germany_raw <- dissimilarity(gt_data_germany_final,12,T)
dissimilarity_germany_modified1 <- dissimilarity(germany_modified1,12,T)
dissimilarity_germany_modified2 <- dissimilarity(germany_modified2,12,T)

all_dissimilarities <- c(dissimilarity_germany_raw,dissimilarity_germany_modified1,
                         dissimilarity_germany_modified2)

plot(dissimilarity_germany_raw,ylim=c(min(all_dissimilarities),max(all_dissimilarities)),
     ylab='dissimilarity',xlab='',main='German GT data')
lines(dissimilarity_germany_modified1,col='red')
lines(dissimilarity_germany_modified2,col='blue')

#international
#load modified and raw gt data
gt_data_international_final <- readRDS("Data/GT_Data/international/gt_data_international_final.RDS")
international_modified1 <- readRDS("Data/GT_Data/modified/international_modified1.RDS")
international_modified2 <- readRDS("Data/GT_Data/modified/international_modified2.RDS")

dissimilarity_international_raw <- dissimilarity(gt_data_international_final,12,T)
dissimilarity_international_modified1 <- dissimilarity(international_modified1,12,T)
dissimilarity_international_modified2 <- dissimilarity(international_modified2,12,T)

all_dissimilarities <- c(dissimilarity_international_raw,dissimilarity_international_modified1,
                         dissimilarity_international_modified2)

plot(dissimilarity_international_raw,ylim=c(min(all_dissimilarities),max(all_dissimilarities)),
     ylab='dissimilarity',xlab='',main='Worldwide GT data')
lines(dissimilarity_international_modified1,col='red')
lines(dissimilarity_international_modified2,col='blue')





