library(R.matlab)
library(tidyverse)

#load data
data <- readMat('M283BarcodeMatrix.mat')

#load raw count matrices (B1 and B2)
b1 = data$B1
b2 = data$B2

#specifically label columns
num_col = 1:41
colnames(b1) = num_col
colnames(b2) = num_col

#load in sorted column values
colVals = read.csv("fragileX_sorted_column_labels_withBrain.csv", header = FALSE)

brain1_info = colVals[colVals$V1 == 1,]
brain2_info = colVals[colVals$V1 == 2,]

#Brain 1
brain1 = b1[,brain1_info$V3] #fully populated, as expected

#Brain 2
brain2 = b2[,brain2_info$V3] #fully populated, as expected

#add column labels to b1_brain1 and b2_brain2
colnames(brain1) = brain1_info$V2
colnames(brain2) = brain2_info$V2

#save data
write.csv(brain1,"brain1_raw.csv", row.names = FALSE)
write.csv(brain2,"brain2_raw.csv", row.names = FALSE)
