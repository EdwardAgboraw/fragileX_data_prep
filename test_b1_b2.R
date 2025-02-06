library(R.matlab)
library(tidyverse)

#load data
data <- readMat('M283BarcodeMatrix.mat')

#load b1norm and b2norm
b1 = data$B1norm
b2 = data$B2norm

#specifically label columns
num_col = 1:41
colnames(b1) = num_col
colnames(b2) = num_col

#load in sorted column values
colVals = read.csv("fragileX_sorted_column_labels_withBrain.csv", header = FALSE)

brain1_info = colVals[colVals$V1 == 1,]
brain2_info = colVals[colVals$V1 == 2,]

#test b1
b1_brain1 = b1[,brain1_info$V3] #fully populated, as expected
b1_brain2 = b1[,brain2_info$V3] #empty, as expected

#test b2
b2_brain1 = b2[,brain1_info$V3] #empty, as expected
b2_brain2 = b2[,brain2_info$V3] #fully populated, as expected.

#add column labels to b1_brain1 and b2_brain2
colnames(b1_brain1) = brain1_info$V2
colnames(b2_brain2) = brain2_info$V2

#save data
write.csv(b1_brain1,"brain1_facility_processed.csv", row.names = FALSE)
write.csv(b2_brain2,"brain2_facility_processed.csv", row.names = FALSE)



