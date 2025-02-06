library(R.matlab)
library(tidyverse)
library(mltools)

#load data
data <- readMat('M283BarcodeMatrix.mat')

#get raw barcode matrix
bc = data$barcodematrix

#format bc

#label columns
num_col = 1:41
colnames(bc) = num_col

#load in sorted column values
colVals = read.csv("fragileX_sorted_column_labels_withBrain.csv", header = FALSE)

brain1_info = colVals[colVals$V1 == 1,]
brain2_info = colVals[colVals$V1 == 2,]

#split raw barcode matrix into individual brains
#for some reason rows are dr
brain1 = bc[,brain1_info$V3]
brain2 = bc[,brain2_info$V3]

#clean up column labels
brain1_colnames_tidy <- make.names(brain1_info$V2, unique=TRUE)
brain2_colnames_tidy <- make.names(brain2_info$V2, unique=TRUE)

#add column labels
colnames(brain1) = brain1_colnames_tidy
colnames(brain2) = brain2_colnames_tidy

write.csv(brain1, "brain1_column_corrected.csv", row.names = FALSE)
write.csv(brain2, "brain2_column_corrected.csv", row.names = FALSE)
