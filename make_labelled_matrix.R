library(R.matlab)
library(tidyverse)

data <- readMat('M283BarcodeMatrix.mat')

#raw, unfiltered barcode count matrix (both mice)
bcmat = data$barcodematrix

#specifically label columns by number
num_col = 1:41
colnames(bcmat) = num_col

#load in sorted column values
colVals = read.csv("fragileX_sorted_column_labels_withBrain.csv", header = FALSE)

brain1_labels = colVals[colVals$V1 == 1,]
b1_cols = brain1_labels$V3

brain2_labels = colVals[colVals$V1 == 2,]
b2_cols = brain2_labels$V3

#split by brain
brain1_raw = bcmat[,b1_cols]
colnames(brain1_raw) = brain1_labels$V2

brain2_raw = bcmat[,b2_cols]
colnames(brain2_raw) = brain2_labels$V2

#repeat for the filtered and normalized datasets
#B1norm dataset
b1norm = data$B1norm
#brain1
brain1_b1norm = b1norm[,b1_cols]
colnames(brain1_b1norm) = brain1_labels$V2
#brain2
brain2_b1norm = b1norm[,b2_cols]
colnames(brain2_b1norm) = brain2_labels$V2

#B2norm dataset
b2norm = data$B2norm
#brain1
brain1_b2norm = b2norm[,b1_cols]
colnames(brain1_b2norm) = brain1_labels$V2
#brain2
brain2_b2norm = b2norm[,b2_cols]
colnames(brain2_b2norm) = brain2_labels$V2













