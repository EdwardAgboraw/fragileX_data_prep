library(R.matlab)
library(tidyverse)

data <- readMat('M283BarcodeMatrix.mat')

#raw, unfiltered barcode count matrix (both mice)
bcmat = data$barcodematrix

#Brain One
#b1columns = c(1,35,2,3,4,5,6,7,8,9,36,10,11,12,13,14,15,16)
b1columns = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,35,36)
#Brain Two
#b2columns = c(37,38,17,18,19,20,21,22,23,24,39,25,26,27,28,29,30,31,40)
b2columns =  c(17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,37,38,39,40)

#Brain 1 and 2 don't project to identical regions

#B1 and B2 both contain the total 41 columns, and therefore contain information from each brain
#B1 = 40700 cells
#B2 = 27822 cells. Is B2 the filtered set?

#32 = target (negative control)
#34 = target (water control)
#41 = injection (water control)

#load in sorted column values
colVals = read.csv("fragileX_sorted_column_labels.csv", header = FALSE)
#index matches column names
b1_colVals = colVals[b1columns,]
b2_colVals = colVals[b2columns,]

brain1_raw <- bcmat[,b1columns]

colnames(bcmat) = labels

b1_test_norm = data$B1norm
b2_test_norm = data$B2norm
