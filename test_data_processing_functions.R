library(R.matlab)
library(tidyverse)
library(mltools)
source("mapseq_preprocess_functions.R")

#load data collection
data <- readMat('M283BarcodeMatrix.mat')

#load spike-in-count values
spikes = readRDS("fragileX_spikes.rds")

#load brain data
brain1 = read.csv("brain1_column_corrected.csv")
brain2 = read.csv("brain2_column_corrected.csv")

#specify injection sites
inj_sites_brain1 = c("PL","RSC")
inj_sites_brain2 = c("X.dACA","PL","RSC","caudal.ACA")

#test facility protocol (works! replicates B1 and B2 objects in data)
brain1_filt = facility_filter(brain1, inj_sites_brain1)
brain2_filt = facility_filter(brain2, inj_sites_brain2)

#test maria's filtering protocol

neg_con = "negative.control..contalateral.olf..Bulb"
brain1_maria_filt = process_MAPseq(brain1, spikes, inj_sites_brain1, neg_con)
brain2_maria_filt = process_MAPseq(brain2, spikes, inj_sites_brain2, neg_con)

#save
write.csv(brain1_maria_filt, "brain1_filtered.csv", row.names = FALSE)
write.csv(brain2_maria_filt, "brain2_filtered.csv", row.names = FALSE)
