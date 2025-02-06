library(R.matlab)
library(tidyverse)
library(mltools)
library(dplyr)
source("mapseq_preprocess_functions.R")

#load data collection
data <- readMat('M283BarcodeMatrix.mat')

#load spike-in-count values
spikes = readRDS("fragileX_spikes.rds")

#load brain data
b1 = read.csv("brain1_column_corrected.csv")
b2 = read.csv("brain2_column_corrected.csv")

#the goal here is to standardize injection sites between brains.

#specify injection sites
inj_sites_b1 = c("PL","RSC")
inj_sites_b2 = c("X.dACA","PL","RSC","caudal.ACA")

#specify negative control
neg_con = "negative.control..contalateral.olf..Bulb"

#filter MAPseq data (Maria Method)
b1_filt = process_MAPseq(b1, spikes, inj_sites_b1, neg_con)
b2_filt = process_MAPseq(b2, spikes, inj_sites_b2, neg_con)

#Injection Site Threshold Filtering (800)
#Brain 1

b1_filt_thresh <- b1_filt %>%
  rowwise() %>%
  filter(get(injection_site) > 800) %>%
  ungroup()
#6163 final neurons

#Brain 2
b2_filt_thresh <- b2_filt %>%
  rowwise() %>%
  filter(get(injection_site) > 800) %>%
  ungroup()
#7664 final neurons

#Only keep neurons in brain 2 for which PL and RSC are the injection site

b2_filt_thresh = b2_filt_thresh[b2_filt_thresh$injection_site %in% c("PL","RSC"), ]

# Brain2: combine x.dACA and caudal ACA into one.
b2_dACA = b2_filt_thresh$X.dACA
b2_caudalACA = b2_filt_thresh$caudal.ACA

b2_new_dACA = b2_dACA + b2_caudalACA

b2_filt_thresh$caudal.ACA = NULL
b2_filt_thresh$X.dACA = NULL

b2_filt_thresh$X.dACA = b2_new_dACA

#now we want to generate proportion bar charts, but we don't want the injection site
#to contribute to projection proportions

#save files
write.csv(b1_filt_thresh, "brain1_filt_inj_thresh_800_test.csv", row.names = FALSE)
write.csv(b2_filt_thresh, "brain2_filt_inj_thresh_800_test.csv", row.names = FALSE)



