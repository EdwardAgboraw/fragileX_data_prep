#exactly how many neurons in Brain 1 and Brain 2 do we lose when apply the
#"must have bc count above 5 in at least 1 target area" filter?

library(R.matlab)
library(tidyverse)
library(mltools)
library(dplyr)
source("mapseq_preprocess_functions.R")

#load brain data
b1 = read.csv("brain1_column_corrected.csv")
b2 = read.csv("brain2_column_corrected.csv")

#load spike-in-count values
spikes = readRDS("fragileX_spikes.rds")

#specify injection sites
inj_sites_b1 = c("PL","RSC")
inj_sites_b2 = c("X.dACA","PL","RSC","caudal.ACA")

#spike in normalize
b1_norm = spike_in_normalize(b1, spikes, inj_sites_b1)
b2_norm = spike_in_normalize(b2, spikes, inj_sites_b2)
#pre-filtering, both = 1196932 neurons

#apply filtering steps
negative_control = "negative.control..contalateral.olf..Bulb"

data_norm = b1_norm
injection_sites = inj_sites_b1
#data_norm = b2_norm
#injection_sites = inj_sites_b2

target_sites = colnames(data_norm)[!(colnames(data_norm) %in% injection_sites)]

# Filter 1: must have injection site count above 0
data_norm$boolean <- apply(data_norm[,injection_sites],1, function(x) all(x==0))
#returns TRUE if all checked regions have a count equal to 0

data_norm <- data_norm %>% filter(boolean == FALSE) %>% select(-boolean)
#brain1: 539109 neurons
#brain2: 400739

# Filter 2: must have no target site count greater than largest inj site count

data_norm_targs = data_norm[,target_sites]
data_norm_inj = data_norm[,injection_sites]

data_norm_targs$max <- apply(data_norm_targs, 1, max, na.rm=TRUE)

#save column name of max injection site
injection_sites = colnames(data_norm_inj)[apply(data_norm_inj,1,which.max)]
#save max value in new column
data_norm_inj$max = apply(data_norm_inj, 1, max, na.rm=TRUE)

data_norm$targ_max = data_norm_targs$max
data_norm$inj_max = data_norm_inj$max

#add injection site information
data_norm$injection_site = injection_sites

data_norm = data_norm %>% filter(inj_max > targ_max)
#brain1: 505912
#brain2: 384098

# Filter 3: must have count greater than 5 in at least 1 target site
data_norm = data_norm %>% filter(targ_max > 5)
#brain1 = 50133 neurons. 488976 neurons lost
#brain2 = 37393 neurons. 346705 neurons lost

# Filter 4: must have no count in the negative control
data_norm$no_nc <- data_norm[,negative_control] == 0
#returns TRUE if negative control count is 0
data_norm <- data_norm %>% filter(no_nc == TRUE) %>% select(-no_nc)
#brain1 = 49239
#brain2 = 36692

#apply injection site threshold filtering (800)
data_norm <- data_norm %>%
  rowwise() %>%
  filter(get(injection_site) > 800) %>%
  ungroup()
#brain1 = 6163 neurons. 43076 neurons lost.
#brain2 = 7664 neurons. 29028 neurons lost.
