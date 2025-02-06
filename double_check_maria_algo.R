library(R.matlab)
library(tidyverse)
library(mltools)

#load brain data
brain1 = read.csv("brain1_column_corrected.csv")
#brain2 = read.csv("brain2_column_corrected.csv")

#load spikes
spike_in_counts = readRDS("fragileX_spikes.rds")

#specify injection sites
injection_sites = c("PL","RSC")

negative_control = "negative.control..contalateral.olf..Bulb"

data = brain1

data = as.data.frame(data)

target_sites = colnames(data)[!(colnames(data) %in% injection_sites)]
#target_sites = target_sites[!target_sites %in% negative_controls]

#format spike info
spike_info = data.frame(matrix(ncol = length(colnames(data)), nrow = 0))
colnames(spike_info) = colnames(data)
spike_info[1,] = spike_in_counts

spikes_targets = as.numeric(spike_info[,target_sites][1,])
spikes_inj = as.numeric(spike_info[,injection_sites][1,])

spikes_target_factor = max(spikes_targets)
spikes_inj_factor = max(spikes_inj)

#Apply Target and Inj Spike-In Normalization

data_targs = data[,target_sites]

for(i in 1:length(spikes_targets)){
  data_targs[,i] <- (data_targs[,i]/spikes_targets[i])*spikes_target_factor
}

data_inj = data[,injection_sites]

for(i in 1:length(spikes_inj)){
  data_inj[,i] <- (data_inj[,i]/spikes_inj[i])*spikes_inj_factor
}

data_norm = cbind(data_targs, data_inj)

# Filter 1: must have injection site count above 0
data_norm$boolean <- apply(data_norm[,injection_sites],1, function(x) all(x==0))
#returns TRUE if all checked regions have a count equal to 0

data_norm <- data_norm %>% filter(boolean == FALSE) %>% select(-boolean)
#only keep rows where at least one injection site count does not equal zero
#remove boolean column

# Filter 2: must have no target site count greater than largest inj site count

data_norm_targs = data_norm[,target_sites]
data_norm_inj = data_norm[,injection_sites]

data_norm_targs$max <- apply(data_norm_targs, 1, max, na.rm=TRUE)
data_norm_inj$max = apply(data_norm_inj, 1, max, na.rm=TRUE)

data_norm$targ_max = data_norm_targs$max
data_norm$inj_max = data_norm_inj$max

data_norm = data_norm %>% filter(inj_max > targ_max)

# Filter 3: must have count greater than 5 in at least 1 target site

data_norm = data_norm %>% filter(targ_max > 5)

# Filter 4: must have no count in the negative control

data_norm$no_nc <- data_norm[,negative_control] == 0
#returns TRUE if negative control count is 0

data_norm <- data_norm %>% filter(no_nc == TRUE) %>% select(-no_nc)
#only keep rows where at least one injection site count does not equal zero
#remove no_nc column

#create final dataset

data_norm$inj_max = NULL
data_norm$targ_max = NULL





