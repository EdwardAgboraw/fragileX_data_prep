library(R.matlab)
library(tidyverse)
library(mltools)

#load data
data <- readMat('M283BarcodeMatrix.mat')

#goal - generate a spike-in normalized data-set according to our in-house protocol
#compare to the normalized data sets provided by the center

#get raw barcode matrix
raw = data$barcodematrix # 49074212 completely unfiltered cells

#remove NAs
bc<- na.omit(raw) # no change in length

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

#specify brain1 and brain2 target and inj_sites

inj_sites_brain1 = c("PL","RSC") #will be user-provided arg in function
target_sites_brain1 = brain1_colnames_tidy[!(brain1_colnames_tidy %in% inj_sites_brain1)]

inj_sites_brain2 = c("X.dACA","PL","RSC","caudal.ACA")
target_sites_brain2 = brain2_colnames_tidy[!(brain2_colnames_tidy %in% inj_sites_brain2)]

#get spike-in counts

spikes = data$spikes
spike_in_counts = c()

for (i in 1:41){

  slice = spikes[,,i]
  countval = length(slice$counts2u)
  spike_in_counts = c(spike_in_counts,countval)

}

saveRDS(spike_in_counts, "fragileX_spikes.rds")

#get brain-specific spikes
spikes_brain1 = spike_in_counts[brain1_info$V3]
spikes_brain2 = spike_in_counts[brain2_info$V3]

#format spike info
brain1_spike_info = data.frame(matrix(ncol = length(brain1_colnames_tidy), nrow = 0))
colnames(brain1_spike_info) = brain1_colnames_tidy
brain1_spike_info[1,] = spikes_brain1

brain2_spike_info = data.frame(matrix(ncol = length(brain2_colnames_tidy), nrow = 0))
colnames(brain2_spike_info) = brain2_colnames_tidy
brain2_spike_info[1,] = spikes_brain2

#split into targets and injection site

#Brain1
spikes_brain1_targets = as.numeric(brain1_spike_info[,target_sites_brain1][1,])
spikes_brain1_inj = as.numeric(brain1_spike_info[,inj_sites_brain1][1,])

spikes_brain1_target_factor = max(spikes_brain1_targets)
spikes_brain1_inj_factor = max(spikes_brain1_inj)

#Brain2
spikes_brain2_targets = as.numeric(brain2_spike_info[,target_sites_brain2][1,])
spikes_brain2_inj = as.numeric(brain2_spike_info[,inj_sites_brain2][1,])

spikes_brain2_target_factor = max(spikes_brain2_targets)
spikes_brain2_inj_factor = max(spikes_brain2_inj)

#Apply Spike-In Normalization to Brain 1

brain1_norm = brain1
for(i in 1:16){
  brain1_norm[,i] <- (brain1[,i]/spikes_brain1[i])*spikes_brain1_target_factor
}

for(i in 17:18){
  brain1_norm[,i] <- (brain1[,i]/spikes_brain1[i])*spikes_brain1_inj_factor
}

#Apply Spike-In Normalization to Brain 2

brain2_norm = brain2
for(i in 1:15){
  brain2_norm[,i] <- (brain2[,i]/spikes_brain2[i])*spikes_brain2_target_factor
}

for(i in 16:19){
  brain2_norm[,i] <- (brain2[,i]/spikes_brain2[i])*spikes_brain2_inj_factor
}
#appears to be working as intended

#filter 1; must have inj count above 0

brain1_norm = as.data.frame(brain1_norm)
brain2_norm = as.data.frame(brain2_norm)

brain1f = brain1_norm %>% filter(PL > 0 | RSC > 0) # 539109 cells
brain2f = brain2_norm %>% filter(X.dACA > 0 | PL > 0 | RSC > 0 | caudal.ACA > 0)
#400739 cells

#filter 2; must have no target site count greater than largest inj site count

brain1f_targs = brain1f[,target_sites_brain1]
brain1f_inj = brain1f[,inj_sites_brain1]

brain1f_targs$max <- apply(brain1f_targs, 1, max, na.rm=TRUE)
brain1f_inj$max = apply(brain1f_inj, 1, max, na.rm=TRUE)

brain1f$targmax = brain1f_targs$max
brain1f$injmax = brain1f_inj$max

brain1f = brain1f %>% filter(injmax > targmax) # 504691 cells


brain2f_targs = brain2f[head(names(brain2f),-4)]
brain2f_inj = brain2f[tail(names(brain2f),4)]

#filter 3; must have no count in negative control

brain1f = brain1f %>% filter(negative.control..contalateral.olf..Bulb == 0)
#503638

#filter 4: must have count in target site greater than 5

brain1f$boolean <- apply(brain1f[,target_sites_brain1],1, function(x) all(x<5))
#returns TRUE if all checked regions have a count less than 5

brain1f_final <- brain1f %>% filter(boolean == FALSE) %>% select(-boolean)
#only keep neurons with at least one target site count above 5,
#then remove the boolean column
#49806 neurons, kinda close to the 40700 neurons reported by the facility
#and this is before injection vs threshold filtering

#save object, --> inj_count vs #target site thresholding




