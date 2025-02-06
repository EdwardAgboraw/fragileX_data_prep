#For both brains, generate barcode count histograms
#use histograms to pick a good barcode count threshold
library(ggplot2)

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

#start with raw data, see what happens, and then see what spike-in norm changes

#Brain2: Combine X.dACA and caudal.ACA into one column

b2_dACA = b2_norm$X.dACA
b2_caudalACA = b2_norm$caudal.ACA

b2_norm$caudal.ACA = NULL
b2_norm$X.dACA = NULL

b2_new_dACA = b2_dACA + b2_caudalACA
b2_norm$X.dACA = b2_new_dACA

#Histograms: Distribution of bar code counts for each region
#Test 1: BC histogram for one region, in one brain
b1_ORB = b1_norm$ORB
#ggplot(data = b1_norm, aes(x = VIS)) +
  #geom_histogram(binwidth = 10)

b1_ORB_nonZero = b1

b1_ORB_log = log10(b1_ORB)

hist(b1_ORB_log, breaks = 5)

#Manually Bin Data
#https://www.statology.org/data-binning-in-r/

#breaks should range from 0 to len(data), with intervals of 10

# Define bin breaks from 0 to max(vec) in steps of 10
breaks <- seq(0, max(b1_ORB+10), by = 10)

# Bin the values
b1_orb_bin <- cut(b1_ORB, breaks, include.lowest = TRUE, right = FALSE)

#b1_orb_bin = cut(b1_ORB, breaks=50)
b1_orb_hist_data = as.data.frame(table(b1_orb_bin))

b1_orb_hist_data$Frequency = log2(b1_orb_hist_data$Freq + 1)

#turn into barchart
b1_orb_chart = ggplot(b1_orb_hist_data, aes(x=b1_orb_bin, y=Frequency)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(b1_orb_chart)

#test function
not_targets_b1 = c("PL","RSC","negative.control..contalateral.olf..Bulb")
targets_b1 = colnames(b1)[!(colnames(b1) %in% not_targets_b1)]

b1_hists = make_barcode_count_histograms(b1_norm, targets_b1)

