
library(tidyverse)
library(mltools)
library(ggplot2)

process_MAPseq <- function(data, spike_in_counts, injection_sites, negative_control){

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

  #save column name of max injection site
  inj_site = colnames(data_norm_inj)[apply(data_norm_inj,1,which.max)]
  #save max value in new column
  data_norm_inj$max = apply(data_norm_inj, 1, max, na.rm=TRUE)

  data_norm$targ_max = data_norm_targs$max
  data_norm$inj_max = data_norm_inj$max

  #add injection site information
  data_norm$injection_site = inj_site

  data_norm = data_norm %>% filter(inj_max > targ_max)

  # Filter 3: must have count greater than 5 in at least 1 target site
  data_norm = data_norm %>% filter(targ_max > 5)

  # Filter 4: must have no count in the negative control
  data_norm$no_nc <- data_norm[,negative_control] == 0
  #returns TRUE if negative control count is 0

  data_norm <- data_norm %>% filter(no_nc == TRUE) %>% select(-no_nc)
  #only keep rows where the negative_control column value is 0
  #remove no_nc column

  #create final dataset

  data_norm$inj_max = NULL
  data_norm$targ_max = NULL

  return(data_norm)

}

spike_in_normalize <- function(data, spike_in_counts, injection_sites){

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

  return(data_norm)

}

facility_filter <- function(data, injection_sites){

  #replicates facility filtering scheme

  data = as.data.frame(data)

  target_sites = colnames(data)[!(colnames(data) %in% injection_sites)]

  #injection site filter

  data_inj = data[,injection_sites]
  data_inj$max = apply(data_inj, 1, max, na.rm=TRUE)

  data$inj_max = data_inj$max

  data = data %>% filter(inj_max > 30)

  #projection area filter

  data_targ = data[,target_sites]
  data_targ$max = apply(data_targ, 1, max, na.rm=TRUE)

  data$targ_max = data_targ$max

  data = data %>% filter(targ_max > 5)

  #create final dataset

  data$inj_max = NULL
  data$targ_max = NULL

  return(data)

}

inj_count_v_target_num <- function(data, targ_num){

  #subset to a specific number of targets
  targ_data = data[data$num.targets == targ_num,]

  #replace count bins with simple continous variable
  targ_data$bin_label = 1:nrow(targ_data)

  #plot graph
  g1 <- ggplot(targ_data, aes(x = bin_label, y = percentage))
  #g1 + geom_point()
  # graph = g1 + geom_smooth(method = "loess") + ggtitle(targ_num)
  graph = g1 + geom_line() + ggtitle(targ_num)

  return(graph)

}

prep_for_graph <- function(brain){

  #remove injection site columns
  brain$injection_site = NULL

  #remove negative control
  brain$negative.control..contalateral.olf..Bulb = NULL

  return(brain)
}

compare_brains <- function(brain1, brain2, title, comp1, comp2){

  shared_targets = colnames(brain1)[(colnames(brain1) %in% colnames(brain2))]

  #binarize data
  brain1_shared_bin = ifelse(brain1 > 0,1,0)
  brain2_shared_bin = ifelse(brain2 > 0,1,0)

  #calculate per-region projection proportions
  brain1_sums = colSums(brain1_shared_bin)
  brain2_sums = colSums(brain2_shared_bin)
  #sums are in the correct order

  brain1_props = c()
  for (x in brain1_sums){

    prop = x/nrow(brain1)

    brain1_props = c(brain1_props, prop)

  }

  brain2_props = c()
  for (x in brain2_sums){

    prop = x/nrow(brain2)

    brain2_props = c(brain2_props, prop)

  }

  #proportions are in the correct order

  # Create wide format table
  prop_frame = tibble(
    Region = shared_targets,
    {{comp1}} := brain1_props,
    {{comp2}} := brain2_props
  )

  print(head(prop_frame))

  #Pivot to long format
  prop_frame <- prop_frame %>%
    pivot_longer(cols = c(comp1, comp2), names_to = "Brain", values_to = "Proportion")

  # Define dodge position with a wider gap
  dodge <- position_dodge(width = 0.8)  # Controls spacing between pairs

  # Create correct bar chart
  chart = ggplot(prop_frame, aes(x = Region, y = Proportion, fill = Brain)) +
    geom_bar(stat = "identity", position = dodge, width = 0.84) +
    labs(title = title, x = "Brain Region", y = "Proportion") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+
    #scale_fill_manual(values = c("#A020F0","#00FF00"))

  print(chart)
}

make_barcode_count_histograms <- function(data, target_sites){


  histograms = list()

  for (region in target_sites){

    counts = pull(data, region)

    # Define bin breaks from 0 to max(vec) in steps of 10
    breaks <- seq(0, max(counts+5), by = 5)

    # Bin the values
    counts_bin <- cut(counts, breaks, include.lowest = TRUE, right = FALSE)

    hist_data = as.data.frame(table(counts_bin))

    hist_data$Frequency = log2(hist_data$Freq + 1)

    #turn into barchart
    chart = ggplot(hist_data, aes(x=counts_bin, y=Frequency)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(x = as.character(region), y = "Log2(Frequency)")

    #add chart to list
    histograms[[region]] = chart


  }

  #return list of histograms
  return(histograms)
}

#injection_site_vs_target_num <- function(data, injection_sites){
 #
  #target_sites = colnames(data)[!(colnames(data) %in% injection_sites)]
  #
  #data = as.data.frame(data)
#
 # data_inj = data[,injection_sites]
  #
  #data_inj$max <- apply(data_inj, 1, max, na.rm=TRUE)
  #
  #data_targets <- data[,target_sites]
  #
  #data_targets_bin = as.data.frame(apply(data_targets[,-ncol(data_targets) -1],
   #                                        2, function(x) {ifelse(x>0, 1, 0)}))
  #
  #data_targets_bin$num_targets = rowSums(data_targets_bin)
  #
  #graph_table = data.frame(data_inj$max, data_targets_bin$num_targets)
  #colnames(graph_table) = c("injection_site_count", "num.targets")
#}

