library(R.matlab)
library(tidyverse)
library(mltools)

#load raw data files
#(replace with spike-in normalized files after running in-house preprocessing)
#brain1 = read.csv("brain1_raw.csv")
#brain2 = read.csv("brain2_raw.csv")

brain1 = read.csv("brain1_filtered.csv")
brain2 = read.csv("brain2_filtered.csv")

#Brain 1 Injection Sites = PL, RSC
#Brain 2 Inection Sites = X.dACA, PL, RSC, caudal.ACA

#Simplest Algorithm
#Extract injection sites from Brain, and then take the max val for each neuron.
#Binarize the remaining matrix, and the calculate the row sum.
#Row sum of binarized data is automatically # target sites.
#Create new Dataframe, injection site sum + number of target sites
#use to create graph

#Extract Injection Sites

injectionSites_brain1 = c("PL", "RSC")
brain1_inj = brain1[,injectionSites_brain1]

injectionSites_brain2 = c("X.dACA", "PL", "RSC", "caudal.ACA")
brain2_inj = brain2[,injectionSites_brain2]

#get max vals ("true" injection site count per neuron)
brain1_inj$max <- apply(brain1_inj, 1, max, na.rm=TRUE)
brain2_inj$max <- apply(brain2_inj, 1, max, na.rm=TRUE)

#Extract Target Region Data
# no injection sites AND no negative control
brain1_not_targets = c(injectionSites_brain1,"negative.control..contalateral.olf..Bulb")
brain2_not_targets = c(injectionSites_brain2,"negative.control..contalateral.olf..Bulb")

brain1_targets <- brain1[,!names(brain1) %in% brain1_not_targets]
brain2_targets = brain2[,!names(brain2) %in% brain2_not_targets]

#binarize data
brain1_targets_bin = as.data.frame(apply(brain1_targets[,-ncol(brain1_targets) -1],
                                         2, function(x) {ifelse(x>0, 1, 0)}))
brain2_targets_bin = as.data.frame(apply(brain2_targets[,-ncol(brain2_targets)-1],
                                         2, function(x) {ifelse(x>0, 1, 0)}))

#calculate rowsums
brain1_targets_bin$num_targets = rowSums(brain1_targets_bin)
brain2_targets_bin$num_targets = rowSums(brain2_targets_bin)

#Create final dataframe for graph
graph_table_brain1 = data.frame(brain1_inj$max, brain1_targets_bin$num_targets)
colnames(graph_table_brain1) = c("injection_site_count", "num.targets")

graph_table_brain2 = data.frame(brain2_inj$max, brain2_targets_bin$num_targets)
colnames(graph_table_brain2) = c("injection_site_count", "num.targets")

#Add Brain information
graph_table_brain1$brain = "one"
graph_table_brain2$brain = "two"

#combine
graph_table = rbind(graph_table_brain1, graph_table_brain2)

#Replicate Maria's Figure

#generate bins
bins_seq <- seq(1, max(graph_table$injection_site_count), by = 20)

graph_table$injection_binned <- bin_data(graph_table$injection_site_count,
                                         bins = bins_seq, boundaryType = "lcro]")

# prepare the dataset for the figure
graph_figure <- graph_table %>% group_by(injection_binned, brain, num.targets) %>% summarise(count = n()) %>%
  pivot_wider(names_from = num.targets, values_from = count) %>% as.data.frame()

graph_figure <- graph_figure %>% rename("one" = "1", "two" = "2", "three" = "3", "four" = "4",
                                          "five" = "5", "six" = "6", "seven" = "7", "eight" = "8", "nine" = "9",
                                          "ten" = "10", "eleven" = "11", "twelve" = "12")

# replace NAs with 0
graph_figure[is.na(graph_figure)] <- 0

# Convert the counts to percentages (normalised by the barcode counts intervals
#e.g., divide each value in an interval [1,51) in a given brain by the sum of the values in that row;
#>>> 62% of barcodes in [1,21) interval in brain 1 project to one target area, 21% to two areas, etc.)
graph_figure_percent <- cbind(injection_binned = graph_figure$injection_binned,
                               brain = graph_figure$brain,
                               as.data.frame(t(apply(graph_figure[,-c(1,2)],
                                                     1, function(row) round(row / sum(row), 5)))))

graph_figure_percent$brain <- as.factor(graph_figure_percent$brain)

#create long version
graph_figure_percent_long <- graph_figure_percent %>% pivot_longer(names_to = "num.targets", values_to = "percentage", cols = -c(injection_binned, brain))
graph_figure_percent_long$num.targets <- factor(graph_figure_percent_long$num.targets,
                                               levels = c("one", "two", "three", "four",
                                                          "five", "six", "seven", "eight",
                                                          "nine", "ten", "eleven",
                                                          "twelve"))

saveRDS(graph_figure_percent_long, "fragileX_inj_site_threshold_data.rds")

#make graphs
targets_fig <- ggplot(graph_figure_percent_long,
                      aes(x = injection_binned, y = percentage,
                          color = brain, group = brain)) +
  geom_line() +
  geom_point(alpha = 0.5) +
  facet_wrap(~num.targets) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title=element_text(size=14,face="bold")) +
  #theme(axis.text=element_text(size=12))
  scale_x_discrete(breaks = unique(graph_figure_percent_long$injection_binned)[seq(1, length(unique(graph_figure_percent_long$injection_binned)), 10)]) +
  xlab("Injection counts (binned every 20 counts)") +
  ylab("Proportion") +
  scale_color_discrete(breaks = c("one", "two"))

targets_fig
#works for raw data, simply need to functionalize and run for processed data
#works for filtered data

#make for only one -> three targets

graph_figure_percent_long_small = graph_figure_percent_long[graph_figure_percent_long$num.targets %in% c("one","two","three"), ]

targets_fig_small <- ggplot(graph_figure_percent_long_small,
                      aes(x = injection_binned, y = percentage,
                          color = brain, group = brain)) +
  geom_line() +
  geom_point(alpha = 0.5) +
  facet_wrap(~num.targets) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title=element_text(size=14,face="bold")) +
  #theme(axis.text=element_text(size=12))
  scale_x_discrete(breaks = unique(graph_figure_percent_long_small$injection_binned)[seq(1, length(unique(graph_figure_percent_long_small$injection_binned)), 10)]) +
  xlab("Injection counts (binned every 20 counts)") +
  ylab("Proportion") +
  scale_color_discrete(breaks = c("one", "two"))

targets_fig_small
