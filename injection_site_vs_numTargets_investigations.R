library(ggplot2)
source("mapseq_preprocess_functions.R")

#load data
plot_data = readRDS("fragileX_inj_site_threshold_data.rds")

#directly calculate Elbow via derivatives?

brain1_data = plot_data[plot_data$brain == "one",]
brain2_data = plot_data[plot_data$brain == "two",]

#one target graphs
b1_one = inj_count_v_target_num(brain1_data, "one")
b2_one = inj_count_v_target_num(brain2_data, "one")
print(b1_one)
print(b2_one)

#two target graphs
b1_two = inj_count_v_target_num(brain1_data, "two")
b2_two = inj_count_v_target_num(brain2_data, "two")
print(b1_two)
print(b2_two)

#three target graphs
b1_three = inj_count_v_target_num(brain1_data, "three")
b2_three = inj_count_v_target_num(brain2_data, "three")
print(b1_three)
print(b2_three)

#smooth regression lines on charts indicate the cutoff should be at bin 300
#
test = brain1_data[brain1_data$num.targets == "one",]
test_b2 = brain2_data[brain2_data$num.targets == "one",]

#injection site threshold = 800

#load filtered data
brain1 = read.csv("brain1_filtered.csv")
brain2 = read.csv("brain2_filtered.csv")

#apply threshold filter to the calculate the maximum injection
inj_sites_brain1 = c("PL","RSC")
inj_sites_brain2 = c("X.dACA","PL","RSC","caudal.ACA")

#Brain 1
brain1_inj = brain1[,inj_sites_brain1]
brain1_inj$max = apply(brain1_inj, 1, max, na.rm=TRUE)

brain1$inj_max = brain1_inj$max
brain1 = brain1 %>% filter(inj_max > 800)
#6163 final neurons

#Brain 2
brain2_inj = brain2[,inj_sites_brain2]
brain2_inj$max = apply(brain2_inj, 1, max, na.rm=TRUE)

brain2$inj_max = brain2_inj$max
brain2 = brain2 %>% filter(inj_max > 800)
#7664 final neurons

#filter out unwanted information
#remove injection sites
b1_target_sites = colnames(brain1)[!(colnames(brain1) %in% inj_sites_brain1)]
brain1 = brain1[,b1_target_sites]
#remove negative control
brain1$negative.control..contalateral.olf..Bulb = NULL
#remove inj max column
brain1$inj_max = NULL
#15 final columns (matches sample sheet)

#remove injection sites
b2_target_sites = colnames(brain2)[!(colnames(brain2) %in% inj_sites_brain2)]
brain2 = brain2[,b2_target_sites]
#remove negative control
brain2$negative.control..contalateral.olf..Bulb = NULL
#remove inj max column
brain2$inj_max = NULL
#14 final columns (matches sample sheet)

#save data
write.csv(brain1, "brain1_filtered_inj_thresh_800.csv", row.names = FALSE)
write.csv(brain2,"brain2_filtered_inj_thresh_800.csv", row.names = FALSE)

