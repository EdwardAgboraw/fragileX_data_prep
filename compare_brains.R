#load data
brain1 = read.csv("brain1_filtered_inj_thresh_800.csv")
brain2 = read.csv("brain2_filtered_inj_thresh_800.csv")

#get shared targets
shared_targets = colnames(brain1)[(colnames(brain1) %in% colnames(brain2))]

#subset
brain1_shared = brain1[,shared_targets]
brain2_shared = brain2[,shared_targets]

#force dfs to share the same column order
brain2_shared <- brain2_shared[, colnames(brain1_shared)]

#Basic Comparison - Proportion of neurons which project to each area

#binarize data
brain1_shared_bin = ifelse(brain1_shared > 0,1,0)
brain2_shared_bin = ifelse(brain2_shared > 0,1,0)
#column order still starts with ORB

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
  Brain1 = brain1_props,
  Brain2 = brain2_props
)

#Pivot to long format
prop_frame <- prop_frame %>%
  pivot_longer(cols = c(Brain1, Brain2), names_to = "Brain", values_to = "Proportion")

# Define dodge position with a wider gap
dodge <- position_dodge(width = 0.8)  # Controls spacing between pairs

# Create correct bar chart
ggplot(prop_frame, aes(x = Region, y = Proportion, fill = Brain)) +
  geom_bar(stat = "identity", position = dodge, width = 0.84) +
  labs(title = "Brain 1 vs Brain 2", x = "Brain Region", y = "Proportion") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))



