source("mapseq_preprocess_functions.R")
#load data
brain1 = read.csv("brain1_filt_inj_thresh_800_test.csv")
brain2 = read.csv("brain2_filt_inj_thresh_800_test.csv")

#get shared targets
shared_targets = colnames(brain1)[(colnames(brain1) %in% colnames(brain2))]

#subset
brain1_shared = brain1[,shared_targets]
brain2_shared = brain2[,shared_targets]

#iterate through/transform data.frame
#if column name matches injection_site value
#set column value to 0

for (i in seq_len(nrow(brain1_shared))) {
  inj_site <- brain1_shared$injection_site[i]  # Get the column name from X
  brain1_shared[i, inj_site] <- 0  # Set the corresponding column to 0
}

for (i in seq_len(nrow(brain2_shared))) {
  inj_site <- brain2_shared$injection_site[i]  # Get the column name from X
  brain2_shared[i, inj_site] <- 0  # Set the corresponding column to 0
}

#generate within-brain datasets for later comparision
brain1_shared_RSC = brain1_shared[brain1_shared$injection_site == "RSC",]
brain1_shared_PL = brain1_shared[brain1_shared$injection_site == "PL",]

brain2_shared_RSC = brain2_shared[brain2_shared$injection_site == "RSC",]
brain2_shared_PL = brain2_shared[brain2_shared$injection_site == "PL",]

#generate brain comparision graphs
brain1_comp = prep_for_graph(brain1_shared)
brain2_comp = prep_for_graph(brain2_shared)
compare_brains(brain1_comp, brain2_comp, "Wildtype vs xKO", "Wildtype", "xKO")

brain1_comp_RSC = prep_for_graph(brain1_shared_RSC)
brain1_comp_PL = prep_for_graph(brain1_shared_PL)
compare_brains(brain1_comp_RSC, brain1_comp_PL, "Wildtype, RSC vs PL", "RSC", "PL")


brain2_comp_RSC = prep_for_graph(brain2_shared_RSC)
brain2_comp_PL = prep_for_graph(brain2_shared_PL)
compare_brains(brain2_comp_RSC, brain2_comp_PL, "xKO, RSC vs PL", "RSC", "PL")

compare_brains(brain1_comp_PL, brain2_comp_PL, "Wildtype vs xKO, PL", "Wildtype", "xKO")
#compare_brains(brain1_comp_PL, brain2_comp_PL, "Brain1 vs Brain2, PL", "Brain1", "Brain2")

compare_brains(brain1_comp_RSC, brain2_comp_RSC, "Wildtype vs xKO, RSC", "Wildtype", "xKO")

#Acb = Acumbens
#CP = caudate putamen

