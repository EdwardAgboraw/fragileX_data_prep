source("mapseq_preprocess_functions.R")
#load data
brain1 = read.csv("brain1_filt_inj_thresh_800_test.csv")
brain2 = read.csv("brain2_filt_inj_thresh_800_test.csv")

#get shared targets
shared_targets = colnames(brain1)[(colnames(brain1) %in% colnames(brain2))]

#subset
brain1_shared = brain1[,shared_targets]
brain2_shared = brain2[,shared_targets]

write.csv(brain1_shared,"brain1_chart_data_fx.csv")
write.csv(brain2_shared,"brain2_chart_data_fx.csv")

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

#Brain1 (WT) vs Brain 2 (KO), all neurons
compare_brains(brain1_shared, brain2_shared, "Wildtype vs xKO, all (Proportions)", "Wildtype", "xKO")

compare_brains_counts(brain1_shared, brain2_shared, "Wildtype vs xKO, all (Total Barcode Counts)", "Wildtype", "xKO")

#Brain 1 vs Brain 2, PL injection site
compare_brains(brain1_shared_PL, brain2_shared_PL, "Wildtype vs xKO, PL (Proportions)", "Wildtype", "xKO")

compare_brains_counts(brain1_shared_PL, brain2_shared_PL, "Wildtype vs xKO, PL (Total Barcode Counts)", "Wildtype", "xKO")

#Brain1 vs Brain2, RSC injection site
compare_brains(brain1_shared_RSC, brain2_shared_RSC, "Wildtype vs xKO, RSC (Proportions)", "Wildtype", "xKO")

compare_brains_counts(brain1_shared_RSC, brain2_shared_RSC, "Wildtype vs xKO, RSC (Total Barcode Counts)", "Wildtype", "xKO")


#within brain comparisions
#Brain 1 RSC vs PL
compare_brains(brain1_shared_RSC, brain1_shared_PL, "Wildtype, RSC vs PL (Proportions)", "RSC", "PL")
compare_brains_counts(brain1_shared_RSC, brain1_shared_PL, "Wildtype, RSC vs PL (Total Barcode Counts)", "RSC", "PL")

#Brain 2 RSC vs PL
compare_brains(brain2_shared_RSC, brain2_shared_PL, "xKO, RSC vs PL (Proportions)", "RSC", "PL")
compare_brains_counts(brain2_shared_RSC, brain2_shared_PL, "xKO, RSC vs PL (Total Barcode Counts)", "RSC", "PL")

#Acb = Acumbens
#CP = caudate putamen

#compare raw counts









