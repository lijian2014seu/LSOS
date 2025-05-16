### description:
# Determine consistency of probe-level changes between two groups:
#  - For each sample well, compare values in group2 vs group1
#  - Classify each feature as “Greater”, “Smaller”, or “Inconsistent”
#  - Returns a vector of classification labels for all rows in `df`

consist_judge_fun = function(group1, group2, df, sample_sheet) {
    samples <- unique(sample_sheet$Sample_Well)
    idx1 <- sapply(samples, function(s) which(sample_sheet$Sample_Well == s & sample_sheet$Sample_Group == group1))
    idx2 <- sapply(samples, function(s) which(sample_sheet$Sample_Well == s & sample_sheet$Sample_Group == group2))
    consist_judge_fun_comp1_to_comp2 <- apply(
    df,
    1,
    function(x) {
            b1 <- x[idx1];  b2 <- x[idx2]
            if (all(b2 > b1))       "Greater"
            else if (all(b2 < b1))  "Smaller"
            else                    "Inconsistent"
        }
    )
    consist_judge_fun_comp1_to_comp2
}