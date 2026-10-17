"""Thresholds registered before any audit result is interpreted.

CFD 0.2 is the activity floor used by several guide-design pipelines when they
separate "report" from "treat as a plausible cut" (Doench et al. 2016 scores a
perfect match as 1). Allele frequency 0.05 is the usual common-variant cutoff.
Sensitivity cutoffs are reported beside the primary cutoff; they are not used
to choose the headline.
"""

ACTIONABLE_CFD = 0.2
COMMON_AF = 0.05
MAX_MISMATCHES = 4
INITIAL_PAM = ("NGG",)
EXPANDED_PAM = ("NAG", "NGA")
SENSITIVITY_CFD = (0.1, 0.5)
SENSITIVITY_AF = (0.01, 0.10)

# Evo 2 calls are the expensive second opinion. The cap is a budget, not a
# scientific filter: every site that can cross the CFD threshold is still scored
# with CFD and gnomAD. Only the Evo 2 prior is capped, and the trace says so.
EVO2_BUDGET = 4

# HBB on GRCh38, 1-based inclusive, padded through the proximal promoter so a
# perfect match just outside the Ensembl gene span is still the on-target.
GENE_LOCI = {
    "HBB": {"chrom": "chr11", "start": 5225464, "end": 5227400},
}
