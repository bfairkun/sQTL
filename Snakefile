include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "sQTL_mapping/MatrixEQTL/Results/images/PCsVsSQTLs.pdf",
        # A plot of number of sQTLs (FDR correction based on nominal
        # snp:phenotype pair P-values). Useful to pick number of PCs to use for
        # more computationally intensive permutations.

        "sQTL_mapping/Chimp.eigenMT.txt.gz",
        # For the ideal number of PCs specified in the config file, this
        # contains FDR corrected qvalues for each phenotype.

        "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/PermutationsCombined.txt"
        # Unfortunately eigenMT doesn't have a method for getting qvalues at
        # the levels of grouped phenotypes which might make more sense for
        # splicing. So this is a table of permutated qvalues to use for a
        # permutation test. (I didn't provide a script to actualy do the
        # permutation test, just scripts to get this file which is a matrix
        # containing a list of permutated pvalues for each phenotype)  More
        # description is in the config file... You will have to edit some
        # scripts to your need if you want to group phenotypes. All of this
        # extra work is already implemented in FastQTL, which you could give a
        # try for comparison... But I have stuck to MatrixEQTL because it
        # provides a way to properly control for the relatedness between some
        # of the chimp individuals.

##### Modules #####
include: "rules/RNASeqMapping.smk"
include: "rules/eqtl_calling.smk"
include: "rules/sqtl_calling.smk"

