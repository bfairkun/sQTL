include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "sQTL_mapping/MatrixEQTL/Results/images/PCsVsSQTLs.pdf",
        # A plot of number of sQTLs (FDR correction based on nominal
        # snp:phenotype pair P-values). Useful to pick number of PCs to use for
        # more computationally intensive permutations.

        "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/PermutationsCombined.txt"
        # this is a table of permutated pvalues to use for a permutation test
        # to get intron level, cluster level, or gene level P-values. (I didn't
        # provide a script to actualy do the permutation test, just scripts to
        # get this file which is a matrix containing a list of the smalest
        # pvalue for each intron after shuffling sample labels. This serves as
        # a null distribution of minimum P-values from which to compute
        # intron/cluster/gene level P-values). Each row in this file is a
        # permutation, and each column is an intron.  More description is in
        # the config file... You will have to create your own script to read in
        # this table, (group phenotypes and retrieve minimum P- within group
        # for each permutation) and compare that to the minimum P from real
        # data. All of this extra work is already implemented in FastQTL, which
        # you could give a try for comparison since it is easy...  But I have
        # stuck to MatrixEQTL because it provides a way to properly control for
        # the relatedness between some of the chimp individuals.

##### Modules #####
include: "rules/RNASeqMapping.smk"
include: "rules/eqtl_calling.smk"
include: "rules/sqtl_calling.smk"

