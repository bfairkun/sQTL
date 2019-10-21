include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "sQTL_mapping/MatrixEQTL/Results/images/PCsVsSQTLs.pdf",
        # A plot of number of sQTLs (FDR correction based on nominal
        # snp:phenotype pair P-values). Useful to pick number of PCs to use for
        # more computationally intensive permutations.

        "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/PermutationsCombined.txt",
        # this is a table of minimum pvalues for each intron from permuted data
        # to use for a permutation test to get intron level, cluster level, or
        # gene level P-values.  Each row is a permutation, and each column is
        # an intron.  More description is in the config file... As I wrote this
        # pipeline, the permutation test gets you cluster level P-values. If
        # you want gene level P-values, or intron level Pvalues, you will have
        # to edit the Rscript to read in this table, (group phenotypes and
        # retrieve minimum P- within group for each permutation) and compare
        # that to the minimum P from real data. All of this extra work is
        # already implemented in FastQTL, which you could give a try for
        # comparison since it is easy...  But I have stuck to MatrixEQTL
        # because it provides a way to properly control for the relatedness
        # between some of the chimp individuals.

        "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/BestPvalueNonPermuted.txt",
        "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/PermutationTestResults.txt"

##### Modules #####

include: "rules/RNASeqMapping.smk"
# contains a bunch of rules for RNA-seq mapping to bam files. You can use the
# bams I already mapped as specified in the config file.
include: "rules/eqtl_calling.smk"
# contains a bunch of rules for pre-processing genotype data for QTL calling
# with MatrixEQTL
include: "rules/sqtl_calling.smk"
# contains a bunch of rules that use Rscripts that use MatrixEQTL library.

