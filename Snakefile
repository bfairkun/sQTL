include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "../../output/ChimpEgenes.eigenMT.txt.gz"

##### Modules #####
include: "rules/RNASeqMapping.smk"
include: "rules/eqtl_calling.smk"
include: "rules/sqtl_calling.smk"

