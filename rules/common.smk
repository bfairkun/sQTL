import pandas as pd
import os
from snakemake.utils import validate
from collections import defaultdict
from itertools import chain


###### Config file and sample sheets #####
configfile: "config.yaml"
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)


RNASeqFastqList = pd.read_table(config["RNASeqFileList"], squeeze=True, dtype=str).set_index(["sample"])
RNASeqSampleToFastq_dict = defaultdict(list)
RNASeqBasenameToFastq = dict()
with open(config["RNASeqFileList"]) as RNASeqFileList_fh:
    RNASeqFileList_fh.readline()
    for line in RNASeqFileList_fh:
        samplename, filepath = line.strip('\n').split('\t')
        RNASeqSampleToFastq_dict[samplename].append(filepath)
        RNASeqBasenameToFastq[os.path.basename(filepath)] = filepath

# contigs in reference genome
contigs = pd.read_table(config["ref"]["genome"] + ".fai",
                        header=None, usecols=[0], squeeze=True, dtype=str)
blacklist_contigs = open(config["sQTL_mapping"]["chromosome_blacklist"], 'r').read().split('\n')

##### Wildcard constraints #####
wildcard_constraints:
    sample="|".join(samples.index)

##### Helper functions #####

def getFastqFromBase(wildcards):
    return(RNASeqBasenameToFastq[wildcards.fastq_basename])

def get_RNASeq_merged_fastq(wildcards):
    """Get merged RNA fastq file for a given {sample} wildcard"""
    return config["temp_files_prefix"] + "RNASeqFastq/{sample}.fastq.gz".format(sample=wildcards.RNASeqSample)
